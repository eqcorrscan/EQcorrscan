"""
Correlation functions for multi-channel cross-correlation of seismic data.

Various routines used mostly for testing, including links to a compiled
routine using FFTW, a Numpy fft routine which uses bottleneck for normalisation
and a compiled time-domain routine. These have varying levels of efficiency,
both in terms of overall speed, and in memory usage.  The time-domain is the
most memory efficient but slowest routine (although fast for small cases of
less than a few hundred correlations), the Numpy routine is fast, but memory
inefficient due to a need to store large double-precision arrays for
normalisation.  The fftw compiled routine is faster and more memory efficient
than the Numpy routine.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import contextlib
import copy
import ctypes
import os
import logging
from multiprocessing import Pool as ProcessPool, cpu_count
from multiprocessing.pool import ThreadPool

import numpy as np
import math

from packaging import version
from timeit import default_timer

from obspy import UTCDateTime

from eqcorrscan.utils.libnames import _load_cdll
from eqcorrscan.utils import FMF_INSTALLED
from eqcorrscan.utils.pre_processing import _stream_quick_select


Logger = logging.getLogger(__name__)

# This is for building docs on readthedocs, which has an old version of
# scipy - without this, this module cannot be imported, which breaks the docs
# Instead we define a dummy function that returns what it is given.
READ_THE_DOCS = os.environ.get('READTHEDOCS', None) == 'True'
if not READ_THE_DOCS:
    from scipy.fftpack.helper import next_fast_len
else:
    def next_fast_len(a):
        return a

XCOR_FUNCS = {}  # cache of functions for doing cross correlations

# methods added to each xcorr func registered
# these implement the stream interface
XCORR_STREAM_METHODS = ('multithread', 'multiprocess', 'concurrent',
                        'stream_xcorr')

# these implement the array interface
XCOR_ARRAY_METHODS = ('array_xcorr')

# Gain shift for low-variance stabilization
MULTIPLIER = 1e8

# Minimum version for compatible correlations from Fast Matched Filter
MIN_FMF_VERSION = version.parse("1.4.0")


class CorrelationError(Exception):
    """ Error handling for correlation functions. """

    def __init__(self, value):
        """ Raise error.

        .. rubric:: Example

        >>> CorrelationError("Something happened")
        Something happened
        """
        self.value = value

    def __repr__(self):
        """ Print error value.

        .. rubric:: Example

        >>> print(CorrelationError("Error").__repr__())
        Error
        """
        return self.value

    def __str__(self):
        """ Print otherwise

        .. rubric:: Example

        >>> print(CorrelationError("Error"))
        Error
        """
        return self.value


# ------------------ Context manager for switching out default
class _Context:
    """ class for permanently or temporarily changing items in a dict """

    def __init__(self, cache, value_to_switch):
        """
        :type cache: dict
        :param cache: A dict to store values in
        :type value_to_switch: str
        :param value_to_switch:
            The key in cache to switch based on different contexts
        """
        self.cache = cache
        self.value_to_switch = value_to_switch
        self.previous_value = None

    def __call__(self, new_value, *args, **kwargs):
        """ # TODO change docs if this ever becomes general use
        Set a new value for the default xcorr function.

        This function can be called directly to permanently change the
        default normxcorr function or it may be used as a context manager
        to only modify it temporarily.

        :param new_value:
        :return:
        """
        self.previous_value = copy.deepcopy(
            self.cache.get(self.value_to_switch))
        self.cache[self.value_to_switch] = get_array_xcorr(new_value)
        return self

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.revert()

    def __repr__(self):
        """ this hides the fact _Context instance are returned after calls """
        name = self.cache[self.value_to_switch].__str__()
        if hasattr(self.cache[self.value_to_switch], '__name__'):
            name = self.cache[self.value_to_switch].__name__
        out_str = ("%s changed to %s" % (self.value_to_switch, name))
        return out_str

    def revert(self):
        """ revert the default xcorr function to previous value """
        # Have to use the previous value as this may contain some custom
        # stream_xcorr functions
        self.cache[self.value_to_switch] = self.previous_value


set_xcorr = _Context(XCOR_FUNCS, 'default')


# ---------------------- generic concurrency functions

@contextlib.contextmanager
def pool_boy(Pool, traces, **kwargs):
    """
    A context manager for handling the setup and cleanup of a pool object.

    :param Pool: any Class (not instance) that implements the multiprocessing
        Pool interface
    :param traces: The number of traces to process
    :type traces: int
    """
    # All parallel processing happens on a per-trace basis, we shouldn't create
    # more workers than there are traces
    n_cores = kwargs.get('cores', cpu_count())
    if n_cores is None:
        n_cores = cpu_count()
    if n_cores > traces:
        n_cores = traces
    pool = Pool(n_cores)
    yield pool
    pool.close()
    pool.join()


def _pool_normxcorr(templates, stream, stack, cc_squared, pool, func, *args, **kwargs):
    chans = [[] for _i in range(len(templates))]
    array_dict_tuple = _get_array_dicts(templates, stream, stack=stack)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    # get parameter iterator
    params = ((template_dict[sid], stream_dict[sid], pad_dict[sid], cc_squared)
              for sid in seed_ids)
    # get cc results and used chans into their own lists
    results = [pool.apply_async(func, param) for param in params]
    try:
        xcorrs, tr_chans = zip(*(res.get() for res in results))
    except KeyboardInterrupt as e:  # pragma: no cover
        pool.terminate()
        raise e
    if stack:
        cccsums = np.sum(xcorrs, axis=0)
    else:
        cccsums = np.asarray(xcorrs).swapaxes(0, 1)
    no_chans = np.sum(np.array(tr_chans).astype(int), axis=0)
    for seed_id, tr_chan in zip(seed_ids, tr_chans):
        for chan, state in zip(chans, tr_chan):
            if state:
                chan.append(seed_id)
    # Need to cope with possibility that earliest channel is unused. In which
    # case we need to pad the ccccsums for that by the pad for that otherwise
    # we get the wrong detection time.
    cccsums, pad_dict = _cope_with_unused_earliest(cccsums, pad_dict, chans)
    if stack:
        cccsums = _zero_invalid_correlation_sums(cccsums, pad_dict, chans)
    chans = [[(seed_id.split('.')[1], seed_id.split('.')[-1].split('_')[0])
              for seed_id in _chans] for _chans in chans]
    return cccsums, no_chans, chans


def _general_multithread(func):
    """ return the general multithreading function using func """

    def multithread(templates, stream, stack=True, cc_squared=False, *args, **kwargs):
        with pool_boy(ThreadPool, len(stream), **kwargs) as pool:
            return _pool_normxcorr(
                templates, stream, stack=stack, pool=pool, cc_squared=cc_squared,
                func=func)

    return multithread


def _general_multiprocess(func):
    def multiproc(templates, stream, stack=True, cc_squared=False, *args, **kwargs):
        with pool_boy(ProcessPool, len(stream), **kwargs) as pool:
            return _pool_normxcorr(
                templates, stream, stack=stack, pool=pool, cc_squared=cc_squared,
                func=func)

    return multiproc


def _general_serial(func):
    def stream_xcorr(templates, stream, stack=True, cc_squared=False, *args, **kwargs):
        no_chans = np.zeros(len(templates))
        chans = [[] for _ in range(len(templates))]
        array_dict_tuple = _get_array_dicts(templates, stream, stack=stack)
        stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
        if stack:
            cccsums = np.zeros([len(templates),
                                len(stream[0]) - len(templates[0][0]) + 1])
        else:
            cccsums = np.zeros([len(templates), len(seed_ids),
                                len(stream[0]) - len(templates[0][0]) + 1])
        for chan_no, seed_id in enumerate(seed_ids):
            tr_cc, tr_chans = func(template_dict[seed_id],
                                   stream_dict[seed_id],
                                   pad_dict[seed_id],
                                   cc_squared)
            if stack:
                cccsums = np.sum([cccsums, tr_cc], axis=0)
            else:
                cccsums[:, chan_no] = tr_cc
            no_chans += tr_chans.astype(int)
            for chan, state in zip(chans, tr_chans):
                if state:
                    chan.append(seed_id)
        # Need to cope with possibility that earliest channel is unused. In which
        # case we need to pad the ccccsums for that by the pad for that otherwise
        # we get the wrong detection time.
        cccsums, pad_dict = _cope_with_unused_earliest(cccsums, pad_dict, chans)
        if stack:
            cccsums = _zero_invalid_correlation_sums(cccsums, pad_dict, chans)
        chans = [[(seed_id.split('.')[1], seed_id.split('.')[-1].split('_')[0])
                  for seed_id in _chans] for _chans in chans]
        return cccsums, no_chans, chans

    return stream_xcorr


def register_array_xcorr(name, func=None, is_default=False):
    """
    Decorator for registering correlation functions.

    Each function must have the same interface as numpy_normxcorr, which is
    *f(templates, stream, pads, *args, **kwargs)* any number of specific kwargs
    can be used.

    Register_normxcorr can be used as a decorator (with or without arguments)
    or as a callable.

    :param name: The name of the function for quick access, or the callable
        that will be wrapped when used as a decorator.
    :type name: str, callable
    :param func: The function to register
    :type func: callable, optional
    :param is_default: True if this function should be marked as default
        normxcorr
    :type is_default: bool

    :return: callable
    """
    valid_methods = set(list(XCOR_ARRAY_METHODS) + list(XCORR_STREAM_METHODS))
    cache = {}

    def register(register_str):
        """
        Register a function as an implementation.

        :param register_str: The registration designation
        :type register_str: str
        """
        if register_str not in valid_methods:
            msg = 'register_name must be in %s' % valid_methods
            raise ValueError(msg)

        def _register(func):
            cache[register_str] = func
            setattr(cache['func'], register_str, func)
            return func

        return _register

    def wrapper(func, func_name=None):
        # register the functions in the XCOR
        fname = func_name or name.__name__ if callable(name) else str(name)
        XCOR_FUNCS[fname] = func
        # if is_default:  # set function as default
        #     XCOR_FUNCS['default'] = func
        # attach some attrs, this is a bit of a hack to avoid pickle problems
        func.register = register
        cache['func'] = func
        func.multithread = _general_multithread(func)
        func.multiprocess = _general_multiprocess(func)
        func.concurrent = _general_multithread(func)
        func.stream_xcorr = _general_serial(func)
        func.array_xcorr = func
        func.registered = True
        if is_default:  # set function as default
            XCOR_FUNCS['default'] = copy.deepcopy(func)
        return func

    # used as a decorator
    if callable(name):
        return wrapper(name)

    # used as a normal function (called and passed a function)
    if callable(func):
        return wrapper(func, func_name=name)

    # called, then used as a decorator
    return wrapper


def _zero_invalid_correlation_sums(cccsums, pad_dict, used_seed_ids):
    """
    Zero the end portion of the correlation sum that does not include the full
    set of channels.

    :param cccsums: 2D numpy array
    :type cccsums: np.ndarray
    :param pad_dict: pad_dict from _get_array_dicts
    :type pad_dict: dict
    :param used_seed_ids: The SEED IDs actually used in correlation
    :type used_seed_ids: list of list of str

    :return:
        Valid correlation stack - end will be zero-ed up to maxmimum moveout
    :rtype: np.ndarray
    """
    # TODO: This is potentially quite a slow way to do this.
    for i, cccsum in enumerate(cccsums):
        moveouts = [value[i] for key, value in pad_dict.items()
                    if key in used_seed_ids[i]]
        max_moveout = 0
        if len(moveouts):
            max_moveout = max(moveouts)
        if max_moveout:
            cccsum[-max_moveout:] = 0.0
        Logger.debug(f"Setting the last {max_moveout} samples to "
                     f"0 for {used_seed_ids[i]}")
    return cccsums

# ------------------ array_xcorr fetching functions


def _get_registerd_func(name_or_func):
    """ get a xcorr function from a str or callable. """
    # get the function or register callable
    if callable(name_or_func):
        func = register_array_xcorr(name_or_func)
    else:
        func = XCOR_FUNCS[name_or_func or 'default']
    assert callable(func), 'func is not callable'
    # ensure func has the added methods
    if not hasattr(func, 'registered'):
        func = register_array_xcorr(func)
    return func


def get_array_xcorr(name_or_func=None):
    """
    Get an normalized cross correlation function that takes arrays as inputs.

    See :func:`eqcorrscan.utils.correlate.array_normxcorr` for expected
    function signature.

    :param name_or_func: Either a name of a registered xcorr function or a
        callable that implements the standard array_normxcorr signature.
    :type name_or_func: str or callable

    :return: callable wth array_normxcorr interface

    see also :func:`eqcorrscan.utils.correlate.get_stream_xcorr`
    """
    func = _get_registerd_func(name_or_func)
    return func.array_xcorr


# ----------------------- registered array_xcorr functions


@register_array_xcorr('numpy')
def numpy_normxcorr(templates, stream, pads, cc_squared, *args, **kwargs):
    """
    Compute the normalized cross-correlation using numpy and bottleneck.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    import bottleneck

    # Generate a template mask
    used_chans = ~np.isnan(templates).any(axis=1)
    # Currently have to use float64 as bottleneck runs into issues with other
    # types: https://github.com/kwgoodman/bottleneck/issues/164
    stream = stream.astype(np.float64)
    templates = templates.astype(np.float64)
    template_length = templates.shape[1]
    stream_length = len(stream)
    assert stream_length > template_length, "Template must be shorter than " \
                                            "stream"
    fftshape = next_fast_len(template_length + stream_length - 1)
    # Set up normalizers
    stream_mean_array = bottleneck.move_mean(
        stream, template_length)[template_length - 1:]
    stream_std_array = bottleneck.move_std(
        stream, template_length)[template_length - 1:]
    # because stream_std_array is in denominator or res, nan all 0s
    stream_std_array[stream_std_array == 0] = np.nan
    # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))
    norm_sum = norm.sum(axis=-1, keepdims=True)
    stream_fft = np.fft.rfft(stream, fftshape)
    template_fft = np.fft.rfft(np.flip(norm, axis=-1), fftshape, axis=-1)
    res = np.fft.irfft(template_fft * stream_fft,
                       fftshape)[:, 0:template_length + stream_length - 1]
    res = ((_centered(res, (templates.shape[0],
                            stream_length - template_length + 1))) -
           norm_sum * stream_mean_array) / stream_std_array
    res[np.isnan(res)] = 0.0

    if cc_squared:
        res *= np.abs(res)

    for i, pad in enumerate(pads):
        res[i] = np.append(res[i], np.zeros(pad))[pad:]
    return res.astype(np.float32), used_chans


def _centered(arr, newshape):
    """
    Hack of scipy.signaltools._centered
    """
    newshape = np.asarray(newshape)
    currshape = np.array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


@register_array_xcorr('time_domain')
def time_multi_normxcorr(templates, stream, pads, cc_squared,
                         threaded=False, *args, **kwargs):
    """
    Compute cross-correlations in the time-domain using C routine.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool
    :param threaded: Whether to use the threaded routine or not
    :type threaded: bool

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    used_chans = ~np.isnan(templates).any(axis=1)

    utilslib = _load_cdll('libutils')

    argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int]
    restype = ctypes.c_int
    if threaded:
        func = utilslib.multi_normxcorr_time_threaded
        argtypes.append(ctypes.c_int)
    else:
        func = utilslib.multi_normxcorr_time
    func.argtypes = argtypes
    func.restype = restype
    # Need to de-mean everything
    templates_means = templates.mean(axis=1).astype(np.float32)[:, np.newaxis]
    stream_mean = stream.mean().astype(np.float32)
    templates = templates.astype(np.float32) - templates_means
    stream = stream.astype(np.float32) - stream_mean
    template_len = templates.shape[1]
    n_templates = templates.shape[0]
    image_len = stream.shape[0]
    ccc_length = image_len - template_len + 1
    assert ccc_length > 0, "Template must be shorter than stream"
    ccc = np.ascontiguousarray(
        np.empty(ccc_length * n_templates), np.float32)
    t_array = np.ascontiguousarray(templates.flatten(), np.float32)
    time_args = [t_array, template_len, n_templates,
                 np.ascontiguousarray(stream, np.float32), image_len, ccc, int(cc_squared)]
    if threaded:
        time_args.append(kwargs.get('cores', cpu_count()))
    func(*time_args)
    ccc[np.isnan(ccc)] = 0.0
    ccc = ccc.reshape((n_templates, image_len - template_len + 1))
    for i in range(len(pads)):
        ccc[i] = np.append(ccc[i], np.zeros(pads[i]))[pads[i]:]
    templates += templates_means
    stream += stream_mean
    return ccc, used_chans


@register_array_xcorr('fftw', is_default=True)
def fftw_normxcorr(templates, stream, pads, cc_squared,
                   threaded=False, *args, **kwargs):
    """
    Normalised cross-correlation using the fftw library.

    Internally this function used double precision numbers, which is definitely
    required for seismic data. Cross-correlations are computed as the
    inverse fft of the dot product of the ffts of the stream and the reversed,
    normalised, templates.  The cross-correlation is then normalised using the
    running mean and standard deviation (not using the N-1 correction) of the
    stream and the sums of the normalised templates.

    This python function wraps the C-library written by C. Chamberlain for this
    purpose.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool
    :param threaded:
        Whether to use the threaded routine or not - note openMP and python
        multiprocessing don't seem to play nice for this.
    :type threaded: bool

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    utilslib = _load_cdll('libutils')

    argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long, ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int]
    restype = ctypes.c_int

    if threaded:
        func = utilslib.normxcorr_fftw_threaded
    else:
        func = utilslib.normxcorr_fftw

    func.argtypes = argtypes
    func.restype = restype

    # Generate a template mask
    used_chans = ~np.isnan(templates).any(axis=1)
    template_length = templates.shape[1]
    stream_length = len(stream)
    n_templates = templates.shape[0]
    fftshape = kwargs.get("fft_len")
    if fftshape is None:
        # In testing, 2**13 consistently comes out fastest - setting to
        # default. https://github.com/eqcorrscan/EQcorrscan/pull/285
        fftshape = min(
            2 ** 13, next_fast_len(template_length + stream_length - 1))
    if fftshape < template_length:
        Logger.warning(
            "FFT length of {0} is shorter than the template, setting to "
            "{1}".format(
                fftshape, next_fast_len(template_length + stream_length - 1)))
        fftshape = next_fast_len(template_length + stream_length - 1)
    # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))

    norm = np.nan_to_num(norm)
    ccc_length = stream_length - template_length + 1
    assert ccc_length > 0, "Template must be shorter than stream"
    ccc = np.zeros((n_templates, ccc_length), np.float32)
    used_chans_np = np.ascontiguousarray(used_chans, dtype=np.intc)
    pads_np = np.ascontiguousarray(pads, dtype=np.intc)
    variance_warning = np.ascontiguousarray([0], dtype=np.intc)
    missed_corr = np.ascontiguousarray([0], dtype=np.intc)

    # Check that stream is non-zero and above variance threshold
    if not np.all(stream == 0) and np.var(stream) < 1e-8:
        # Apply gain
        stream *= MULTIPLIER
        Logger.warning("Low variance found for, applying gain "
                       "to stabilise correlations")
        multiplier = MULTIPLIER
    else:
        multiplier = 1
    ret = func(
        np.ascontiguousarray(norm.flatten(order='C'), np.float32),
        template_length, n_templates,
        np.ascontiguousarray(stream, np.float32), stream_length,
        np.ascontiguousarray(ccc, np.float32), fftshape,
        used_chans_np, pads_np, variance_warning, missed_corr,
        int(cc_squared))
    if ret < 0:
        raise MemoryError()
    elif ret > 0:
        Logger.critical('Error in C code (possible normalisation error)')
        Logger.critical('Maximum ccc %f at %i' % (ccc.max(), ccc.argmax()))
        Logger.critical('Minimum ccc %f at %i' % (ccc.min(), ccc.argmin()))
        Logger.critical('Recommend checking your data for spikes, clipping '
                        'or artefacts')
        raise CorrelationError("Internal correlation error")
    if missed_corr[0]:
        Logger.warning(
            "{0} correlations not computed, are there gaps in the "
            "data? If not, consider increasing gain".format(
                missed_corr[0]))
    if variance_warning[0] and variance_warning[0] > template_length:
        Logger.warning(
            "Low variance found in {0} positions, check result.".format(
                variance_warning[0]))
    # Remove variance correction
    stream /= multiplier
    return ccc, used_chans


# The time-domain routine can be sped up massively on large machines (many
# threads) using the openMP threaded functions.

@time_multi_normxcorr.register('concurrent')
def _time_threaded_normxcorr(templates, stream, stack=True, cc_squared=False, *args, **kwargs):
    """
    Use the threaded time-domain routine for concurrency

    :type templates: list
    :param templates:
        A list of templates, where each one should be an obspy.Stream object
        containing multiple traces of seismic data and the relevant header
        information.
    :type stream: obspy.core.stream.Stream
    :param stream:
        A single Stream object to be correlated with the templates.
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool

    :returns:
        New list of :class:`numpy.ndarray` objects.  These will contain
        the correlation sums for each template for this day of data.
    :rtype: list
    :returns:
        list of ints as number of channels used for each cross-correlation.
    :rtype: list
    :returns:
        list of list of tuples of station, channel for all cross-correlations.
    :rtype: list
    """
    no_chans = np.zeros(len(templates))
    chans = [[] for _ in range(len(templates))]
    array_dict_tuple = _get_array_dicts(templates, stream, stack=stack)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    ccc_length = max(
        len(stream[0]) - len(templates[0][0]) + 1,
        len(templates[0][0]) - len(stream[0]) + 1)
    if stack:
        cccsums = np.zeros([len(templates), ccc_length])
    else:
        cccsums = np.zeros([len(templates), len(seed_ids), ccc_length])
    for chan_no, seed_id in enumerate(seed_ids):
        tr_cc, tr_chans = time_multi_normxcorr(
            template_dict[seed_id], stream_dict[seed_id], pad_dict[seed_id],
            cc_squared=cc_squared, threaded=True)
        if stack:
            cccsums = np.sum([cccsums, tr_cc], axis=0)
        else:
            cccsums[:, chan_no] = tr_cc
        no_chans += tr_chans.astype(int)
        for chan, state in zip(chans, tr_chans):
            if state:
                chan.append(seed_id)
    # Need to cope with possibility that earliest channel is unused. In which
    # case we need to pad the ccccsums for that by the pad for that otherwise
    # we get the wrong detection time.
    cccsums, pad_dict = _cope_with_unused_earliest(cccsums, pad_dict, chans)
    if stack:
        cccsums = _zero_invalid_correlation_sums(cccsums, pad_dict, chans)
    chans = [[(seed_id.split('.')[1], seed_id.split('.')[-1].split('_')[0])
              for seed_id in _chans] for _chans in chans]
    return cccsums, no_chans, chans


@fftw_normxcorr.register('stream_xcorr')
@fftw_normxcorr.register('multithread')
@fftw_normxcorr.register('concurrent')
def _fftw_stream_xcorr(templates, stream, stack=True, cc_squared=False, *args, **kwargs):
    """
    Apply fftw normxcorr routine concurrently.

    :type templates: list
    :param templates:
        A list of templates, where each one should be an obspy.Stream object
        containing multiple traces of seismic data and the relevant header
        information.
    :type stream: obspy.core.stream.Stream
    :param stream:
        A single Stream object to be correlated with the templates.
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool

    :returns:
        New list of :class:`numpy.ndarray` objects.  These will contain
        the correlation sums for each template for this day of data.
    :rtype: list
    :returns:
        list of ints as number of channels used for each cross-correlation.
    :rtype: list
    :returns:
        list of list of tuples of station, channel for all cross-correlations.
    :rtype: list
    """
    array_dict_tuple = _get_array_dicts(templates, stream, stack=stack)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    assert set(seed_ids)
    # number of threads:
    #   default to using inner threads
    #   if `cores` or `cores_outer` passed in then use that
    #   else if OMP_NUM_THREADS set use that
    #   otherwise use all available
    inner_kwargs = copy.copy(kwargs)
    num_cores_inner, num_cores_outer = _set_inner_outer_threading(
        inner_kwargs.pop('cores', None),
        inner_kwargs.pop("cores_outer", None),
        len(stream))

    chans = [[] for _i in range(len(templates))]
    cccsums, tr_chans = fftw_multi_normxcorr(
        template_array=template_dict, stream_array=stream_dict,
        pad_array=pad_dict, seed_ids=seed_ids, cores_inner=num_cores_inner,
        cores_outer=num_cores_outer, stack=stack, cc_squared=cc_squared,
        *args, **inner_kwargs)
    no_chans = np.sum(np.array(tr_chans).astype(int), axis=0)
    for seed_id, tr_chan in zip(seed_ids, tr_chans):
        for chan, state in zip(chans, tr_chan):
            if state:
                chan.append(seed_id)
    # Need to cope with possibility that earliest channel is unused. In which
    # case we need to pad the ccccsums for that by the pad for that otherwise
    # we get the wrong detection time.
    cccsums, pad_dict = _cope_with_unused_earliest(cccsums, pad_dict, chans)

    if stack:
        cccsums = _zero_invalid_correlation_sums(cccsums, pad_dict, chans)

    # Reshape to (station, channel)
    chans = [[(seed_id.split('.')[1], seed_id.split('.')[-1].split('_')[0])
              for seed_id in _chans] for _chans in chans]

    return cccsums, no_chans, chans


def _rm_unused_pads(pads, sids, offset=0):
    pad_sids = list(pads.keys())
    pad_values = [pads[k] for k in pad_sids]
    # Find the lowest pad in the set and get the key for it
    min_pad_sid = pad_sids[pad_values.index(min(pad_values))]
    # If the key is in the sids, we are done
    if min_pad_sid in sids:
        return offset, pads
    Logger.debug(f"Min pad: {min(pad_values)} on {min_pad_sid}, but not "
                 f"used, moving to next channel")
    # Otherwise, we need to remove that from pads and reduce the pads
    # by the delta between the min pad and the new min pad
    min_pad_value = pads.pop(min_pad_sid)
    if len(pads) == 0:
        return offset, pads
    new_min_pad_value = min(pads.values())
    # Add offset to total offset
    offset += new_min_pad_value - min_pad_value
    pads = {sid: v - (new_min_pad_value - min_pad_value)
            for sid, v in pads.items()}
    # Need to recurse
    offset, pads = _rm_unused_pads(pads=pads, sids=sids, offset=offset)
    return offset, pads


def _cope_with_unused_earliest(cccsums, pad_dict, chans):
    pad_dict_out = {k: [] for k in pad_dict.keys()}
    for i, ccsum in enumerate(cccsums):
        # Check if the earliest pad is in chans
        pads = {sid: pads[i] for sid, pads in pad_dict.items()}
        Logger.debug(f"pads {pads}")
        offset, edited_pads = _rm_unused_pads(pads, chans[i])
        # Put shuffled pads back into pads out
        for key, value in pad_dict_out.items():
            value.append(edited_pads.get(key, 0))
        if offset:
            Logger.debug(f"Rolling cccsum array by {offset}")
            # Roll along the axis by offset
            cccsums[i] = np.roll(ccsum, offset)
            # Set offset samples to zero
            cccsums[i][0:offset] = 0
    return cccsums, pad_dict_out


def _set_inner_outer_threading(num_cores_inner, num_cores_outer, n_chans):
    max_threads = int(os.getenv("OMP_NUM_THREADS", cpu_count()))
    if num_cores_inner is None:
        num_cores_inner = max_threads
    num_cores_outer = num_cores_outer or 1
    if num_cores_outer > 1:
        if num_cores_outer > n_chans:
            Logger.info(
                "More outer cores than channels, setting to {0}".format(
                    n_chans))
            num_cores_outer = n_chans
        if num_cores_outer * num_cores_inner > max_threads:
            Logger.info("More threads requested than exist, falling back to "
                        "outer-loop parallelism")
            num_cores_outer = min(max_threads, num_cores_outer)
            if 2 * num_cores_outer < max_threads:
                num_cores_inner = max_threads // num_cores_outer
            else:
                num_cores_inner = 1
    return num_cores_inner, num_cores_outer


def fftw_multi_normxcorr(
    template_array,
    stream_array,
    pad_array,
    seed_ids,
    cores_inner,
    cores_outer,
    stack=True,
    cc_squared=False,
    *args,
    **kwargs
):
    """
    Use a C loop rather than a Python loop - in some cases this will be fast.

    :type template_array: dict
    :param template_array:
    :type stream_array: dict
    :param stream_array:
    :type pad_array: dict
    :param pad_array:
    :type seed_ids: list
    :param seed_ids:
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool

    rtype: np.ndarray, list
    :return: 3D Array of cross-correlations and list of used channels.
    """
    utilslib = _load_cdll('libutils')

    utilslib.multi_normxcorr_fftw.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long, ctypes.c_long, ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_int]
    utilslib.multi_normxcorr_fftw.restype = ctypes.c_int
    '''
    Arguments are:
        templates (stacked [ch_1-t_1, ch_1-t_2, ..., ch_2-t_1, ch_2-t_2, ...])
        number of templates
        template length
        number of channels
        image (stacked [ch_1, ch_2, ..., ch_n])
        image length
        cross-correlations (stacked as per image)
        fft-length
        used channels (stacked as per templates)
        pad array (stacked as per templates)
        num threads inner
        num threads outer
        variance warnings
        missed correlation warnings (usually due to gaps)
        stack option
        cc squared option
    '''

    # pre processing
    used_chans = []
    template_len = template_array[seed_ids[0]].shape[1]
    for seed_id in seed_ids:
        used_chans.append(~np.isnan(template_array[seed_id]).any(axis=1))
        template_array[seed_id] = (
            (template_array[seed_id] -
             template_array[seed_id].mean(axis=-1, keepdims=True)) / (
                template_array[seed_id].std(axis=-1, keepdims=True) *
                template_len))
        template_array[seed_id] = np.nan_to_num(template_array[seed_id])
    n_channels = len(seed_ids)
    n_templates = template_array[seed_ids[0]].shape[0]
    image_len = stream_array[seed_ids[0]].shape[0]
    fft_len = kwargs.get("fft_len")
    if fft_len is None:
        # In testing, 2**13 consistently comes out fastest - setting to
        # default. https://github.com/eqcorrscan/EQcorrscan/pull/285
        # But this results in lots of chunks - 2 ** 17 is also good.
        fft_len = min(2 ** 17, next_fast_len(template_len + image_len - 1))
    if fft_len < template_len:
        Logger.warning(
            f"FFT length of {fft_len} is shorter than the template, setting to"
            f" {next_fast_len(template_len + image_len - 1)}")
        fft_len = next_fast_len(template_len + image_len - 1)
    template_array = np.ascontiguousarray(
        [template_array[x] for x in seed_ids], dtype=np.float32)
    multipliers = {}
    for x in seed_ids:
        # Check that stream is non-zero and above variance threshold
        if not np.all(stream_array[x] == 0) and np.var(stream_array[x]) < 1e-8:
            # Apply gain
            stream_array[x] *= MULTIPLIER
            Logger.warning(f"Low variance found for {x}, applying gain "
                           "to stabilise correlations")
            multipliers.update({x: MULTIPLIER})
        else:
            multipliers.update({x: 1})
    stream_array = np.ascontiguousarray([stream_array[x] for x in seed_ids],
                                        dtype=np.float32)
    ccc_length = image_len - template_len + 1
    assert ccc_length > 0, "Template must be shorter than stream"
    if stack:
        cccs = np.zeros((n_templates, ccc_length), np.float32)
    else:
        cccs = np.zeros(
            (n_templates, n_channels, ccc_length), dtype=np.float32)
    used_chans_np = np.ascontiguousarray(used_chans, dtype=np.intc)
    pad_array_np = np.ascontiguousarray(
        [pad_array[seed_id] for seed_id in seed_ids], dtype=np.intc)
    variance_warnings = np.ascontiguousarray(
        np.zeros(n_channels), dtype=np.intc)
    missed_correlations = np.ascontiguousarray(
        np.zeros(n_channels), dtype=np.intc)

    # call C function
    ret = utilslib.multi_normxcorr_fftw(
        template_array, n_templates, template_len, n_channels, stream_array,
        image_len, cccs, fft_len, used_chans_np, pad_array_np,
        cores_inner, cores_outer, variance_warnings, missed_correlations,
        int(stack), int(cc_squared))
    if ret < 0:
        raise MemoryError("Memory allocation failed in correlation C-code")
    elif ret > 0:
        Logger.critical(
            'Out-of-range correlation in C-code, see WARNING from C-code.'
            'You are STRONGLY RECOMMENDED to check your data for spikes, '
            'clipping or non-physical artifacts')
        # raise CorrelationError("Internal correlation error")
    for i, missed_corr in enumerate(missed_correlations):
        if missed_corr:
            Logger.debug(
                f"{missed_corr} correlations not computed on {seed_ids[i]}, "
                f"are there gaps in the data? If not, consider "
                "increasing gain")
    for i, variance_warning in enumerate(variance_warnings):
        if variance_warning and variance_warning > template_len:
            Logger.warning(
                f"Low variance found in {variance_warning} places for "
                f"{seed_ids[i]}, check result.")
    # Remove gain
    for i, x in enumerate(seed_ids):
        stream_array[i] *= multipliers[x]
    return cccs, used_chans


# ------------------------------- FastMatchedFilter Wrapper

def _run_fmf_xcorr(template_arr, data_arr, weights, pads, arch, step=1):
    template_arr, data_arr, multipliers = _fmf_stabilisation(
        template_arr=template_arr, data_arr=data_arr)

    return _stabalised_fmf(
        template_arr=template_arr, data_arr=data_arr, weights=weights,
        pads=pads, arch=arch, multipliers=multipliers, step=step)


def _stabalised_fmf(template_arr, data_arr, weights, pads, arch, multipliers,
                    step):
    if not FMF_INSTALLED:
        raise ImportError("FastMatchedFilter is not available")
    import fast_matched_filter

    if version.parse(fast_matched_filter.__version__) >= MIN_FMF_VERSION:
        from fast_matched_filter import matched_filter as fmf
    else:
        raise ImportError(f"FMF version {fast_matched_filter.__version__} "
                          f"must be >= {MIN_FMF_VERSION}")

    Logger.info("Handing off to FMF")
    cccsums = fmf(
        templates=template_arr, weights=weights, moveouts=pads,
        data=data_arr, step=step, arch=arch, normalize="full")
    Logger.info("FMF returned")
    # Remove gain
    if np.any(multipliers != 1):
        data_arr /= multipliers
    return cccsums


def _fmf_stabilisation(template_arr, data_arr):
    """ FMF doesn't do the internal check that EQC correlations do. """
    # Demean
    tic = default_timer()
    template_arr -= template_arr.mean(axis=-1, keepdims=True)
    data_arr -= data_arr.mean(axis=-1, keepdims=True)
    toc = default_timer()
    Logger.info(f"Removing mean took {toc - tic:.4f} s")

    # Stability checking
    tic = default_timer()
    # var is fairly slow, var = mean(abs(a - a.mean()) ** 2) - mean is zero,
    # so we can skip a step
    stability_issues = np.logical_and(
        # data_arr.var(axis=1, keepdims=True) < 1e-8,
        np.mean(np.abs(data_arr) ** 2, axis=1, keepdims=True) < 1e-8,
        ~np.all(data_arr == 0, axis=1, keepdims=True))
    multipliers = np.ones_like(stability_issues, dtype=float)
    multipliers[stability_issues] = MULTIPLIER
    if np.any(stability_issues):
        Logger.warning(
            f"Low variance found for channels {np.where(stability_issues)},"
            f"applying gain to stabilise correlations")
        data_arr *= multipliers
    toc = default_timer()
    Logger.info(f"Checking stability took {toc - tic:.4f} s")
    return template_arr, data_arr, multipliers


@register_array_xcorr("fmf")
def fmf_xcorr(templates, stream, pads, arch="precise", cc_squared=False, *args, **kwargs):
    """
    Compute cross-correlations in the time-domain using the FMF routine.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list
    :param arch: "gpu" or "precise" to run on GPU or CPU respectively
    type arch: str
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    assert templates.ndim == 2, "Templates must be 2D"
    assert stream.ndim == 1, "Stream must be 1D"

    used_chans = ~np.isnan(templates).any(axis=1)

    # We have to reshape to an extra dimension for FMF
    ccc = _run_fmf_xcorr(
        template_arr=templates.reshape(
            (1, templates.shape[0], templates.shape[1])).swapaxes(0, 1),
        data_arr=stream.reshape((1, stream.shape[0])),
        weights=np.ones((1, templates.shape[0])),
        pads=np.array([pads]),
        arch=arch.lower())
    if cc_squared:
        ccc *= np.abs(ccc)

    return ccc, used_chans


@fmf_xcorr.register("stream_xcorr")
@fmf_xcorr.register("concurrent")
def _fmf_gpu(templates, stream, cc_squared=False, *args, **kwargs):
    """
    Thin wrapper of fmf_multi_xcorr setting arch to gpu.
    """
    from fast_matched_filter import GPU_LOADED
    if not GPU_LOADED:
        Logger.warning("FMF reports GPU not loaded, reverting to CPU")
        return _fmf_cpu(templates=templates, stream=stream,
                        cc_squared=cc_squared*args, **kwargs)
    return _fmf_multi_xcorr(templates, stream, cc_squared=cc_squared, arch="gpu")


@fmf_xcorr.register("multithread")
@fmf_xcorr.register("multiprocess")
def _fmf_cpu(templates, stream, cc_squared=False, *args, **kwargs):
    """
    Thin wrapper of fmf_multi_xcorr setting arch to cpu.
    """
    from fast_matched_filter import CPU_LOADED
    if not CPU_LOADED:
        raise NotImplementedError(
            "FMF reports CPU not loaded - try rebuilding FMF")
    return _fmf_multi_xcorr(templates, stream, cc_squared=cc_squared,
                            arch="precise")


def _fmf_reshape(template_dict, stream_dict, pad_dict, seed_ids):
    tic = default_timer()
    # Reshape templates into [templates x traces x time]
    t_arr = np.array([template_dict[seed_id]
                      for seed_id in seed_ids]).swapaxes(0, 1)
    # Reshape stream into [traces x time]
    d_arr = np.array([stream_dict[seed_id] for seed_id in seed_ids])
    # Moveouts should be [templates x traces]
    pads = np.array([pad_dict[seed_id] for seed_id in seed_ids]).swapaxes(0, 1)
    # Weights should be shaped like pads
    weights = np.ones_like(pads)
    toc = default_timer()
    Logger.info(f"Reshaping for FMF took {toc - tic:.4f} s")
    return t_arr, d_arr, weights, pads


def _fmf_multi_xcorr(templates, stream, cc_squared=False, *args, **kwargs):
    """
    Apply FastMatchedFilter routine concurrently.

    :type templates: list
    :param templates:
        A list of templates, where each one should be an obspy.Stream object
        containing multiple traces of seismic data and the relevant header
        information.
    :type stream: obspy.core.stream.Stream
    :param stream:
        A single Stream object to be correlated with the templates.
    :param cc_squared: Whether to output cc-squared or not
    :type cc_squared: bool

    :returns:
        New list of :class:`numpy.ndarray` objects.  These will contain
        the correlation sums for each template for this day of data.
    :rtype: list
    :returns:
        list of ints as number of channels used for each cross-correlation.
    :rtype: list
    :returns:
        list of list of tuples of station, channel for all cross-correlations.
    :rtype: list
    """
    if kwargs.get("stack", False):
        raise NotImplementedError(
            "FMF does not support unstacked correlations, use a different "
            "backend")
    arch = kwargs.get("arch", "gpu")
    Logger.info(f"Running FMF targeting the {arch}")

    chans = [[] for _i in range(len(templates))]
    array_dict_tuple = _get_array_dicts(templates, stream, stack=True)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    assert set(seed_ids)

    t_arr, d_arr, weights, pads = _fmf_reshape(
        template_dict=template_dict, stream_dict=stream_dict,
        pad_dict=pad_dict, seed_ids=seed_ids)

    cccsums = _run_fmf_xcorr(
        template_arr=t_arr, weights=weights, pads=pads,
        data_arr=d_arr, step=1, arch=arch)

    tr_chans = np.array([~np.isnan(template_dict[seed_id]).any(axis=1)
                         for seed_id in seed_ids])
    no_chans = np.sum(np.array(tr_chans).astype(int), axis=0)
    # Note: FMF already returns the zeroed end of correlations - we don't
    # need to call _get_valid_correlation_sum
    used_sids = [[] for _i in range(len(templates))]
    for seed_id, tr_chan in zip(seed_ids, tr_chans):
        for chan, used_sid, state in zip(chans, used_sids, tr_chan):
            if state:
                chan.append((seed_id.split('.')[1],
                             seed_id.split('.')[-1].split('_')[0]))
                used_sid.append(seed_id)
    # Need to cope with possibility that earliest channel is unused. In which
    # case we need to pad the ccccsums for that by the pad for that otherwise
    # we get the wrong detection time.
    cccsums, pad_dict = _cope_with_unused_earliest(cccsums, pad_dict, used_sids)
    if cc_squared:
        cccsums *= np.abs(cccsums)
    return cccsums, no_chans, chans


# ------------------------------- stream_xcorr functions


def get_stream_xcorr(name_or_func=None, concurrency=None):
    """
    Return a function for performing normalized cross correlation on lists of
    streams.

    :param name_or_func:
        Either a name of a registered function or a callable that implements
        the standard array_normxcorr signature.
    :param concurrency:
        Optional concurrency strategy, options are below.

    :return: A callable with the interface of stream_normxcorr

    :Concurrency options:
        - multithread - use a threadpool for concurrency;
        - multiprocess - use a process pool for concurrency;
        - concurrent - use a customized concurrency strategy for the function,
          if not defined threading will be used.
    """
    func = _get_registerd_func(name_or_func)

    concur = concurrency or 'stream_xcorr'
    if not hasattr(func, concur):
        msg = '%s does not support concurrency %s' % (func.__name__, concur)
        raise ValueError(msg)
    outfunc = getattr(func, concur)
    Logger.info(f"Found {outfunc} for {name_or_func} - {concurrency}")
    return outfunc


# --------------------------- stream prep functions


def _get_array_dicts(templates, stream, stack, *args, **kwargs):
    """ prepare templates and stream, return dicts """
    tic = default_timer()
    # Do some reshaping
    # init empty structures for data storage
    template_dict = {}
    stream_dict = {}
    pad_dict = {}
    t_starts = []

    stream.sort(['network', 'station', 'location', 'channel'])
    for template in templates:
        template.sort(['network', 'station', 'location', 'channel'])
        t_starts.append(
            UTCDateTime(ns=min([tr.stats.starttime.__dict__['_UTCDateTime__ns']
                                for tr in template])))
    stream_start = min([tr.stats.starttime for tr in stream])
    # get seed ids, make sure these are collected on sorted streams
    seed_ids = [tr.id + '_' + str(i) for i, tr in enumerate(templates[0])]
    # pull common channels out of streams and templates and put in dicts
    for i, seed_id in enumerate(seed_ids):
        temps_with_seed = [template.traces[i].data for template in templates]
        t_ar = np.array(temps_with_seed).astype(np.float32)
        template_dict.update({seed_id: t_ar})
        stream_channel = _stream_quick_select(stream, seed_id.split('_')[0])[0]
        # Normalize data to ensure no float overflow
        stream_data = stream_channel.data / (np.max(
            np.abs(stream_channel.data)) / 1e5)
        stream_dict.update(
            {seed_id: stream_data.astype(np.float32)})
        # PROBLEM - DEBUG: if two traces start just before / just after a
        # "full-sample-time", then stream_offset can become 1, while a value in
        # pad_list can become 0. 0-1 = -1; which is problematic.
        stream_offset = int(
            math.floor(stream_channel.stats.sampling_rate *
                       (stream_channel.stats.starttime - stream_start)))
        if stack:
            pad_list = [
                int(round(
                    template.traces[i].stats.__dict__['sampling_rate'] *
                    (template.traces[i].stats.starttime.__dict__[
                        '_UTCDateTime__ns'] -
                     t_starts[j].__dict__['_UTCDateTime__ns']) / 1e9)) -
                stream_offset
                for j, template in enumerate(templates)]
        else:
            pad_list = [0 for _ in templates]
        pad_dict.update({seed_id: pad_list})
    toc = default_timer()
    Logger.info(f"Making array dicts for {len(seed_ids)} seed ids "
                f"took {toc - tic:.4f} s")

    return stream_dict, template_dict, pad_dict, seed_ids


# Remove fmf if it isn't installed
if not FMF_INSTALLED:
    XCOR_FUNCS.pop("fmf")
# a dict of built in xcorr functions, used to distinguish from user-defined
XCORR_FUNCS_ORIGINAL = copy.copy(XCOR_FUNCS)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
