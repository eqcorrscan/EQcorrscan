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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import contextlib
import copy
import ctypes
import os
from multiprocessing import Pool as ProcessPool, cpu_count
from multiprocessing.pool import ThreadPool

import numpy as np
from future.utils import native_str

from eqcorrscan.utils.libnames import _load_cdll

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
        self.previous_value = self.cache.get(self.value_to_switch)
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
        new_value = self.previous_value
        self(new_value)


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


def _pool_normxcorr(templates, stream, pool, func, *args, **kwargs):
    chans = [[] for _i in range(len(templates))]
    array_dict_tuple = _get_array_dicts(templates, stream)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    # get parameter iterator
    params = ((template_dict[sid], stream_dict[sid], pad_dict[sid])
              for sid in seed_ids)
    # get cc results and used chans into their own lists
    results = [pool.apply_async(func, param) for param in params]
    xcorrs, tr_chans = zip(*(res.get() for res in results))
    cccsums = np.sum(xcorrs, axis=0)
    no_chans = np.sum(np.array(tr_chans).astype(np.int), axis=0)
    for seed_id, tr_chan in zip(seed_ids, tr_chans):
        for chan, state in zip(chans, tr_chan):
            if state:
                chan.append((seed_id.split('.')[1],
                             seed_id.split('.')[-1].split('_')[0]))
    return cccsums, no_chans, chans


def _general_multithread(func):
    """ return the general multithreading function using func """

    def multithread(templates, stream, *args, **kwargs):
        with pool_boy(ThreadPool, len(stream), **kwargs) as pool:
            return _pool_normxcorr(templates, stream, pool=pool, func=func)

    return multithread


def _general_multiprocess(func):
    def multiproc(templates, stream, *args, **kwargs):
        with pool_boy(ProcessPool, len(stream), **kwargs) as pool:
            return _pool_normxcorr(templates, stream, pool=pool, func=func)

    return multiproc


def _general_serial(func):
    def stream_xcorr(templates, stream, *args, **kwargs):
        no_chans = np.zeros(len(templates))
        chans = [[] for _ in range(len(templates))]
        array_dict_tuple = _get_array_dicts(templates, stream)
        stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
        cccsums = np.zeros([len(templates),
                            len(stream[0]) - len(templates[0][0]) + 1])
        for seed_id in seed_ids:
            tr_cc, tr_chans = func(template_dict[seed_id],
                                   stream_dict[seed_id],
                                   pad_dict[seed_id])
            cccsums = np.sum([cccsums, tr_cc], axis=0)
            no_chans += tr_chans.astype(np.int)
            for chan, state in zip(chans, tr_chans):
                if state:
                    chan.append((seed_id.split('.')[1],
                                 seed_id.split('.')[-1].split('_')[0]))
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
        if is_default:  # set function as default
            XCOR_FUNCS['default'] = func
        # attach some attrs, this is a bit of a hack to avoid pickle problems
        func.register = register
        cache['func'] = func
        func.multithread = _general_multithread(func)
        func.multiprocess = _general_multiprocess(func)
        func.concurrent = _general_multithread(func)
        func.stream_xcorr = _general_serial(func)
        func.array_xcorr = func
        func.registered = True
        return func

    # used as a decorator
    if callable(name):
        return wrapper(name)

    # used as a normal function (called and passed a function)
    if callable(func):
        return wrapper(func, func_name=name)

    # called, then used as a decorator
    return wrapper


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
def numpy_normxcorr(templates, stream, pads, *args, **kwargs):
    """
    Compute the normalized cross-correlation using numpy and bottleneck.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    import bottleneck
    from scipy.signal.signaltools import _centered

    # Generate a template mask
    used_chans = ~np.isnan(templates).any(axis=1)
    # Currently have to use float64 as bottleneck runs into issues with other
    # types: https://github.com/kwgoodman/bottleneck/issues/164
    stream = stream.astype(np.float64)
    templates = templates.astype(np.float64)
    template_length = templates.shape[1]
    stream_length = len(stream)
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
    res = ((_centered(res, stream_length - template_length + 1)) -
           norm_sum * stream_mean_array) / stream_std_array
    res[np.isnan(res)] = 0.0
    # res[np.isinf(res)] = 0.0
    for i, pad in enumerate(pads):  # range(len(pads)):
        res[i] = np.append(res[i], np.zeros(pad))[pad:]
    return res.astype(np.float32), used_chans


@register_array_xcorr('time_domain')
def time_multi_normxcorr(templates, stream, pads, threaded=False, *args,
                         **kwargs):
    """
    Compute cross-correlations in the time-domain using C routine.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list
    :param threaded: Whether to use the threaded routine or not
    :type threaded: bool

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    used_chans = ~np.isnan(templates).any(axis=1)

    utilslib = _load_cdll('libutils')

    argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS'))]
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
    ccc = np.ascontiguousarray(
        np.empty((image_len - template_len + 1) * n_templates), np.float32)
    t_array = np.ascontiguousarray(templates.flatten(), np.float32)
    time_args = [t_array, template_len, n_templates,
                 np.ascontiguousarray(stream, np.float32), image_len, ccc]
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
def fftw_normxcorr(templates, stream, pads, threaded=False, *args, **kwargs):
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
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long, ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags=native_str('C_CONTIGUOUS'))]
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
    fftshape = next_fast_len(template_length + stream_length - 1)

    # # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))

    norm = np.nan_to_num(norm)
    ccc = np.zeros((n_templates, stream_length - template_length + 1),
                   np.float32)
    used_chans_np = np.ascontiguousarray(used_chans, dtype=np.intc)
    pads_np = np.ascontiguousarray(pads, dtype=np.intc)

    ret = func(
        np.ascontiguousarray(norm.flatten(order='C'), np.float32),
        template_length, n_templates,
        np.ascontiguousarray(stream, np.float32), stream_length,
        np.ascontiguousarray(ccc, np.float32), fftshape,
        used_chans_np, pads_np)
    if ret < 0:
        raise MemoryError()
    elif ret not in [0, 999]:
        print('Error in C code (possible normalisation error)')
        print('Maximum ccc %f at %i' % (ccc.max(), ccc.argmax()))
        print('Minimum ccc %f at %i' % (ccc.min(), ccc.argmin()))
        raise CorrelationError("Internal correlation error")
    elif ret == 999:
        msg = ("Some correlations not computed, are there "
               "zeros in data? If not, consider increasing gain.")
        print(msg)

    return ccc, used_chans


# The time-domain routine can be sped up massively on large machines (many
# threads) using the openMP threaded functions.

@time_multi_normxcorr.register('concurrent')
def _time_threaded_normxcorr(templates, stream, *args, **kwargs):
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
    array_dict_tuple = _get_array_dicts(templates, stream)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    cccsums = np.zeros([len(templates),
                        len(stream[0]) - len(templates[0][0]) + 1])
    for seed_id in seed_ids:
        tr_cc, tr_chans = time_multi_normxcorr(
            template_dict[seed_id], stream_dict[seed_id], pad_dict[seed_id],
            True)
        cccsums = np.sum([cccsums, tr_cc], axis=0)
        no_chans += tr_chans.astype(np.int)
        for chan, state in zip(chans, tr_chans):
            if state:
                chan.append((seed_id.split('.')[1],
                             seed_id.split('.')[-1].split('_')[0]))
    return cccsums, no_chans, chans


@fftw_normxcorr.register('stream_xcorr')
@fftw_normxcorr.register('multithread')
@fftw_normxcorr.register('concurrent')
def _fftw_stream_xcorr(templates, stream, *args, **kwargs):
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
    # number of threads:
    #   default to using inner threads
    #   if `cores` or `cores_outer` passed in then use that
    #   else if OMP_NUM_THREADS set use that
    #   otherwise use all available
    num_cores_inner = kwargs.get('cores')
    num_cores_outer = kwargs.get('cores_outer')
    if num_cores_inner is None and num_cores_outer is None:
        num_cores_inner = int(os.getenv("OMP_NUM_THREADS", cpu_count()))
        num_cores_outer = 1
    elif num_cores_inner is not None and num_cores_outer is None:
        num_cores_outer = 1
    elif num_cores_outer is not None and num_cores_inner is None:
        num_cores_inner = 1

    chans = [[] for _i in range(len(templates))]
    array_dict_tuple = _get_array_dicts(templates, stream)
    stream_dict, template_dict, pad_dict, seed_ids = array_dict_tuple
    assert set(seed_ids)
    cccsums, tr_chans = fftw_multi_normxcorr(
        template_array=template_dict, stream_array=stream_dict,
        pad_array=pad_dict, seed_ids=seed_ids, cores_inner=num_cores_inner,
        cores_outer=num_cores_outer)
    no_chans = np.sum(np.array(tr_chans).astype(np.int), axis=0)
    for seed_id, tr_chan in zip(seed_ids, tr_chans):
        for chan, state in zip(chans, tr_chan):
            if state:
                chan.append((seed_id.split('.')[1],
                             seed_id.split('.')[-1].split('_')[0]))
    return cccsums, no_chans, chans


def fftw_multi_normxcorr(template_array, stream_array, pad_array, seed_ids,
                         cores_inner, cores_outer):
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

    rtype: np.ndarray, list
    :return: 3D Array of cross-correlations and list of used channels.
    """
    utilslib = _load_cdll('libutils')

    utilslib.multi_normxcorr_fftw.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long, ctypes.c_long, ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=np.intc,
                               flags=native_str('C_CONTIGUOUS')),
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
    fft_len = next_fast_len(template_len + image_len - 1)
    template_array = np.ascontiguousarray([template_array[x]
                                           for x in seed_ids],
                                          dtype=np.float32)
    stream_array = np.ascontiguousarray([stream_array[x] for x in seed_ids],
                                        dtype=np.float32)
    cccs = np.zeros((n_templates, image_len - template_len + 1),
                    np.float32)
    used_chans_np = np.ascontiguousarray(used_chans, dtype=np.intc)
    pad_array_np = np.ascontiguousarray([pad_array[seed_id]
                                         for seed_id in seed_ids],
                                        dtype=np.intc)

    # call C function
    ret = utilslib.multi_normxcorr_fftw(
        template_array, n_templates, template_len, n_channels, stream_array,
        image_len, cccs, fft_len, used_chans_np, pad_array_np, cores_outer,
        cores_inner)
    if ret < 0:
        raise MemoryError("Memory allocation failed in correlation C-code")
    elif ret not in [0, 999]:
        print('Error in C code (possible normalisation error)')
        print('Maximum cccs %f at %s' %
              (cccs.max(), np.unravel_index(cccs.argmax(), cccs.shape)))
        print('Minimum cccs %f at %s' %
              (cccs.min(), np.unravel_index(cccs.argmin(), cccs.shape)))
        raise CorrelationError("Internal correlation error")
    elif ret == 999:
        msg = ("Some correlations not computed, are there "
               "zeros in data? If not, consider increasing gain.")
        print(msg)

    return cccs, used_chans


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
    return getattr(func, concur)


# --------------------------- stream prep functions


def _get_array_dicts(templates, stream, copy_streams=True):
    """ prepare templates and stream, return dicts """
    # Do some reshaping
    # init empty structures for data storage
    template_dict = {}
    stream_dict = {}
    pad_dict = {}
    t_starts = []

    stream.sort(['network', 'station', 'location', 'channel'])
    for template in templates:
        template.sort(['network', 'station', 'location', 'channel'])
        t_starts.append(min([tr.stats.starttime for tr in template]))
    # get seed ids, make sure these are collected on sorted streams
    seed_ids = [tr.id + '_' + str(i) for i, tr in enumerate(templates[0])]
    # pull common channels out of streams and templates and put in dicts
    for i, seed_id in enumerate(seed_ids):
        temps_with_seed = [template[i].data for template in templates]
        t_ar = np.array(temps_with_seed).astype(np.float32)
        template_dict.update({seed_id: t_ar})
        stream_dict.update(
            {seed_id: stream.select(
                id=seed_id.split('_')[0])[0].data.astype(np.float32)})
        pad_list = [
            int(round(template[i].stats.sampling_rate *
                      (template[i].stats.starttime - t_starts[j])))
            for j, template in zip(range(len(templates)), templates)]
        pad_dict.update({seed_id: pad_list})

    return stream_dict, template_dict, pad_dict, seed_ids


# a dict of built in xcorr functions, used to distinguish from user-defined
XCORR_FUNCS_ORIGINAL = copy.copy(XCOR_FUNCS)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
