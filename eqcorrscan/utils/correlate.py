"""
Correlation functions for multi-channel cross-correlation of seismic data.

Various routines used mostly for testing, including links to a compiled
routine using FFTW, a Numpy fft routine which uses bottleneck for normalisation
and a compiled time-domain routine. These have varying levels of efficiency,
both in terms of overall speed, and in memory usage.  The time-domain is the
most memory efficient but slowest routine (although fastest for small cases of
less than a few hundred correlations), the Numpy routine is fast, but memory
inefficient due to a need to store large double-precision arrays for
normalisation.  The fftw compiled routine is fastest and more memory efficient
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

import numpy as np
import ctypes
from future.utils import native_str

from multiprocessing import Pool
from scipy.fftpack.helper import next_fast_len
from eqcorrscan.utils.libnames import _load_cdll


def numpy_normxcorr(templates, stream, pads):
    """
    Compute the normalized cross-correlation of multiple templates with data.

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
    for i in range(len(pads)):
        res[i] = np.append(res[i], np.zeros(pads[i]))[pads[i]:]
    return res.astype(np.float32), used_chans


def multichannel_normxcorr(templates, stream, cores=1, time_domain=False,
                           openmp=False):
    """
    Cross-correlate multiple channels either in parallel or not

    :type templates: list
    :param templates:
        A list of templates, where each one should be an obspy.Stream object
        containing multiple traces of seismic data and the relevant header
        information.
    :type stream: obspy.core.stream.Stream
    :param stream:
        A single Stream object to be correlated with the templates.
    :type cores: int
    :param cores:
        Number of processes to use, if set to None, no Python multiprocessing
        will be done.
    :type cores: int
    :param cores: Number of cores to loop over
    :type time_domain: bool
    :param time_domain:
        Whether to compute in the time-domain using the compiled openMP
        parallel cross-correlation routine.
    :type openmp: bool
    :param openmp: Whether to use the openmp, compiled threaded loop or not.

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

    .. Note::
        Each template must contain the same channels as every other template,
        the stream must also contain the same channels (note that if there
        are duplicate channels in the template you do not need duplicate
        channels in the stream).
    """
    no_chans = np.zeros(len(templates))
    chans = [[] for _i in range(len(templates))]
    # Do some reshaping
    stream.sort(['network', 'station', 'location', 'channel'])
    t_starts = []
    for template in templates:
        template.sort(['network', 'station', 'location', 'channel'])
        t_starts.append(min([tr.stats.starttime for tr in template]))
    seed_ids = [tr.id + '_' + str(i) for i, tr in enumerate(templates[0])]
    template_array = {}
    stream_array = {}
    pad_array = {}
    for i, seed_id in enumerate(seed_ids):
        t_ar = np.array([template[i].data
                         for template in templates]).astype(np.float32)
        template_array.update({seed_id: t_ar})
        stream_array.update(
            {seed_id: stream.select(
                id=seed_id.split('_')[0])[0].data.astype(np.float32)})
        pad_list = [
            int(round(template[i].stats.sampling_rate *
                      (template[i].stats.starttime - t_starts[j])))
            for j, template in zip(range(len(templates)), templates)]
        pad_array.update({seed_id: pad_list})
    if cores is None and not openmp:
        cccsums = np.zeros([len(templates),
                            len(stream[0]) - len(templates[0][0]) + 1])
        for seed_id in seed_ids:
            if time_domain:
                tr_xcorrs, tr_chans = time_multi_normxcorr(
                    templates=template_array[seed_id],
                    stream=stream_array[seed_id], pads=pad_array[seed_id])
            else:
                tr_xcorrs, tr_chans = fftw_normxcorr(
                    templates=template_array[seed_id],
                    stream=stream_array[seed_id], pads=pad_array[seed_id],
                    threaded=True)
            cccsums = np.sum([cccsums, tr_xcorrs], axis=0)
            no_chans += tr_chans.astype(np.int)
            for chan, state in zip(chans, tr_chans):
                if state:
                    chan.append((seed_id.split('.')[1],
                                 seed_id.split('.')[-1].split('_')[0]))
    elif not openmp:
        pool = Pool(processes=cores)
        if time_domain:
            results = [pool.apply_async(time_multi_normxcorr, (
                template_array[seed_id], stream_array[seed_id],
                pad_array[seed_id])) for seed_id in seed_ids]
        else:
            results = [pool.apply_async(fftw_normxcorr, (
                template_array[seed_id], stream_array[seed_id],
                pad_array[seed_id], False)) for seed_id in seed_ids]
        pool.close()
        results = [p.get() for p in results]
        xcorrs = [p[0] for p in results]
        tr_chans = np.array([p[1] for p in results])
        pool.join()
        cccsums = np.sum(xcorrs, axis=0)
        no_chans = np.sum(tr_chans.astype(np.int), axis=0)
        for seed_id, tr_chan in zip(seed_ids, tr_chans):
            for chan, state in zip(chans, tr_chan):
                if state:
                    chan.append((seed_id.split('.')[1],
                                 seed_id.split('.')[-1].split('_')[0]))
    else:
        xcorrs, tr_chans = fftw_multi_normxcorr(
            template_array=template_array, stream_array=stream_array,
            pad_array=pad_array, seed_ids=seed_ids)
        cccsums = np.sum(xcorrs, axis=0)
        no_chans = np.sum(np.array(tr_chans).astype(np.int), axis=0)
        for seed_id, tr_chan in zip(seed_ids, tr_chans):
            for chan, state in zip(chans, tr_chan):
                if state:
                    chan.append((seed_id.split('.')[1],
                                 seed_id.split('.')[-1].split('_')[0]))
    return cccsums, no_chans, chans


def time_multi_normxcorr(templates, stream, pads):
    """
    Compute cross-correlations in the time-domain using C routine.

    :param templates: 2D Array of templates
    :type templates: np.ndarray
    :param stream: 1D array of continuous data
    :type stream: np.ndarray
    :param pads: List of ints of pad lengths in the same order as templates
    :type pads: list

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    used_chans = ~np.isnan(templates).any(axis=1)

    utilslib = _load_cdll('libutils')

    utilslib.multi_normxcorr_time.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.multi_normxcorr_time.restype = ctypes.c_int

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
    utilslib.multi_normxcorr_time(
        t_array, template_len, n_templates,
        np.ascontiguousarray(stream, np.float32), image_len, ccc)
    ccc[np.isnan(ccc)] = 0.0
    ccc = ccc.reshape((n_templates, image_len - template_len + 1))
    for i in range(len(pads)):
        ccc[i] = np.append(ccc[i], np.zeros(pads[i]))[pads[i]:]
    templates += templates_means
    stream += stream_mean
    return ccc, used_chans


def fftw_normxcorr(templates, stream, pads, threaded=True):
    """
    Normalised cross-correlation using the fftw library.

    Internally this function used double precision numbers, which is definitely
    required for seismic data. Cross-correlations are computed as the
    inverse fft of the dot product of the ffts of the stream and the reversed,
    normalised, templates.  The cross-correlation is then normalised using the
    running mean and standard deviation (not using the N-1 correction) of the
    stream and the sums of the normalised templates.

    This python fucntion wraps the C-library written by C. Chamberlain for this
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
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
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
    fftshape = next_fast_len(template_length + stream_length - 1)
    # # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))

    norm = np.nan_to_num(norm)
    ccc = np.empty((n_templates, stream_length - template_length + 1),
                   np.float32).flatten(order='C')
    ret = func(
        np.ascontiguousarray(norm.flatten(order='C'), np.float32),
        template_length, n_templates,
        np.ascontiguousarray(stream, np.float32), stream_length,
        np.ascontiguousarray(ccc, np.float32), fftshape)
    if ret != 0:
        print(ret)
        raise MemoryError()
    ccc = ccc.reshape((n_templates, stream_length - template_length + 1))
    for i in range(n_templates):
        if not used_chans[i]:
            ccc[i] = np.zeros(stream_length - template_length + 1)
    ccc[np.isnan(ccc)] = 0.0
    if np.any(np.abs(ccc) > 1.01):
        print('Normalisation error in C code')
        print(ccc.max())
        print(ccc.min())
        raise MemoryError()
    ccc[ccc > 1.0] = 1.0
    ccc[ccc < -1.0] = -1.0
    for i in range(len(pads)):
        ccc[i] = np.append(ccc[i], np.zeros(pads[i]))[pads[i]:]
    return ccc, used_chans


def fftw_multi_normxcorr(template_array, stream_array, pad_array, seed_ids):
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
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int]
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
    '''
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
    template_array = np.array(list(template_array.values())).flatten(order='C')
    stream_array = np.array(list(stream_array.values())).flatten(order='C')
    cccs = np.empty((n_channels, n_templates, image_len - template_len + 1),
                    np.float32).flatten(order='C')
    ret = utilslib.multi_normxcorr_fftw(
        template_array, n_templates, template_len, n_channels, stream_array,
        image_len, cccs, fft_len)
    if ret != 0:
        raise MemoryError()
    cccs = cccs.reshape((n_channels, n_templates,
                         image_len - template_len + 1))
    for j in range(n_channels):
        for i in range(n_templates):
            if not used_chans[j][i]:
                cccs[j][i] = np.zeros(image_len - template_len + 1)
    cccs[np.isnan(cccs)] = 0.0
    if np.any(np.abs(cccs) > 1.01):
        print('Normalisation error in C code')
        print(cccs.max())
        print(cccs.min())
        raise MemoryError()
    cccs[cccs > 1.0] = 1.0
    cccs[cccs < -1.0] = -1.0
    for j, seed_id in enumerate(seed_ids):
        for i in range(len(pad_array[seed_id])):
            cccs[j][i] = np.append(
                cccs[j][i],
                np.zeros(pad_array[seed_id][i]))[pad_array[seed_id][i]:]
    cccs = cccs.reshape(n_channels, n_templates, image_len - template_len + 1)
    return cccs, used_chans


if __name__ == '__main__':
    import doctest
    doctest.testmod()
