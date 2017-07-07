"""
Correlation functions for multi-channel cross-correlation of seismic data.

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

from multiprocessing import Pool
from scipy.signal.signaltools import _centered
from scipy.fftpack.helper import next_fast_len
from eqcorrscan.utils.normalise import multi_norm
from eqcorrscan.utils.libnames import _load_cdll
from eqcorrscan.core.match_filter import MatchFilterError


def multi_normxcorr(templates, stream, pads):
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
    # Generate a template mask
    used_chans = ~np.isnan(templates).any(axis=1)
    template_length = templates.shape[1]
    stream_length = len(stream)
    fftshape = next_fast_len(template_length + stream_length - 1)
    # # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))
    # # CPU bound
    norm_sum = norm.sum(axis=-1, keepdims=True)
    # CPU bound
    res = np.fft.irfft(
        np.fft.rfft(np.flip(norm, axis=-1), fftshape, axis=-1) *
        np.fft.rfft(stream, fftshape),
        fftshape)[:, 0:template_length + stream_length - 1]
    del norm
    res = res.astype(np.float32)
    # CPU bound
    res = multi_norm(_centered(res, stream_length - template_length + 1),
                     stream, norm_sum, template_length)
    for i in range(len(pads)):
        res[i] = np.append(res[i], np.zeros(pads[i]))[pads[i]:]
    return res, used_chans


def multichannel_xcorr(templates, stream, cores=1, time_domain=False):
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
        Number of processed to use, if set to None, and dask==False, no
        multiprocessing will be done.
    :type cores: int
    :param cores: Number of cores to loop over
    :type time_domain: bool
    :param time_domain:
        Whether to compute in the time-domain using the compiled openMP parallel
        cross-correlation routine.

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
    if cores is None and not time_domain:
        cccsums = np.zeros([len(templates),
                            len(stream[0]) - len(templates[0][0]) + 1])
        for seed_id in seed_ids:
            # tr_xcorrs, tr_chans = multi_normxcorr(
            tr_xcorrs, tr_chans = fftw_compiled_xcorr(
                templates=template_array[seed_id],
                stream=stream_array[seed_id], pads=pad_array[seed_id])
            cccsums = np.sum([cccsums, tr_xcorrs], axis=0)
            no_chans += tr_chans.astype(np.int)
            for chan, state in zip(chans, tr_chans):
                if state:
                    chan.append((seed_id.split('.')[1],
                                 seed_id.split('.')[-1].split('_')[0]))
    elif not time_domain:
        pool = Pool(processes=cores)
        # results = [pool.apply_async(multi_normxcorr, (
        results = [pool.apply_async(fftw_compiled_xcorr, (
            template_array[seed_id], stream_array[seed_id],
            pad_array[seed_id])) for seed_id in seed_ids]
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
        cccsums = np.zeros([len(templates),
                            len(stream[0]) - len(templates[0][0]) + 1])
        for seed_id in seed_ids:
            tr_xcorrs, tr_chans = time_multi_normxcorr(
                templates=template_array[seed_id],
                stream=stream_array[seed_id], pads=pad_array[seed_id])
            cccsums = np.sum([cccsums, tr_xcorrs], axis=0)
            no_chans += tr_chans.astype(np.int)
            for chan, state in zip(chans, tr_chans):
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
    from future.utils import native_str

    used_chans = ~np.isnan(templates).any(axis=1)

    utilslib = _load_cdll('libutils')

    utilslib.multi_corr.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.multi_corr.restype = ctypes.c_int

    template_len = templates.shape[1]
    n_templates = templates.shape[0]
    image_len = stream.shape[0]
    ccc = np.ascontiguousarray(
        np.empty((image_len - template_len + 1) * n_templates), np.float32)
    t_array = np.ascontiguousarray(templates.flatten(), np.float32)
    utilslib.multi_corr(t_array, template_len, n_templates,
                        np.ascontiguousarray(stream, np.float32), image_len,
                        ccc)
    ccc[np.isnan(ccc)] = 0.0
    ccc = ccc.reshape((n_templates, image_len - template_len + 1))
    for i in range(len(pads)):
        ccc[i] = np.append(ccc[i], np.zeros(pads[i]))[pads[i]:]
    return ccc, used_chans


def fftw_compiled_xcorr(templates, stream, pads):
    """

    :param templates:
    :param stream:
    :param pads:
    :return:
    """
    from future.utils import native_str

    utilslib = _load_cdll('libutils')

    utilslib.normxcorr_fftw_loop.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int]
    utilslib.normxcorr_fftw_loop.restype = ctypes.c_int
    # Generate a template mask
    used_chans = ~np.isnan(templates).any(axis=1)
    template_length = templates.shape[1]
    stream_length = len(stream)
    n_templates = templates.shape[0]
    fftshape = next_fast_len(template_length + stream_length - 1)
    # # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))
    norm = np.ascontiguousarray(norm.flatten(), np.float32)
    ccc = np.empty((n_templates, stream_length - template_length + 1),
                   np.float32)
    ccc = np.ascontiguousarray(ccc.flatten())
    ret = utilslib.normxcorr_fftw_loop(
        norm, template_length, np.ascontiguousarray(stream, np.float32),
        stream_length, ccc, fftshape, n_templates)
    if ret:
        raise MemoryError()
    ccc = ccc.reshape((n_templates, stream_length - template_length + 1))
    ccc[np.isnan(ccc)] = 0.0
    if ret != 0:
        raise MemoryError()
    if np.any(ccc > 1.001):
        raise MatchFilterError('Normalisation error in C code')
    ccc[ccc > 1.0] = 1.0
    for i in range(len(pads)):
        ccc[i] = np.append(ccc[i], np.zeros(pads[i]))[pads[i]:]
    return ccc, used_chans


if __name__ == '__main__':
    import doctest
    doctest.testmod()