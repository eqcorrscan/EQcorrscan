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
from scipy.fftpack.helper import next_fast_len
from eqcorrscan.utils.libnames import _load_cdll


def scipy_normxcorr(templates, stream, pads):
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
        Whether to compute in the time-domain using the compiled openMP
        parallel cross-correlation routine.

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
    if cores is None:
        cccsums = np.zeros([len(templates),
                            len(stream[0]) - len(templates[0][0]) + 1])
        for seed_id in seed_ids:
            if time_domain:
                tr_xcorrs, tr_chans = time_multi_normxcorr(
                    templates=template_array[seed_id],
                    stream=stream_array[seed_id], pads=pad_array[seed_id])
            else:
                tr_xcorrs, tr_chans = fftw_xcorr(
                    templates=template_array[seed_id],
                    stream=stream_array[seed_id], pads=pad_array[seed_id])
            cccsums = np.sum([cccsums, tr_xcorrs], axis=0)
            no_chans += tr_chans.astype(np.int)
            for chan, state in zip(chans, tr_chans):
                if state:
                    chan.append((seed_id.split('.')[1],
                                 seed_id.split('.')[-1].split('_')[0]))
    else:
        pool = Pool(processes=cores)
        if time_domain:
            results = [pool.apply_async(time_multi_normxcorr, (
                template_array[seed_id], stream_array[seed_id],
                pad_array[seed_id])) for seed_id in seed_ids]
        else:
            results = [pool.apply_async(fftw_xcorr, (
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


def fftw_xcorr(templates, stream, pads):
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

    :return: np.ndarray of cross-correlations
    :return: np.ndarray channels used
    """
    from future.utils import native_str

    utilslib = _load_cdll('libutils')

    utilslib.normxcorr_fftw_1d.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int]
    utilslib.normxcorr_fftw_1d.restype = ctypes.c_int
    # Generate a template mask
    used_chans = ~np.isnan(templates).any(axis=1)
    template_length = templates.shape[1]
    stream_length = len(stream)
    n_templates = templates.shape[0]
    fftshape = next_fast_len(template_length + stream_length - 1)
    # # Normalize and flip the templates
    norm = ((templates - templates.mean(axis=-1, keepdims=True)) / (
        templates.std(axis=-1, keepdims=True) * template_length))

    ccc = np.empty((n_templates, stream_length - template_length + 1),
                   np.float32)
    for i in range(n_templates):
        if np.all(np.isnan(norm[i])):
            ccc[i] = np.zeros(stream_length - template_length + 1)
        else:
            ret = utilslib.normxcorr_fftw_1d(
                np.ascontiguousarray(norm[i], np.float32), template_length,
                np.ascontiguousarray(stream, np.float32), stream_length,
                np.ascontiguousarray(ccc[i], np.float32), fftshape)
            if ret != 0:
                raise MemoryError()
    ccc = ccc.reshape((n_templates, stream_length - template_length + 1))
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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
