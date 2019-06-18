"""
Functions for despiking seismic data.

.. warning:: Despike functions are in beta, they do not work that well.

.. todo:: Deconvolve spike to remove it, find peaks in the f-domain.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import logging
import numpy as np

from multiprocessing import Pool, cpu_count
from obspy import Trace

from eqcorrscan.utils.timer import Timer
from eqcorrscan.utils.correlate import get_array_xcorr
from eqcorrscan.utils.findpeaks import find_peaks2_short


Logger = logging.getLogger(__name__)


def median_filter(tr, multiplier=10, windowlength=0.5,
                  interp_len=0.05):
    """
    Filter out spikes in data above a multiple of MAD of the data.

    Currently only has the ability to replaces spikes with linear
    interpolation.  In the future we would aim to fill the gap with something
    more appropriate.  Works in-place on data.

    :type tr: obspy.core.trace.Trace
    :param tr: trace to despike
    :type multiplier: float
    :param multiplier:
        median absolute deviation multiplier to find spikes above.
    :type windowlength: float
    :param windowlength: Length of window to look for spikes in in seconds.
    :type interp_len: float
    :param interp_len: Length in seconds to interpolate around spikes.

    :returns: :class:`obspy.core.trace.Trace`

    .. warning::
        Not particularly effective, and may remove earthquake signals, use with
        caution.
    """
    num_cores = cpu_count()
    # Note - might be worth finding spikes in filtered data
    filt = tr.copy()
    filt.detrend('linear')
    try:
        filt.filter('bandpass', freqmin=10.0,
                    freqmax=(tr.stats.sampling_rate / 2) - 1)
    except Exception as e:
        Logger.error("Could not filter due to error: {0}".format(e))
    data = filt.data
    del filt
    # Loop through windows
    _windowlength = int(windowlength * tr.stats.sampling_rate)
    _interp_len = int(interp_len * tr.stats.sampling_rate)
    peaks = []
    with Timer() as t:
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(_median_window,
                                    args=(data[chunk * _windowlength:
                                               (chunk + 1) * _windowlength],
                                          chunk * _windowlength, multiplier,
                                          tr.stats.starttime + windowlength,
                                          tr.stats.sampling_rate))
                   for chunk in range(int(len(data) / _windowlength))]
        pool.close()
        for p in results:
            peaks += p.get()
        pool.join()
        for peak in peaks:
            tr.data = _interp_gap(tr.data, peak[1], _interp_len)
    Logger.debug("Despiking took: %s s" % t.secs)
    return tr


def _median_window(window, window_start, multiplier, starttime, sampling_rate):
    """
    Internal function to aid parallel processing

    :type window: numpy.ndarry
    :param window: Data to look for peaks in.
    :type window_start: int
    :param window_start: Index of window start point in larger array, used \
        for peak indexing.
    :type multiplier: float
    :param multiplier: Multiple of MAD to use as threshold
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Starttime of window, used in debug plotting.
    :type sampling_rate: float
    :param sampling_rate in Hz, used for debug plotting

    :returns: peaks
    :rtype: list
    """
    MAD = np.median(np.abs(window))
    thresh = multiplier * MAD
    Logger.debug(
        'Threshold for window is: ' + str(thresh) + '\nMedian is: ' +
        str(MAD) + '\nMax is: ' + str(np.max(window)))
    peaks = find_peaks2_short(arr=window, thresh=thresh, trig_int=5)
    if peaks:
        peaks = [(peak[0], peak[1] + window_start) for peak in peaks]
    else:
        peaks = []
    return peaks


def _interp_gap(data, peak_loc, interp_len):
    """
    Internal function for filling gap with linear interpolation

    :type data: numpy.ndarray
    :param data: data to remove peak in
    :type peak_loc: int
    :param peak_loc: peak location position
    :type interp_len: int
    :param interp_len: window to interpolate

    :returns: Trace works in-place
    :rtype: :class:`obspy.core.trace.Trace`
    """
    start_loc = peak_loc - int(0.5 * interp_len)
    end_loc = peak_loc + int(0.5 * interp_len)
    if start_loc < 0:
        start_loc = 0
    if end_loc > len(data) - 1:
        end_loc = len(data) - 1
    fill = np.linspace(data[start_loc], data[end_loc], end_loc - start_loc)
    data[start_loc:end_loc] = fill
    return data


def template_remove(tr, template, cc_thresh, windowlength, interp_len):
    """
    Looks for instances of template in the trace and removes the matches.

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to remove spikes from.
    :type template: osbpy.core.trace.Trace
    :param template: Spike template to look for in data.
    :type cc_thresh: float
    :param cc_thresh: Cross-correlation threshold (-1 - 1).
    :type windowlength: float
    :param windowlength: Length of window to look for spikes in in seconds.
    :type interp_len: float
    :param interp_len: Window length to remove and fill in seconds.

    :returns: tr, works in place.
    :rtype: :class:`obspy.core.trace.Trace`
    """
    _interp_len = int(tr.stats.sampling_rate * interp_len)
    if _interp_len < len(template.data):
        Logger.warning('Interp_len is less than the length of the template, '
                       'will used the length of the template!')
        _interp_len = len(template.data)
    if isinstance(template, Trace):
        template = np.array([template.data])
    with Timer() as t:
        normxcorr = get_array_xcorr("fftw")
        cc, _ = normxcorr(stream=tr.data.astype(np.float32),
                          templates=template.astype(np.float32), pads=[0])
        peaks = find_peaks2_short(
            arr=cc.flatten(), thresh=cc_thresh,
            trig_int=windowlength * tr.stats.sampling_rate)
        for peak in peaks:
            tr.data = _interp_gap(
                data=tr.data, peak_loc=peak[1] + int(0.5 * _interp_len),
                interp_len=_interp_len)
    Logger.info("Despiking took: {0:.4f} s".format(t.secs))
    return tr


if __name__ == '__main__':
    import doctest
    doctest.testmod()
