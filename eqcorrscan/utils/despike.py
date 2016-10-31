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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool, cpu_count

from eqcorrscan.utils.timer import Timer


def median_filter(tr, multiplier=10, windowlength=0.5,
                  interp_len=0.05, debug=0):
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
    :type debug: int
    :param debug: Debug output level between 0 and 5, higher is more output.

    :returns: :class:`obspy.core.trace.Trace`

    .. warning::
        Not particularly effective, and may remove earthquake signals, use with
        caution.
    """
    num_cores = cpu_count()
    if debug >= 1:
        data_in = tr.copy()
    # Note - might be worth finding spikes in filtered data
    filt = tr.copy()
    filt.detrend('linear')
    filt.filter('bandpass', freqmin=10.0,
                freqmax=(tr.stats.sampling_rate / 2) - 1)
    data = filt.data
    del(filt)
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
                                          tr.stats.sampling_rate,
                                          debug))
                   for chunk in range(int(len(data) / _windowlength))]
        pool.close()
        for p in results:
            peaks += p.get()
        pool.join()
        for peak in peaks:
            tr.data = _interp_gap(tr.data, peak[1], _interp_len)
    print("Despiking took: %s s" % t.secs)
    if debug >= 1:
        plt.plot(data_in.data, 'r', label='raw')
        plt.plot(tr.data, 'k', label='despiked')
        plt.legend()
        plt.show()
    return tr


def _median_window(window, window_start, multiplier, starttime, sampling_rate,
                   debug=0):
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
    :type debug: int
    :param debug: debug level, if want plots, >= 4.

    :returns: peaks
    :rtype: list
    """
    from eqcorrscan.utils.findpeaks import find_peaks2_short
    from eqcorrscan.utils.plotting import peaks_plot

    MAD = np.median(np.abs(window))
    thresh = multiplier * MAD
    if debug >= 2:
        print('Threshold for window is: ' + str(thresh) +
              '\nMedian is: ' + str(MAD) +
              '\nMax is: ' + str(np.max(window)))
    peaks = find_peaks2_short(arr=window,
                              thresh=thresh, trig_int=5, debug=0)
    if debug >= 4 and peaks:
        peaks_plot(window, starttime, sampling_rate,
                   save=False, peaks=peaks)
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


def template_remove(tr, template, cc_thresh, windowlength,
                    interp_len, debug=0):
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
    :type debug: int
    :param debug: Debug level.

    :returns: tr, works in place.
    :rtype: :class:`obspy.core.trace.Trace`
    """
    from eqcorrscan.core.match_filter import normxcorr2
    from eqcorrscan.utils.findpeaks import find_peaks2_short
    from obspy import Trace
    from eqcorrscan.utils.timer import Timer
    import matplotlib.pyplot as plt
    import warnings

    data_in = tr.copy()
    _interp_len = int(tr.stats.sampling_rate * interp_len)
    if _interp_len < len(template.data):
        warnings.warn('Interp_len is less than the length of the template,'
                      'will used the length of the template!')
        _interp_len = len(template.data)
    if isinstance(template, Trace):
        template = template.data
    with Timer() as t:
        cc = normxcorr2(tr.data.astype(np.float32),
                        template.astype(np.float32))
        if debug > 3:
            plt.plot(cc.flatten(), 'k', label='cross-correlation')
            plt.legend()
            plt.show()
        peaks = find_peaks2_short(arr=cc.flatten(), thresh=cc_thresh,
                                  trig_int=windowlength * tr.stats.
                                  sampling_rate)
        for peak in peaks:
            tr.data = _interp_gap(data=tr.data,
                                  peak_loc=peak[1] + int(0.5 * _interp_len),
                                  interp_len=_interp_len)
    print("Despiking took: %s s" % t.secs)
    if debug > 2:
        plt.plot(data_in.data, 'r', label='raw')
        plt.plot(tr.data, 'k', label='despiked')
        plt.legend()
        plt.show()
    return tr


if __name__ == '__main__':
    import doctest
    doctest.testmod()
