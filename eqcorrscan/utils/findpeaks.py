"""Functions to find peaks in data above a certain threshold.

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

import ctypes
import random
import numpy as np

from obspy import UTCDateTime
from scipy import ndimage
from multiprocessing import Pool
from future.utils import native_str
from itertools import compress

from eqcorrscan.utils.correlate import pool_boy
from eqcorrscan.utils.libnames import _load_cdll
from eqcorrscan.utils.debug_log import debug_print


def is_prime(number):
    """
    Function to test primality of a number. Function lifted from online
    resource:
        http://www.codeproject.com/Articles/691200/Primality-test-algorithms-Prime-test-The-fastest-w

    This function is distributed under a separate licence:
        This article, along with any associated source code and files, is \
        licensed under The Code Project Open License (CPOL)

    :type number: int
    :param number: Integer to test for primality

    :returns: bool

    >>> is_prime(4)
    False
    >>> is_prime(3)
    True
    """
    ''' if number != 1 '''
    if number > 1:
        ''' repeat the test few times '''
        for time in range(3):
            ''' Draw a RANDOM number in range of number ( Z_number )  '''
            randomNumber = random.randint(2, number - 1)
            ''' Test if a^(n-1) = 1 mod n '''
            if pow(randomNumber, number - 1, number) != 1:
                return False
        return True
    else:
        ''' case number == 1 '''
        return False


def find_peaks2_short(arr, thresh, trig_int, debug=0, starttime=False,
                      samp_rate=1.0, full_peaks=False):
    """
    Determine peaks in an array of data above a certain threshold.

    Uses a mask to remove data below threshold and finds peaks in what is left.

    :type arr: numpy.ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh:
        The threshold below which will be considered noise and peaks will
        not be found in.
    :type trig_int: int
    :param trig_int:
        The minimum difference in samples between triggers, if multiple
        peaks within this window this code will find the highest.
    :type debug: int
    :param debug: Optional, debug level 0-5
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Starttime for plotting, only used if debug > 2.
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz, only used for plotting if debug > 2.
    :type full_peaks: bool
    :param full_peaks:
        If True, will remove the issue eluded to below, by declustering within
        data-sections above the threshold, rather than just taking the peak
        within that section. This will take more time. This defaults to True
        for match_filter.

    :return: peaks: Lists of tuples of peak values and locations.
    :rtype: list


    >>> import numpy as np
    >>> arr = np.random.randn(100)
    >>> threshold = 10
    >>> arr[40] = 20
    >>> arr[60] = 100
    >>> find_peaks2_short(arr, threshold, 3)
    [(20.0, 40), (100.0, 60)]

    .. note::
        peak-finding is optimised for zero-mean cross-correlation data where
        fluctuations are frequent.  Because of this, in certain cases some
        peaks may be missed if the trig_int is short and the threshold is low.
        Consider the following case:

        >>> arr = np.array([1, .2, .2, .2, .2, 1, .2, .2, .2, .2, 1])
        >>> find_peaks2_short(arr, thresh=.2, trig_int=3)
        [(1.0, 0)]

        Whereas you would expect the following:

        >>> arr = np.array([1, .2, .2, .2, .2, 1, .2, .2, .2, .2, 1])
        >>> find_peaks2_short(arr, thresh=.2, trig_int=3, full_peaks=True)
        [(1.0, 0), (1.0, 5), (1.0, 10)]

        This is rare and unlikely to happen for correlation cases, where
        trigger intervals are usually large and thresholds high.

    """
    if not starttime:
        starttime = UTCDateTime(0)
    # Set everything below the threshold to zero
    image = np.copy(arr)
    image = np.abs(image)
    debug_print("Threshold: {0}\tMax: {1}".format(thresh, max(image)),
                2, debug)
    image[image < thresh] = 0
    if len(image[image > thresh]) == 0:
        debug_print("No values over threshold {0}".format(thresh), 0, debug)
        return []
    debug_print('Found {0} samples above the threshold'.format(
        len(image[image > thresh])), 0, debug)
    initial_peaks = []
    # Find the peaks
    labeled_image, number_of_objects = ndimage.label(image)
    peak_slices = ndimage.find_objects(labeled_image)
    for peak_slice in peak_slices:
        window = arr[peak_slice[0].start: peak_slice[0].stop]
        if peak_slice[0].stop - peak_slice[0].start > trig_int and full_peaks:
            peaks = decluster(
                peaks=window, trig_int=trig_int,
                index=np.arange(peak_slice[0].start, peak_slice[0].stop))
        else:
            peaks = [(window[np.argmax(abs(window))],
                      int(peak_slice[0].start + np.argmax(abs(window))))]
        initial_peaks.extend(peaks)
    peaks = decluster(peaks=np.array(list(zip(*initial_peaks))[0]),
                      index=np.array(list(zip(*initial_peaks))[1]),
                      trig_int=trig_int)
    if initial_peaks:
        if debug >= 3:
            from eqcorrscan.utils import plotting
            _fname = ''.join([
                'peaks_', starttime.datetime.strftime('%Y-%m-%d'), '.pdf'])
            plotting.peaks_plot(
                data=image, starttime=starttime, samp_rate=samp_rate,
                save=True, peaks=peaks, savefile=_fname)
        peaks = sorted(peaks, key=lambda time: time[1], reverse=False)
        return peaks
    else:
        print('No peaks for you!')
        return []


def multi_find_peaks(arr, thresh, trig_int, debug=0, starttime=False,
                     samp_rate=1.0, parallel=True, full_peaks=False,
                     cores=None):
    """
    Wrapper for find-peaks for multiple arrays.

    :type arr: numpy.ndarray
    :param arr: 2-D numpy array is required
    :type thresh: list
    :param thresh:
        The threshold below which will be considered noise and peaks will not
        be found in. One threshold per array.
    :type trig_int: int
    :param trig_int:
        The minimum difference in samples between triggers, if multiple
        peaks within this window this code will find the highest.
    :type debug: int
    :param debug: Optional, debug level 0-5
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Starttime for plotting, only used if debug > 2.
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz, only used for plotting if debug > 2.
    :type parallel: bool
    :param parallel:
        Whether to compute in parallel or not - will use multiprocessing
    :type full_peaks: bool
    :param full_peaks: See `eqcorrscan.utils.findpeaks.find_peaks2_short`
    :type cores: int
    :param cores:
        Maximum number of processes to spin up for parallel peak-finding

    :returns:
        List of list of tuples of (peak, index) in same order as input arrays
    """
    peaks = []
    if not parallel:
        for sub_arr, arr_thresh in zip(arr, thresh):
            peaks.append(find_peaks2_short(
                arr=sub_arr, thresh=arr_thresh, trig_int=trig_int, debug=debug,
                starttime=starttime, samp_rate=samp_rate,
                full_peaks=full_peaks))
    else:
        if cores is None:
            cores = arr.shape[0]
        with pool_boy(Pool=Pool, traces=cores) as pool:
            params = ((sub_arr, arr_thresh, trig_int, debug,
                       False, 1.0, full_peaks)
                      for sub_arr, arr_thresh in zip(arr, thresh))
            results = [pool.apply_async(find_peaks2_short, param)
                       for param in params]
            peaks = [res.get() for res in results]
    return peaks


def decluster(peaks, index, trig_int):
    """
    Decluster peaks based on an enforced minimum separation.
    :type peaks: np.array
    :param peaks: array of peak values
    :type index: np.ndarray
    :param index: locations of peaks
    :type trig_int: int
    :param trig_int: Minimum trigger interval in samples

    :return: list of tuples of (value, sample)
    """
    utilslib = _load_cdll('libutils')

    length = np.int32(len(peaks))
    utilslib.find_peaks.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_float, ctypes.c_float,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.find_peaks.restype = ctypes.c_int
    peaks_sort = sorted(zip(peaks, index),
                        key=lambda amplitude: abs(amplitude[0]),
                        reverse=True)
    arr, inds = zip(*peaks_sort)
    arr = np.ascontiguousarray(arr, dtype=np.float32)
    inds = np.array(inds, dtype=np.float32) / trig_int
    inds = np.ascontiguousarray(inds, dtype=np.float32)
    out = np.zeros(len(arr), dtype=np.uint32)
    ret = utilslib.find_peaks(
        arr, inds, length, 0, np.float32(1), out)
    if ret != 0:
        raise MemoryError("Issue with c-routine, returned %i" % ret)
    peaks_out = list(compress(peaks_sort, out))
    return peaks_out


def coin_trig(peaks, stachans, samp_rate, moveout, min_trig, trig_int):
    """
    Find network coincidence triggers within peaks of detection statistics.

    Useful for finding network detections from sets of detections on individual
    stations.

    :type peaks: list
    :param peaks: List of lists of tuples of (peak, index) for each \
        station-channel.  Index should be in samples.
    :type stachans: list
    :param stachans: List of tuples of (station, channel) in the order of \
        peaks.
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz
    :type moveout: float
    :param moveout: Allowable network moveout in seconds.
    :type min_trig: int
    :param min_trig: Minimum station-channels required to declare a trigger.
    :type trig_int: float
    :param trig_int:
        Minimum allowable time between network triggers in seconds.

    :return:
        List of tuples of (peak, index), for the earliest detected station.
    :rtype: list

    >>> peaks = [[(0.5, 100), (0.3, 800)], [(0.4, 120), (0.7, 850)]]
    >>> triggers = coin_trig(peaks, [('a', 'Z'), ('b', 'Z')], 10, 3, 2, 1)
    >>> print(triggers)
    [(0.45, 100)]
    """
    triggers = []
    for stachan, _peaks in zip(stachans, peaks):
        for peak in _peaks:
            trigger = (peak[1], peak[0], '.'.join(stachan))
            triggers.append(trigger)
    coincidence_triggers = []
    for i, master in enumerate(triggers):
        slaves = triggers[i + 1:]
        coincidence = 1
        trig_time = master[0]
        trig_val = master[1]
        for slave in slaves:
            if abs(slave[0] - master[0]) <= (moveout * samp_rate) and \
               slave[2] != master[2]:
                coincidence += 1
                if slave[0] < master[0]:
                    trig_time = slave[0]
                trig_val += slave[1]
        if coincidence >= min_trig:
            coincidence_triggers.append((trig_val / coincidence,
                                         trig_time))
    # Sort by trigger-value, largest to smallest - remove duplicate detections
    if coincidence_triggers:
        coincidence_triggers.sort(key=lambda tup: tup[0], reverse=True)
        output = [coincidence_triggers[0]]
        for coincidence_trigger in coincidence_triggers[1:]:
            add = True
            for peak in output:
                # If the event occurs within the trig_int time then do not add
                # it, and break out of the inner loop.
                if abs(coincidence_trigger[1] - peak[1]) < (trig_int *
                                                            samp_rate):
                    add = False
                    break
            if add:
                output.append((coincidence_trigger[0],
                               coincidence_trigger[1]))
        output.sort(key=lambda tup: tup[1])
        return output
    else:
        return []


if __name__ == "__main__":
    import doctest
    doctest.testmod()
