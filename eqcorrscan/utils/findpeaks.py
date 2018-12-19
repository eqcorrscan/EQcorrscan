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
import logging
import numpy as np

from multiprocessing import Pool
from scipy import ndimage
from future.utils import native_str
from itertools import compress

from eqcorrscan.utils.correlate import pool_boy
from eqcorrscan.utils.libnames import _load_cdll


Logger = logging.getLogger(__name__)


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


def find_peaks_compiled(arr, thresh, trig_int, full_peaks=False):
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
    :type full_peaks: bool
    :param full_peaks:
        If True, will remove the issue eluded to below, by declustering within
        data-sections above the threshold, rather than just taking the peak
        within that section. This will take more time. This defaults to True
        for match_filter.

    :return: peaks: Lists of tuples of peak values and locations.
    :rtype: list
    """
    if not np.any(np.abs(arr) > thresh):
        # Fast fail
        return []
    if not full_peaks:
        peak_vals, peak_indices = _find_peaks_c(array=arr, threshold=thresh)
    else:
        peak_vals = arr
        peak_indices = np.arange(arr.shape[0])
    if len(peak_vals) > 0:
        peaks = decluster(peaks=np.array(peak_vals),
                          index=np.array(peak_indices), trig_int=trig_int + 1)
        peaks = sorted(peaks, key=lambda time: time[1], reverse=False)
        return peaks
    else:
        print('No peaks for you!')
        return []


def find_peaks2_short(arr, thresh, trig_int, full_peaks=False):
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
        trigger intervals are usually relatively large and thresholds high.

    """
    # Set everything below the threshold to zero
    image = np.copy(arr)
    Logger.debug("Threshold: {0}\tMax: {1}".format(thresh, max(image)))
    image[np.abs(image) < thresh] = 0
    if len(image[np.abs(image) > thresh]) == 0:
        Logger.debug("No values over threshold {0}".format(thresh))
        return []
    if np.all(np.abs(arr) > thresh):
        full_peaks = True
    Logger.debug('Found {0} samples above the threshold'.format(
        len(image[image > thresh])))
    initial_peaks = []
    # Find the peaks
    labeled_image, number_of_objects = ndimage.label(image)
    peak_slices = ndimage.find_objects(labeled_image)
    for peak_slice in peak_slices:
        window = arr[peak_slice[0].start: peak_slice[0].stop]
        if peak_slice[0].stop - peak_slice[0].start >= trig_int and full_peaks:
            window_peaks, window_peak_indexes = ([], [])
            for i in np.arange(peak_slice[0].start, peak_slice[0].stop):
                if i == peak_slice[0].start:
                    prev_value = 0
                else:
                    prev_value = arr[i - 1]
                if i == peak_slice[0].stop - 1:
                    next_value = 0
                else:
                    next_value = arr[i + 1]
                # Check for consistent sign - either both greater or
                # both smaller.
                if (next_value - arr[i]) * (prev_value - arr[i]) > 0:
                    window_peaks.append(arr[i])
                    window_peak_indexes.append(i)
            peaks = decluster(
                peaks=np.array(window_peaks), trig_int=trig_int + 1,
                index=np.array(window_peak_indexes))
        else:
            peaks = [(window[np.argmax(abs(window))],
                      int(peak_slice[0].start + np.argmax(abs(window))))]
        initial_peaks.extend(peaks)
    peaks = decluster(peaks=np.array(list(zip(*initial_peaks))[0]),
                      index=np.array(list(zip(*initial_peaks))[1]),
                      trig_int=trig_int + 1)
    if initial_peaks:
        peaks = sorted(peaks, key=lambda time: time[1], reverse=False)
        return peaks
    else:
        print('No peaks for you!')
        return []


def multi_find_peaks(arr, thresh, trig_int, parallel=True, full_peaks=False,
                     cores=None, internal_func=find_peaks_compiled):
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
    :type parallel: bool
    :param parallel:
        Whether to compute in parallel or not - will use multiprocessing if not
        using the compiled internal_func
    :type full_peaks: bool
    :param full_peaks: See `eqcorrscan.utils.findpeaks.find_peaks2_short`
    :type cores: int
    :param cores:
        Maximum number of processes to spin up for parallel peak-finding
    :type internal_func: callable
    :param internal_func:
        Function to use for peak finding - defaults to the compiled version.

    :returns:
        List of list of tuples of (peak, index) in same order as input arrays
    """
    peaks = []
    if not parallel:
        for sub_arr, arr_thresh in zip(arr, thresh):
            peaks.append(internal_func(
                arr=sub_arr, thresh=arr_thresh, trig_int=trig_int,
                full_peaks=full_peaks))
    else:
        if cores is None:
            cores = arr.shape[0]
        if internal_func.__name__ != 'find_peaks_compiled':
            with pool_boy(Pool=Pool, traces=arr.shape[0], cores=cores) as pool:
                params = ((sub_arr, arr_thresh, trig_int, full_peaks)
                          for sub_arr, arr_thresh in zip(arr, thresh))
                results = [
                    pool.apply_async(internal_func, param) for param in params]
                peaks = [res.get() for res in results]
        else:
            peaks = _multi_find_peaks_compiled(
                arr, thresh, trig_int, full_peaks=full_peaks, cores=cores)
    return peaks


def _multi_find_peaks_compiled(arrays, thresholds, trig_int, full_peaks, cores):
    """

    :param arrays:
    :param thresholds:
    :param trig_int:
    :param full_peaks:
    :param cores:
    :return:
    """
    peaks = []
    #TODO: multi-decluster compiled.
    return peaks


def _multi_decluster(peaks, indices, trig_int, thresholds):
    """
    Decluster peaks based on an enforced minimum separation.

    Only works when peaks and indices are all the same shape.

    :type peaks: np.ndarray
    :param peaks: array of peak values
    :type indices: np.ndarray
    :param indices: locations of peaks
    :type trig_int: list of int
    :param trig_int: Minimum trigger interval in samples

    :return: list of tuples of (value, sample)
    """
    utilslib = _load_cdll('libutils')
    # TODO Write with mixed sizes

    length = np.int64(peaks.shape[1])
    n = np.int32(peaks.shape[0])
    utilslib.multi_decluster.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=np.int64, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_float, ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.multi_decluster.restype = ctypes.c_int
    peaks_sort = sorted(zip(peaks, index),
                        key=lambda amplitude: abs(amplitude[0]),
                        reverse=True)
    arr, inds = zip(*peaks_sort)
    arr = np.ascontiguousarray(arr, dtype=np.float32)
    inds = np.ascontiguousarray(inds, dtype=np.int64)
    out = np.zeros(len(arr), dtype=np.uint32)
    ret = utilslib.decluster(
        arr, inds, length, np.float32(threshold), np.int64(trig_int), out)
    if ret != 0:
        raise MemoryError("Issue with c-routine, returned %i" % ret)
    peaks_out = list(compress(peaks_sort, out))
    return peaks_out


def decluster(peaks, index, trig_int, threshold=0):
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

    length = np.int64(peaks.shape[0])
    utilslib.decluster.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=np.int64, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_float, ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.decluster.restype = ctypes.c_int
    peaks_sort = sorted(zip(peaks, index),
                        key=lambda amplitude: abs(amplitude[0]),
                        reverse=True)
    arr, inds = zip(*peaks_sort)
    arr = np.ascontiguousarray(arr, dtype=np.float32)
    inds = np.ascontiguousarray(inds, dtype=np.int64)
    out = np.zeros(len(arr), dtype=np.uint32)
    ret = utilslib.decluster(
        arr, inds, length, np.float32(threshold), np.int64(trig_int), out)
    if ret != 0:
        raise MemoryError("Issue with c-routine, returned %i" % ret)
    peaks_out = list(compress(peaks_sort, out))
    return peaks_out


def _find_peaks_c(array, threshold):
    """
    Use a C func to find peaks in the array.
    """
    utilslib = _load_cdll('libutils')

    length = np.int64(array.shape[0])
    utilslib.find_peaks.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length, ),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long, ctypes.c_float]
    utilslib.find_peaks.restype = ctypes.c_int
    array = np.ascontiguousarray(array, np.float32)
    ret = utilslib.find_peaks(array, length, threshold)
    if ret != 0:
        raise MemoryError("Internal error")
    peaks_locations = np.nonzero(array)
    return array[peaks_locations], peaks_locations[0]


def _multi_find_peaks_c(arrays, thresholds, threads):
    """
    Wrapper for multi-find peaks C-func
    """
    utilslib = _load_cdll('libutils')

    length = np.int64(arrays.shape[1])
    n = np.int32(arrays.shape[0])
    thresholds = np.ascontiguousarray(thresholds, np.float32)
    arrays = np.ascontiguousarray(arrays, np.float32)
    utilslib.multi_find_peaks.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(n, length),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(n, ),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int]
    utilslib.multi_find_peaks.restype = ctypes.c_int
    ret = utilslib.multi_find_peaks(arrays, length, n, thresholds, threads)
    if ret != 0:
        raise MemoryError("Internal error")
    peaks = []
    peak_locations = []
    for i in range(n):
        peak_locs = np.nonzero(arrays[i])
        peaks.append(arrays[i][peak_locs])
        peak_locations.append(peak_locs)

    return peaks, peak_locations


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
