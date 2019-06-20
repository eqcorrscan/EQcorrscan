"""Functions to find peaks in data above a certain threshold.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import ctypes
import random
import logging
import numpy as np

from multiprocessing import Pool, cpu_count
from scipy import ndimage
from future.utils import native_str

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
        If True, will declustering within data-sections above the threshold,
        rather than just taking the peak within that section. This will take
        more time. This defaults to True for match_filter.

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
        peaks = decluster(
            peaks=np.array(peak_vals), index=np.array(peak_indices),
            trig_int=trig_int + 1, threshold=thresh)
        peaks = sorted(peaks, key=lambda peak: peak[1], reverse=False)
        return peaks
    else:
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
        Logger.info('No peaks for you!')
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
            cores = min(arr.shape[0], cpu_count())
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


def _multi_find_peaks_compiled(arrays, thresholds, trig_int, full_peaks,
                               cores):
    """
    Determine peaks in an array or arrays of data above a certain threshold.

    :type arrays: numpy.ndarray
    :param arrays: 2-D numpy array is required
    :type thresholds: list
    :param thresholds:
        Minimum value for peaks.
    :type trig_int: int
    :param trig_int:
        The minimum difference in samples between triggers, if multiple
        peaks within this window this code will find the highest.
    :type full_peaks: bool
    :param full_peaks:
        If True, will decluster within data-sections above the threshold,
        rather than just taking the peak within that section. This will take
        more time. This defaults to True for match_filter.
    :type cores: int
    :param cores: Number of threads to parallel across

    :return: peaks: List of List of tuples of peak values and locations.
    :rtype: list
    """
    if not full_peaks:
        peak_vals, peak_indices = _multi_find_peaks_c(
            arrays=arrays, thresholds=thresholds, threads=cores)
        # Remove empty arrays
        peak_mapper = {}
        map_index = 0
        _peak_vals = []
        _peak_indices = []
        _thresholds = []
        for i in range(arrays.shape[0]):
            if len(peak_vals[i]) > 0:
                peak_mapper.update({i: map_index})
                _peak_vals.append(peak_vals[i])
                _peak_indices.append(peak_indices[i])
                _thresholds.append(thresholds[i])
                map_index += 1
        peak_vals = _peak_vals
        peak_indices = _peak_indices
        thresholds = _thresholds
    else:
        peak_vals = arrays
        peak_indices = [np.arange(arr.shape[0]) for arr in arrays]
        peak_mapper = {i: i for i in range(len(peak_indices))}
    if len(peak_indices) > 0:
        peaks = _multi_decluster(
            peaks=peak_vals, indices=peak_indices, trig_int=trig_int,
            thresholds=thresholds, cores=cores)
        peaks = [sorted(_peaks, key=lambda peak: peak[1], reverse=False)
                 for _peaks in peaks]
    out_peaks = []
    for i in range(arrays.shape[0]):
        if i in peak_mapper.keys():
            out_peaks.append(peaks[peak_mapper[i]])
        else:
            out_peaks.append([])
    return out_peaks


def _multi_decluster(peaks, indices, trig_int, thresholds, cores):
    """
    Decluster peaks based on an enforced minimum separation.

    Only works when peaks and indices are all the same shape.

    :type peaks: list
    :param peaks: list of arrays of peak values
    :type indices: list
    :param indices: list of arrays of locations of peaks
    :type trig_int: int
    :param trig_int: Minimum trigger interval in samples
    :type thresholds: list
    :param thresholds: list of float of threshold values

    :return: list of lists of tuples of (value, sample)
    """
    utilslib = _load_cdll('libutils')

    lengths = np.array([peak.shape[0] for peak in peaks], dtype=int)
    trig_int = int(trig_int)
    n = np.int32(len(peaks))
    cores = min(cores, n)

    total_length = lengths.sum()

    max_indexes = [_indices.max() for _indices in indices]
    max_index = max(max_indexes)
    for var in [trig_int, lengths.max(), max_index]:
        if var == ctypes.c_long(var).value:
            long_type = ctypes.c_long
            func = utilslib.multi_decluster
        elif var == ctypes.c_longlong(var).value:
            long_type = ctypes.c_longlong
            func = utilslib.multi_decluster_ll
        else:
            # Note, could use numpy.gcd to try and find greatest common
            # divisor and make numbers smaller
            raise OverflowError("Maximum index larger than internal long long")

    func.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(total_length,),
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=long_type, shape=(total_length,),
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=long_type, shape=(n,),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(n,),
                               flags=native_str('C_CONTIGUOUS')),
        long_type,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(total_length,),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int]
    func.restype = ctypes.c_int

    peaks_sorted = np.empty(total_length, dtype=np.float32)
    indices_sorted = np.empty_like(peaks_sorted, dtype=np.float32)

    # TODO: When doing full decluster from match-filter, all lengths will be
    # TODO: the same - would be more efficient to use numpy sort on 2D matrix
    start_ind = 0
    end_ind = 0
    for _peaks, _indices, length in zip(peaks, indices, lengths):
        end_ind += length
        sorted_indices = np.abs(_peaks).argsort()
        peaks_sorted[start_ind: end_ind] = _peaks[sorted_indices[::-1]]
        indices_sorted[start_ind: end_ind] = _indices[sorted_indices[::-1]]
        start_ind += length

    peaks_sorted = np.ascontiguousarray(peaks_sorted, dtype=np.float32)
    indices_sorted = np.ascontiguousarray(
        indices_sorted, dtype=long_type)
    lengths = np.ascontiguousarray(lengths, dtype=long_type)
    thresholds = np.ascontiguousarray(thresholds, dtype=np.float32)
    out = np.zeros(total_length, dtype=np.uint32)
    ret = func(
        peaks_sorted, indices_sorted, lengths, np.int32(n), thresholds,
        long_type(trig_int + 1), out, np.int32(cores))
    if ret != 0:
        raise MemoryError("Issue with c-routine, returned %i" % ret)

    peaks_out = []
    slice_start = 0
    for length in lengths:
        slice_end = slice_start + length
        out_mask = out[slice_start: slice_end].astype(bool)
        declustered_peaks = peaks_sorted[slice_start: slice_end][out_mask]
        declustered_indices = indices_sorted[slice_start: slice_end][out_mask]
        peaks_out.append(list(zip(declustered_peaks, declustered_indices)))
        slice_start = slice_end
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

    length = peaks.shape[0]
    trig_int = int(trig_int)

    for var in [index.max(), trig_int]:
        if var == ctypes.c_long(var).value:
            long_type = ctypes.c_long
            func = utilslib.decluster
        elif var == ctypes.c_longlong(var).value:
            long_type = ctypes.c_longlong
            func = utilslib.decluster_ll
        else:
            raise OverflowError("Maximum index larger than internal long long")

    func.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=long_type, shape=(length,),
                               flags=native_str('C_CONTIGUOUS')),
        long_type, ctypes.c_float, long_type,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length,),
                               flags=native_str('C_CONTIGUOUS'))]
    func.restype = ctypes.c_int

    sorted_inds = np.abs(peaks).argsort()
    arr = peaks[sorted_inds[::-1]]
    inds = index[sorted_inds[::-1]]
    arr = np.ascontiguousarray(arr, dtype=np.float32)
    inds = np.ascontiguousarray(inds, dtype=long_type)
    out = np.zeros(len(arr), dtype=np.uint32)

    ret = func(
        arr, inds, long_type(length), np.float32(threshold),
        long_type(trig_int), out)
    if ret != 0:
        raise MemoryError("Issue with c-routine, returned %i" % ret)

    peaks_out = list(zip(arr[out.astype(bool)], inds[out.astype(bool)]))
    return peaks_out


def _find_peaks_c(array, threshold):
    """
    Use a C func to find peaks in the array.
    """
    utilslib = _load_cdll('libutils')

    length = array.shape[0]
    utilslib.find_peaks.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length, ),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long, ctypes.c_float,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length, ),
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.find_peaks.restype = ctypes.c_int
    arr = np.ascontiguousarray(array, np.float32)
    out = np.ascontiguousarray(np.zeros((length, ), dtype=np.uint32))
    ret = utilslib.find_peaks(arr, ctypes.c_long(length), threshold, out)

    if ret != 0:
        raise MemoryError("Internal error")

    peaks_locations = np.nonzero(out)
    return array[peaks_locations], peaks_locations[0]


def _multi_find_peaks_c(arrays, thresholds, threads):
    """
    Wrapper for multi-find peaks C-func
    """
    utilslib = _load_cdll('libutils')

    length = arrays.shape[1]
    n = np.int32(arrays.shape[0])
    thresholds = np.ascontiguousarray(thresholds, np.float32)
    arr = np.ascontiguousarray(arrays.flatten(), np.float32)
    utilslib.multi_find_peaks.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(n * length,),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_long, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(n, ),
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(n * length, ),
                               flags=native_str('C_CONTIGUOUS'))]
    utilslib.multi_find_peaks.restype = ctypes.c_int

    out = np.ascontiguousarray(np.zeros((n * length, ), dtype=np.uint32))
    ret = utilslib.multi_find_peaks(
        arr, ctypes.c_long(length), n, thresholds, threads, out)
    # Copy data to avoid farking the users data
    if ret != 0:
        raise MemoryError("Internal error")
    peaks = []
    peak_locations = []
    out = out.reshape(n, length)
    for i in range(n):
        peak_locs = np.nonzero(out[i])
        peaks.append(arrays[i][peak_locs])
        peak_locations.append(peak_locs[0])

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
