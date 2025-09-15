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
import warnings

import numpy as np

from typing import List, Union

from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor
from scipy import ndimage
from obspy import UTCDateTime
from obspy.core.event import Catalog, Event

from eqcorrscan.utils.correlate import pool_boy
from eqcorrscan.utils.libnames import _load_cdll
from eqcorrscan.utils.clustering import dist_mat_km


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
        If True, will decluster within data-sections above the threshold,
        rather than just taking the peak within that section. This will take
        more time. This defaults to False for match_filter.

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
    peaks = []
    if len(peak_vals) > 0:
        peaks = decluster(
            peaks=np.array(peak_vals), index=np.array(peak_indices),
            trig_int=trig_int + 1, threshold=thresh)
        peaks = sorted(peaks, key=lambda peak: peak[1], reverse=False)
    return peaks


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
        If True will by decluster within
        data-sections above the threshold, rather than just taking the peak
        within that section. This will take more time. This defaults to False
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
        Logger.debug("All values above threshold, running full peak finding")
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
        # Convert peaks to python types
        peaks = [(p[0].item(), p[1].item()) for p in peaks]
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
    if not parallel or cores == 1:
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
            to_run = ((arr[i], thresh[i], trig_int)
                      for i in range(len(thresh)))
            with ThreadPoolExecutor(cores) as executor:
                results = executor.map(
                    lambda args: find_peaks_compiled(*args), to_run)
            peaks = [r for r in results]
    return peaks


def decluster_distance_time(peaks, index, trig_int, catalog,
                            hypocentral_separation, threshold=0,
                            num_threads=None):
    """
    Decluster based on time between peaks, and distance between events.

    Peaks, index and catalog must all be sorted the same way, e.g. peak[i]
    corresponds to index[i] and catalog[i]. Peaks that are within the time
    threshold of one-another, but correspond to events separated by more than
    the hypocentral_separation threshold will not be removed.

    :type peaks: np.array
    :param peaks: array of peak values
    :type index: np.ndarray
    :param index: locations of peaks
    :type trig_int: int
    :param trig_int: Minimum trigger interval in samples
    :type catalog: obspy.core.event.Catalog
    :param catalog:
        Catalog of events with origins to use to measure inter-event distances
    :type hypocentral_separation: float
    :param hypocentral_separation:
        Maximum inter-event distance to decluster over in km
    :type threshold: float
    :param threshold: Minimum absolute peak value to retain it
    :type num_threads: int
    :param num_threads:
        Number of threads to use for distance matrix calculation.

    :return: list of tuples of (value, sample)
    """
    utilslib = _load_cdll("libutils")
    length = peaks.shape[0]
    trig_int = int(trig_int)

    for var in [index.max(), trig_int]:
        if var == ctypes.c_long(var).value:
            long_type = ctypes.c_long
            func = utilslib.decluster_dist_time
        elif var == ctypes.c_longlong(var).value:
            long_type = ctypes.c_longlong
            func = utilslib.decluster_dist_time_ll
        else:
            raise OverflowError("Maximum index larger than internal long long")

    func.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length,),
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=long_type, shape=(length,),
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length * length,),
                               flags='C_CONTIGUOUS'),
        long_type, ctypes.c_float, long_type, ctypes.c_float,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length,),
                               flags='C_CONTIGUOUS')]
    func.restype = ctypes.c_int

    sorted_inds = np.abs(peaks).argsort()
    # Sort everything in the same way.
    arr = peaks[sorted_inds[::-1]]
    inds = index[sorted_inds[::-1]]
    sorted_events = [catalog[i] for i in sorted_inds[::-1]]
    distance_matrix = dist_mat_km(
        catalog=sorted_events, num_threads=num_threads)

    arr = np.ascontiguousarray(arr, dtype=np.float32)
    inds = np.ascontiguousarray(inds, dtype=long_type)
    distance_matrix = np.ascontiguousarray(
        distance_matrix.flatten(order="C"), dtype=np.float32)
    out = np.zeros(len(arr), dtype=np.uint32)

    ret = func(
        arr, inds, distance_matrix, long_type(length), np.float32(threshold),
        long_type(trig_int), hypocentral_separation, out)
    if ret != 0:
        raise MemoryError("Issue with c-routine, returned %i" % ret)

    peaks_out = list(zip(arr[out.astype(bool)], inds[out.astype(bool)]))
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
    :type threshold: float
    :param threshold: Minimum absolute peak value to retain it.

    :return: list of tuples of (value, sample)
    """
    utilslib = _load_cdll("libutils")
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
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=long_type, shape=(length,),
                               flags='C_CONTIGUOUS'),
        long_type, ctypes.c_float, long_type,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length,),
                               flags='C_CONTIGUOUS')]
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
    utilslib = _load_cdll("libutils")
    length = array.shape[0]
    utilslib.find_peaks.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, shape=(length, ),
                               flags='C_CONTIGUOUS'),
        ctypes.c_long, ctypes.c_float,
        np.ctypeslib.ndpointer(dtype=np.uint32, shape=(length, ),
                               flags='C_CONTIGUOUS')]
    utilslib.find_peaks.restype = ctypes.c_int
    arr = np.ascontiguousarray(array, np.float32)
    out = np.ascontiguousarray(np.zeros((length, ), dtype=np.uint32))
    ret = utilslib.find_peaks(arr, ctypes.c_long(length), threshold, out)

    if ret != 0:
        raise MemoryError("Internal error")

    peaks_locations = np.nonzero(out)
    return array[peaks_locations], peaks_locations[0]


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


# ##################### Declustering based on pick-time #######################


def _pick_time(
    event: Event,
    station: str,
    phase_hint: str,
    zero_time: UTCDateTime
) -> float:
    picks = [p for p in event.picks
             if p.waveform_id.station_code == station
             and p.phase_hint == phase_hint]
    if len(picks) == 0:
        return np.nan
    pick_times = [p.time - zero_time for p in picks]
    if len(pick_times) == 1:
        return pick_times[0]
    return sum(pick_times) / len(pick_times)


def average_pick_time_diff_matrix(
    catalog: Union[Catalog, List[Event]],
    compiled: bool = True
) -> np.ndarray:
    # Make arrays of pick times for each phase and each station
    sta_phases = list(
        {f"{p.waveform_id.station_code}_{p.phase_hint}"
         for ev in catalog for p in ev.picks if p.phase_hint in "PS"})
    # Using a list for sta_phases to retain order
    pick_arrays = dict()
    zero_time = (catalog[0].preferred_origin() or catalog[0].origins[-1]).time
    # This bit is not very fast, but gets the arrays in C shape
    for sta_ph in sta_phases:
        pick_arrays.update({
            sta_ph: np.array(
                [_pick_time(ev, *sta_ph.split('_'), zero_time=zero_time)
                 for ev in catalog])})
    if compiled:
        return _compute_matrix_c(
            pick_arrays=pick_arrays, n_events=len(catalog))
    else:
        return _compute_matrix(
            pick_arrays=pick_arrays, n_events=len(catalog))


def _compute_matrix(pick_arrays: dict, n_events: int) -> np.ndarray:
    # This bit could be ported to C
    out = np.zeros((n_events, n_events))
    for i in range(n_events):
        # Matrix is symmetric, so we can just do upper or lower half.
        for j in range(i + 1, n_events):
            # We get occasional RuntimeWarnings about mean of empty slice
            # when no stations are shared. nan is returned for this case.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pick_delta = np.nanmean(
                    [pick_arrays[sta_ph][i] - pick_arrays[sta_ph][j]
                     for sta_ph in pick_arrays.keys()])
            out[i, j] = pick_delta
    out -= out.T
    return np.nan_to_num(out, nan=np.inf)


def _compute_matrix_c(pick_arrays: dict, n_events: int) -> np.ndarray:
    out = np.ascontiguousarray(
        np.zeros(n_events * n_events), np.float64)
    pick_array = np.ascontiguousarray(
        np.concatenate([value for value in pick_arrays.values()]), np.float64)

    utilslib = _load_cdll("libutils")
    utilslib.average_pick_time_diff_matrix.argtypes = [
        ctypes.c_int, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, shape=(len(pick_array), ),
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, shape=(len(out), ),
                               flags='C_CONTIGUOUS')]
    utilslib.average_pick_time_diff_matrix.restype = ctypes.c_int

    ret = utilslib.average_pick_time_diff_matrix(
        ctypes.c_int(len(pick_arrays.keys())), ctypes.c_int(n_events),
        pick_array, out)

    out = out.reshape(n_events, n_events)
    out -= out.T
    return np.nan_to_num(out, nan=np.inf)


def get_det_val(event: Event) -> Union[float, None]:
    """
    Get the detection value for an event. Prefers sum of pick correlations.


    """
    det_val = get_comment_val(value_name="detect_val", event=event)
    pick_sum = 0
    for pick in event.picks:
        pick_val = get_comment_val("cc_max", event=pick)
        if pick_val:
            pick_sum += pick_val
    if det_val:
        if pick_sum > det_val:
            return pick_sum
        else:
            return det_val
    elif pick_sum:
        return pick_sum
    return None


def get_comment_val(value_name: str, event: Event) -> Union[float, None]:
    value = None
    for comment in event.comments:
        if value_name in comment.text:
            if "=" in comment.text:
                # Should be a number
                try:
                    value = float(comment.text.split('=')[-1])
                except ValueError:
                    # Leave as a string
                    break
                except Exception as e:
                    Logger.exception(
                        f"Could not get {value_name} {comment.text} "
                        f"due to {e}")
                else:
                    break
            elif ":" in comment.text:
                # Should be a string
                try:
                    value = comment.text.split(": ")[-1]
                except Exception as e:
                    Logger.exception(
                        f"Could not get {value_name} from {comment.text} "
                        f"due to {e}")
                else:
                    break
    return value


def decluster_pick_times(
    catalog: Union[Catalog, List[Event]],
    trig_int: float
) -> Catalog:
    detect_vals = np.array([abs(get_det_val(ev) or len(ev.picks))
                            for ev in catalog])
    # Sort catalog in order of detect_vals - from largest to smallest
    order = np.argsort(detect_vals)[::-1]
    # detect_vals = detect_vals[order]
    if isinstance(catalog, Catalog):
        events = catalog.events
    else:
        events = catalog

    events = [events[ind] for ind in order]

    # This is obviously very slow at the moment
    pick_matrix = average_pick_time_diff_matrix(catalog=events)
    # We only care about the absolute value (+/- trig_int)
    pick_matrix = np.abs(pick_matrix)

    drop_events = np.zeros(len(events), dtype=bool)
    # Loop through from highest detection value to lowest - this way we keep
    # the best detections
    for i, ev in enumerate(events):
        # If we are dropping this event then carry on
        if drop_events[i]:
            continue
        # Find any events within trig int and set their drop status to True
        drop_mask = pick_matrix[i] <= trig_int
        # Make sure we retain the "core" event in this loop
        drop_mask[i] = False
        drop_events[drop_mask] = True
        # _dropped = [ev for j, ev in _dropped(events) if drop_mask[j]]

    # Get the events that we want to keep!
    return Catalog([ev for i, ev in enumerate(events) if not drop_events[i]])


if __name__ == "__main__":
    import doctest
    doctest.testmod()
