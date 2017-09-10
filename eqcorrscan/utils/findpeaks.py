r"""Functions to find peaks in data above a certain threshold.

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

import random
import numpy as np

from obspy import UTCDateTime
from scipy import ndimage


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
                      samp_rate=1.0):
    """
    Determine peaks in an array of data above a certain threshold.

    Uses a mask to remove data below threshold and finds peaks in what is left.

    :type arr: numpy.ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh: The threshold below which will be considered noise and \
        peaks will not be found in.
    :type trig_int: int
    :param trig_int: The minimum difference in samples between triggers,\
        if multiple peaks within this window this code will find the highest.
    :type debug: int
    :param debug: Optional, debug level 0-5
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Starttime for plotting, only used if debug > 2.
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz, only used for plotting if debug > 2.

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
    if not starttime:
        starttime = UTCDateTime(0)
    # Set everything below the threshold to zero
    image = np.copy(arr)
    image = np.abs(image)
    image[image < thresh] = 0
    if len(image[image > thresh]) == 0:
        if debug > 0:
            print('No values over threshold found')
        return []
    if debug > 0:
        print(' '.join(['Found', str(len(image[image > thresh])),
                        'samples above the threshold']))
    initial_peaks = []
    peaks = []
    # Find the peaks
    labeled_image, number_of_objects = ndimage.label(image)
    peak_slices = ndimage.find_objects(labeled_image)
    for peak_slice in peak_slices:
        # print('Width of peak='+str(peak_slice[0].stop-peak_slice[0].start)
        window = arr[peak_slice[0].start: peak_slice[0].stop]
        initial_peaks.append((max(window),
                              int(peak_slice[0].start + np.argmax(window))))
    if initial_peaks:
        peaks = decluster(peaks=initial_peaks, trig_int=trig_int, debug=debug)
        if debug >= 3:
            from eqcorrscan.utils import plotting
            _fname = ''.join(['peaks_',
                              starttime.datetime.strftime('%Y-%m-%d'),
                              '.pdf'])
            plotting.peaks_plot(data=image, starttime=starttime,
                                samp_rate=samp_rate, save=True,
                                peaks=peaks, savefile=_fname)
        peaks = sorted(peaks, key=lambda time: time[1], reverse=False)
        return peaks
    else:
        print('No peaks for you!')
        return peaks


def decluster(peaks, trig_int, return_ind=False, debug=0):
    """
    Decluster peaks based on an enforced minimum separation.

    :type peaks: list
    :param peaks: list of tuples of (value, sample)
    :type trig_int: int
    :param trig_int: Minimum trigger interval in samples
    :type return_ind: bool
    :param return_ind:
        Whether to also return the indices of the original peaks or not.

    :return: list of tuples of (value, sample)
    """
    peaks_sort = sorted(zip(peaks, np.arange(len(peaks))),
                        key=lambda amplitude: amplitude[0][0],
                        reverse=True)
    peaks_sort, inds = zip(*peaks_sort)
    # Debugging
    if debug >= 4:
        for peak in peaks:
            print(peak)
    peaks_out = []
    ind_out = []
    add = False
    ind_out.append(inds[0])
    peaks_out.append(peaks_sort[0])  # Definitely take the biggest peak
    if debug > 3:
        print(' '.join(['Added the biggest peak of', str(peaks_out[0][0]),
                        'at sample', str(peaks_out[0][1])]))
    if len(peaks) > 1:
        if debug > 3:
            msg = ' '.join(['Multiple peaks found, checking',
                            'them now to see if they overlap'])
            print(msg)
        for next_peak, next_ind in zip(peaks_sort, inds):
            # i in range(1,len(peaks_sort)):
            # Loop through the amplitude sorted peaks
            # if the next highest amplitude peak is within trig_int of any
            # peak already in peaks then we don't want it, else, add it
            # next_peak=peaks_sort[i]
            if debug > 3:
                print(next_peak)
            for peak in peaks_out:
                # Use add as a switch for whether or not to append
                # next peak to peaks, if once gone through all the peaks
                # it is True, then we will add it, otherwise we won't!
                if abs(next_peak[1] - peak[1]) < trig_int:
                    if debug > 3:
                        msg = ' '.join(['Difference in time is',
                                        str(next_peak[1] - peak[1]), '\n',
                                        'Which is less than',
                                        str(trig_int)])
                        print(msg)
                    add = False
                    # Need to exit the loop here if false
                    break
                else:
                    add = True
            if add:
                if debug > 3:
                    msg = ' '.join(['Adding peak of', str(next_peak[0]),
                                    'at sample', str(next_peak[1])])
                    print(msg)
                ind_out.append(next_ind)
                peaks_out.append(next_peak)
            elif debug > 3:
                msg = ' '.join(['I did not add peak of', str(next_peak[0]),
                                'at sample', str(next_peak[1])])
                print(msg)
    if return_ind:
        return peaks_out, ind_out
    else:
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
