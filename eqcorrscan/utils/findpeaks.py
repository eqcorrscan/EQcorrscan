r"""Functions to find peaks in data above a certain threshold.

:copyright:
    Calum Chamberlain, Chet Hopp.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


def is_prime(number):
    r"""Function to test primality of a number. Function lifted from online \
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
    import random
    ''' if number != 1 '''
    if (number > 1):
        ''' repeat the test few times '''
        for time in range(3):
            ''' Draw a RANDOM number in range of number ( Z_number )  '''
            randomNumber = random.randint(2, number)-1
            ''' Test if a^(n-1) = 1 mod n '''
            if (pow(randomNumber, number-1, number) != 1):
                return False
        return True
    else:
        ''' case number == 1 '''
        return False

# Note this function doesn't give the expected results.
# def find_peaks2(arr, thresh, trig_int, debug=0, maxwidth=10,
#                 starttime=False, samp_rate=1.0):
#     r"""Function to determine peaks in an array of data using scipy \
#     find_peaks_cwt, works fast in certain cases, but for match_filter cccsum \
#     peak finding, find_peaks2_short works better.  Test it out and see which \
#     works best for your application.
#
#     :type arr: ndarray
#     :param arr: 1-D numpy array is required
#     :type thresh: float
#     :param thresh: The threshold below which will be considered noise and \
#     peaks will not be found in.
#     :type trig_int: int
#     :param trig_int: The minimum difference in samples between triggers, \
#     if multiple peaks within this window this code will find the highest.
#     :type debug: int
#     :param debug: Optional, debug level 0-5
#     :type maxwidth: int
#     :param maxwidth: Maximum peak width to look for in samples
#     :type starttime: osbpy.UTCDateTime
#     :param starttime: Starttime for plotting, only used if debug > 2.
#     :type samp_rate: float
#     :param samp_rate: Sampling rate in Hz, only used for plotting if debug > 2.
#
#     :return: peaks: Lists of tuples of peak values and locations.
#
#     >>> import numpy as np
#     >>> arr = np.random.randn(100)
#     >>> threshold = arr.max()
#     >>> arr[40] = 20
#     >>> arr[60] = 100
#     >>> find_peaks2(arr, threshold, 3)
#     Finding peaks
#     sorting peaks
#     [(20.0, 40), (100.0, 60)]
#     """
#     from scipy.signal import find_peaks_cwt
#     import numpy as np
#     from obspy import UTCDateTime
#     if not starttime:
#         starttime = UTCDateTime(0)
#     # Set everything below the threshold to zero
#     image = np.copy(arr)
#     image = np.abs(image)
#     image[image < thresh] = thresh
#     # We need to check if the number of samples in the image is prime, if it
#     # is this method will be really slow, so we add a pad to the end to make
#     # it not of prime length!
#     if is_prime(len(image)):
#         image = np.append(image, 0.0)
#         print('Input array has a prime number of samples, appending a zero')
#         print(len(image))
#     if len(image[image > thresh]) == 0:
#         print('No values over threshold found')
#         return []
#     if debug > 0:
#         msg = ' '.join(['Found', str(len(image[image > thresh])),
#                         'samples above the threshold'])
#         print(msg)
#     initial_peaks = []
#     peaks = []
#     # Find the peaks
#     print('Finding peaks')
#     peakinds = find_peaks_cwt(image, np.arange(1, maxwidth))
#     initial_peaks = [(image[peakind], peakind) for peakind in peakinds]
#     # Sort initial peaks according to amplitude
#     print('sorting peaks')
#     peaks_sort = sorted(initial_peaks, key=lambda amplitude: amplitude[0],
#                         reverse=True)
#     if debug >= 4:
#         for peak in initial_peaks:
#             print(peak)
#     if initial_peaks:
#         peaks.append(peaks_sort[0])  # Definitely take the biggest peak
#         if debug > 3:
#             msg = ' '.join(['Added the biggest peak of', str(peaks[0][0]),
#                             'at sample', str(peaks[0][1])])
#             print(msg)
#         if len(initial_peaks) > 1:
#             if debug > 3:
#                 msg = ' '.join(['Multiple peaks found, checking them',
#                                 'now to see if they overlap'])
#                 print(msg)
#             for next_peak in peaks_sort:
#                 # i in range(1,len(peaks_sort)):
#                 # Loop through the amplitude sorted peaks
#                 # if the next highest amplitude peak is within trig_int of any
#                 # peak already in peaks then we don't want it, else, add it
#                 # next_peak = peaks_sort[i]
#                 if debug > 3:
#                     print(next_peak)
#                 for peak in peaks:
#                     add = False
#                     # Use add as a switch for whether or not to append
#                     # next peak to peaks, if once gone through all the peaks
#                     # it is True, then we will add it, otherwise we won't!
#                     if abs(next_peak[1] - peak[1]) < trig_int:
#                         if debug > 3:
#                             msg = ' '.join(['Difference in time is',
#                                             str(next_peak[1] - peak[1]), '\n'
#                                             'Which is less than',
#                                             str(trig_int)])
#                             print(msg)
#                         add = False
#                         # Need to exit the loop here if false
#                         break
#                     else:
#                         add = True
#                 if add:
#                     if debug > 3:
#                         msg = ' '.join(['Adding peak of', str(next_peak[0]),
#                                         'at sample', str(next_peak[1])])
#                         print(msg)
#                     peaks.append(next_peak)
#                 elif debug > 3:
#                     msg = ' '.join(['I did not add peak of',
#                                     str(next_peak[0]), 'at sample',
#                                     str(next_peak[1])])
#                     print(msg)
#
#         if debug >= 3:
#             from eqcorrscan.utils import plotting
#             _fname = ''.join(['peaks_',
#                               starttime.datetime.strftime('%Y-%m-%d'),
#                               '.pdf'])
#             print(' '.join(['Saving plot to', _fname]))
#             plotting.peaks_plot(image, starttime, samp_rate, True,
#                                 peaks, _fname)
#         peaks = sorted(peaks, key=lambda time: time[1], reverse=False)
#         return peaks
#     else:
#         print('No peaks for you!')
#         return peaks


def find_peaks2_short(arr, thresh, trig_int, debug=0, starttime=False,
                      samp_rate=1.0):
    r"""Determine peaks in an array of data above a certain \
    threshold. Uses a mask to remove data below threshold and finds peaks in \
    what is left.

    :type arr: ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh: The threshold below which will be considered noise and \
        peaks will not be found in.
    :type trig_int: int
    :param trig_int: The minimum difference in samples between triggers,\
        if multiple peaks within this window this code will find the highest.
    :type debug: int
    :param debug: Optional, debug level 0-5
    :type starttime: osbpy.UTCDateTime
    :param starttime: Starttime for plotting, only used if debug > 2.
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz, only used for plotting if debug > 2.

    :return: peaks: Lists of tuples of peak values and locations.


    >>> import numpy as np
    >>> arr = np.random.randn(100)
    >>> threshold = 10
    >>> arr[40] = 20
    >>> arr[60] = 100
    >>> find_peaks2_short(arr, threshold, 3)
    [(20.0, 40), (100.0, 60)]
    """
    from scipy import ndimage
    import numpy as np
    from obspy import UTCDateTime
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
    # Sort initial peaks according to amplitude
    peaks_sort = sorted(initial_peaks, key=lambda amplitude: amplitude[0],
                        reverse=True)
    # Debugging
    if debug >= 4:
        for peak in initial_peaks:
            print(peak)
    if initial_peaks:
        peaks.append(peaks_sort[0])  # Definitely take the biggest peak
        if debug > 3:
            print(' '.join(['Added the biggest peak of', str(peaks[0][0]),
                            'at sample', str(peaks[0][1])]))
        if len(initial_peaks) > 1:
            if debug > 3:
                msg = ' '.join(['Multiple peaks found, checking',
                                'them now to see if they overlap'])
                print(msg)
            for next_peak in peaks_sort:
                # i in range(1,len(peaks_sort)):
                # Loop through the amplitude sorted peaks
                # if the next highest amplitude peak is within trig_int of any
                # peak already in peaks then we don't want it, else, add it
                # next_peak=peaks_sort[i]
                if debug > 3:
                    print(next_peak)
                for peak in peaks:
                    add = False
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
                    peaks.append(next_peak)
                elif debug > 3:
                    msg = ' '.join(['I did not add peak of', str(next_peak[0]),
                                    'at sample', str(next_peak[1])])
                    print(msg)

        if debug >= 3:
            from eqcorrscan.utils import plotting
            _fname = ''.join(['peaks_',
                              starttime.datetime.strftime('%Y-%m-%d'),
                              '.pdf'])
            plotting.peaks_plot(image, starttime, samp_rate, True,
                                peaks, _fname)
        peaks = sorted(peaks, key=lambda time: time[1], reverse=False)
        return peaks
    else:
        print('No peaks for you!')
        return peaks


def find_peaks_dep(arr, thresh, trig_int, debug=0, starttime=False,
                   samp_rate=1.0):
    r"""Determine peaks in an array of data above a certain \
    threshold: depreciated.

    Depreciated peak-finding routine, very slow, but accurate.  If all else \
    fails this one should work.

    :type arr: ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh: The threshold below which will be considered noise and \
        peaks will not be found in.
    :type trig_int: int
    :param trig_int: The minimum difference in samples between triggers,\
        if multiple peaks within this window this code will find the highest.
    :type starttime: osbpy.UTCDateTime
    :param starttime: Starttime for plotting, only used if debug > 2.
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz, only used for plotting if debug > 2.

    :return: peaks: Lists of tuples of peak values and locations.

    >>> import numpy as np
    >>> arr = np.random.randn(100)
    >>> threshold = 10
    >>> arr[40] = 20
    >>> arr[60] = 100
    >>> find_peaks_dep(arr, threshold, 3)
    [(20.0, 40), (100.0, 60)]
    """
    import numpy as np
    from obspy import UTCDateTime
    if not starttime:
        starttime = UTCDateTime(0)
    # Perform some checks
    if trig_int < 3:
        msg = 'Trigger interval must be greater than 2 samples to find maxima'
        raise IOError(msg)
    # from joblib import Parallel, delayed
    # Will find peaks in the absolute then transfer these to the true values
    sig = np.abs(arr) - thresh
    true_peaks = []
    for i in range(int(trig_int), int(len(sig) - trig_int), int(trig_int)):
        window = sig[i - trig_int: i + trig_int]
        # Define a moving window containing data from +/- the trigger iterval
        peaks = []
        locs = []
        for j in range(1, len(window) - 1):
            # Find all turning points within the window
            if window[j] > 0.0 and window[j] > window[j + 1] and\
               window[j] > window[j - 1]:
                peaks.append(window[j])
                locs.append(i - trig_int + j)
        # Find maximum peak in window
        if peaks:
            true_peaks.append((np.max(np.array(peaks)),
                               locs[np.argmax(np.array(peaks))]))
    # Get unique values
    peaks = sorted(list(set(true_peaks)), key=lambda loc: loc[1])
    # Find highest peak in peaks within trig_int of each other
    for i in range(1, len(peaks) - 1):
        if peaks[i + 1][1] - peaks[i][1] < trig_int:
            if peaks[i][0] < peaks[i + 1][0]:
                peaks[i] = peaks[i + 1]
            else:
                peaks[i + 1] = peaks[i]
        elif peaks[i][1] - peaks[i - 1][1] < trig_int:
            if peaks[i][0] < peaks[i - 1][0]:
                peaks[i] = peaks[i - 1]
            else:
                peaks[i - 1] = peaks[i]
    peaks = sorted(list(set(peaks)), key=lambda loc: loc[1])
    if debug >= 3:
        from eqcorrscan.utils import plotting
        _fname = ''.join(['peaks_',
                          starttime.datetime.strftime('%Y-%m-%d'),
                          '.pdf'])
        plotting.peaks_plot(arr, starttime, samp_rate, True, peaks,
                            _fname)
    peaks = [(peak[0] + thresh, peak[1]) for peak in peaks]
    return peaks


if __name__ == "__main__":
    import doctest
    doctest.testmod()
