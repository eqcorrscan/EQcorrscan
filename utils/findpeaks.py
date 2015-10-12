#!/usr/bin/python
"""
Function to find peaks in data above a certain threshold as part of the EQcorr
package written by Calum Chamberlain of Victoria University of Wellington in
early 2015.

Copyright 2015 Calum Chamberlain

This file is part of EQcorrscan.

    EQcorrscan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EQcorrscan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EQcorrscan.  If not, see <http://www.gnu.org/licenses/>.

"""
from obspy import UTCDateTime
import numpy as np

def is_prime(number):
    """
    Function to test primality of a number. Function lifted from online resource:
        http://www.codeproject.com/Articles/691200/Primality-test-algorithms-Prime-test-The-fastest-w

    This function is distributed under a seperate licence:
        This article, along with any associated source code and files, is
        licensed under The Code Project Open License (CPOL)

    :type number: int
    :param number: Integer to test for primality

    :returns: bool
    """
    import random
    ''' if number != 1 '''
    if (number > 1):
        ''' repeat the test few times '''
        for time in range(3):
            ''' Draw a RANDOM number in range of number ( Z_number )  '''
            randomNumber = random.randint(2, number)-1
            ''' Test if a^(n-1) = 1 mod n '''
            if ( pow(randomNumber, number-1, number) != 1 ):
                return False
        return True
    else:
        ''' case number == 1 '''
        return False


def find_peaks2(arr,thresh, trig_int, debug=0, maxwidth=10,\
                starttime=UTCDateTime('1970-01-01'), samp_rate=1.0):
    """
    Function to determine peaks in an array of data using scipy find_peaks_cwt,
    works fast in certain cases, but for match_filter cccsum peak finding,
    find_peaks2_short works better.  Test it out and see which works best for
    your application.

    :type arr: ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh: The threshold below which will be considered noise and peaks\
    will not be found in.
    :type trig_int: int
    :param trig_int: The minimum difference in samples between triggers,\
    if multiple peaks within this window this code will find the highest.
    :type debug: int
    :param debug: Optional, debug level 0-5
    :type maxwidth: int
    :param maxwidth: Maximum peak width to look for in samples

    :return: peaks, locs: Lists of peak values and locations.

    """
    from scipy.signal import find_peaks_cwt
    # Set everything below the threshold to zero
    image=np.copy(arr)
    image=np.abs(image)
    image[image<thresh]=thresh
    # We need to check if the number of samples in the image is prime, if it
    # is this method will be really slow, so we add a pad to the end to make
    # it not of prime length!
    if is_prime(len(image)):
        image=np.append(image, 0.0)
        print 'Input array has a prime number of samples, appending a zero'
        print len(image)
    if len(image[image>thresh])==0:
        print 'No values over threshold found'
        return []
    if debug > 0:
        print 'Found '+str(len(image[image>thresh]))+' samples above the threshold'
    initial_peaks=[]
    peaks=[]
    # Find the peaks
    print 'Finding peaks'
    peakinds = find_peaks_cwt(image, np.arange(1,maxwidth))
    initial_peaks=[(image[peakind], peakind) for peakind in peakinds]
    # Sort initial peaks according to amplitude
    print 'sorting peaks'
    peaks_sort=sorted(initial_peaks, key=lambda amplitude:amplitude[0],\
                      reverse=True)
    # Debugging
    if debug>=4:
        for peak in initial_peaks:
            print peak
    if initial_peaks:
        peaks.append(peaks_sort[0]) # Definitely take the biggest peak
        if debug > 3:
            print 'Added the biggest peak of '+str(peaks[0][0])+' at sample '+\
                    str(peaks[0][1])
        if len(initial_peaks) > 1:
            if debug>3:
                print 'Multiple peaks found, checking them now to see if they overlap'
            for next_peak in peaks_sort:#i in xrange(1,len(peaks_sort)): # Loop through the amplitude sorted peaks
                # if the next highest amplitude peak is within trig_int of any
                # peak already in peaks then we don't want it, else, add it
                #next_peak=peaks_sort[i]
                if debug>3:
                    print next_peak
                for peak in peaks:
                    add=False # Use add as a switch for whether or not to append
                    # next peak to peaks, if once gone through all the peaks
                    # it is True, then we will add it, otherwise we won't!
                    if abs(next_peak[1]-peak[1]) < trig_int:
                        if debug>3:
                            print 'Difference in time is '+str(next_peak[1]-peak[1])
                            print 'Which is less than '+str(trig_int)
                        add=False
                        # Need to exit the loop here if false
                        break
                    else:
                        add=True
                if add:
                    if debug>3:
                        print 'Adding peak of '+str(next_peak[0])+' at sample '+\
                                str(next_peak[1])
                    peaks.append(next_peak)
                elif debug >3:
                    print 'I did not add peak of '+str(next_peak[0])+\
                            ' at sample '+str(next_peak[1])

        if debug >= 3:
            from utils import EQcorrscan_plotting
            EQcorrscan_plotting.peaks_plot(image, starttime, samp_rate, True, peaks,
                                            'debug_output/peaks_'+\
                                              str(starttime.year)+'-'+\
                                              str(starttime.month)+'-'+\
                                              str(starttime.day)+'.pdf')
        peaks=sorted(peaks, key=lambda time:time[1], reverse=False)
        return peaks
    else:
        print 'No peaks for you!'
        return peaks



def find_peaks2_short(arr,thresh, trig_int, debug=0, \
                starttime=UTCDateTime('1970-01-01'), samp_rate=1.0):
    """
    Function to determine peaks in an array of data above a certain threshold.
    Uses a mask to remove data below threshold and finds peaks in what is left.

    :type arr: ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh: The threshold below which will be considered noise and peaks\
    will not be found in.
    :type trig_int: int
    :param trig_int: The minimum difference in samples between triggers,\
    if multiple peaks within this window this code will find the highest.
    :type debug: int
    :param debug: Optional, debug level 0-5

    :return: peaks, locs: Lists of peak values and locations.

    """
    from scipy import ndimage
    # Set everything below the threshold to zero
    image=np.copy(arr)
    image=np.abs(image)
    image[image<thresh]=0
    if len(image[image>thresh])==0:
        print 'No values over threshold found'
        return []
    if debug > 0:
        print 'Found '+str(len(image[image>thresh]))+' samples above the threshold'
    initial_peaks=[]
    peaks=[]
    # Find the peaks
    labeled_image, number_of_objects= ndimage.label(image)
    peak_slices = ndimage.find_objects(labeled_image)
    for peak_slice in peak_slices:
        #print 'Width of peak='+str(peak_slice[0].stop-peak_slice[0].start)
        window=arr[peak_slice[0].start:peak_slice[0].stop]
        initial_peaks.append((max(window),peak_slice[0].start+np.argmax(window)))
    # Sort initial peaks according to amplitude
    peaks_sort=sorted(initial_peaks, key=lambda amplitude:amplitude[0],\
                      reverse=True)
    # Debugging
    if debug>=4:
        for peak in initial_peaks:
            print peak
    if initial_peaks:
        peaks.append(peaks_sort[0]) # Definitely take the biggest peak
        if debug > 3:
            print 'Added the biggest peak of '+str(peaks[0][0])+' at sample '+\
                    str(peaks[0][1])
        if len(initial_peaks) > 1:
            if debug>3:
                print 'Multiple peaks found, checking them now to see if they overlap'
            for next_peak in peaks_sort:#i in xrange(1,len(peaks_sort)): # Loop through the amplitude sorted peaks
                # if the next highest amplitude peak is within trig_int of any
                # peak already in peaks then we don't want it, else, add it
                #next_peak=peaks_sort[i]
                if debug>3:
                    print next_peak
                for peak in peaks:
                    add=False # Use add as a switch for whether or not to append
                    # next peak to peaks, if once gone through all the peaks
                    # it is True, then we will add it, otherwise we won't!
                    if abs(next_peak[1]-peak[1]) < trig_int:
                        if debug>3:
                            print 'Difference in time is '+str(next_peak[1]-peak[1])
                            print 'Which is less than '+str(trig_int)
                        add=False
                        # Need to exit the loop here if false
                        break
                    else:
                        add=True
                if add:
                    if debug>3:
                        print 'Adding peak of '+str(next_peak[0])+' at sample '+\
                                str(next_peak[1])
                    peaks.append(next_peak)
                elif debug >3:
                    print 'I did not add peak of '+str(next_peak[0])+\
                            ' at sample '+str(next_peak[1])

        if debug >= 3:
            from utils import EQcorrscan_plotting
            EQcorrscan_plotting.peaks_plot(image, starttime, samp_rate, True, peaks,
                                            'debug_output/peaks_'+\
                                              str(starttime.year)+'-'+\
                                              str(starttime.month)+'-'+\
                                              str(starttime.day)+'.pdf')
        peaks=sorted(peaks, key=lambda time:time[1], reverse=False)
        return peaks
    else:
        print 'No peaks for you!'
        return peaks


def find_peaks_dep(arr, thresh, trig_int, debug=0,\
               starttime=UTCDateTime('1970-01-01'), samp_rate=1.0):
    """
    Function to determine peaks in an array of data above a certain threshold.

    Depreciated peak-finding routine, very slow, but accurate.  If all else fails
    this one should work.

    :type arr: ndarray
    :param arr: 1-D numpy array is required
    :type thresh: float
    :param thresh: The threshold below which will be considered noise and peaks\
    will not be found in.
    :type trig_int: int
    :param trig_int: The minimum difference in samples between triggers,\
    if multiple peaks within this window this code will find the highest.

    :return: peaks, locs: Lists of peak values and locations.
    """

    # Perform some checks
    if trig_int < 3:
        import sys
        print 'Trigger interval must be greater than two samples to find maxima'
        sys.exit()
    #from joblib import Parallel, delayed
    # Will find peaks in the absolute then transfer these to the true values
    sig=np.abs(arr)-thresh
    true_peaks=[]
    for i in xrange(int(trig_int),int(len(sig)-trig_int), int(trig_int)):
        window=sig[i-trig_int:i+trig_int] # Define a moving window containing
                                          # data from +/- the trigger iterval
        peaks=[]
        locs=[]
        for j in xrange(1,len(window)-1):
            # Find all turning points within the window
            if window[j] > 0.0 and window[j] > window[j+1] and window[j] > window[j-1]:
                peaks.append(window[j])
                locs.append(i-trig_int+j)
        # Find maximum peak in window
        if peaks:
            true_peaks.append((np.max(np.array(peaks)),\
                               locs[np.argmax(np.array(peaks))]))
    # Get unique values
    peaks=sorted(list(set(true_peaks)), key=lambda loc: loc[1])
    # Find highest peak in peaks within trig_int of each other
    for i in xrange(1,len(peaks)-1):
        if peaks[i+1][1]-peaks[i][1] < trig_int:
            if peaks[i][0] < peaks[i+1][0]:
                peaks[i]=peaks[i+1]
            else:
                peaks[i+1]=peaks[i]
        elif peaks[i][1]-peaks[i-1][1] < trig_int:
            if peaks[i][0] < peaks[i-1][0]:
                peaks[i]=peaks[i-1]
            else:
                peaks[i-1]=peaks[i]
    peaks=sorted(list(set(peaks)), key=lambda loc: loc[1])
    if debug >= 3:
        from utils import EQcorrscan_plotting
        EQcorrscan_plotting.peaks_plot(arr, starttime, samp_rate, True, peaks,
                                        'debug_output/peaks_'+\
                                        str(starttime.year)+'-'+\
                                        str(starttime.month)+'-'+\
                                        str(starttime.day)+'.pdf')
    return peaks
