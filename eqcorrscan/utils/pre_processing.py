#!/usr/bin/python
"""
Utilities module for the EQcorrscan package written by Calum Chamberlain of
Victoria University Wlelington.  These functions are designed to do the basic
processing of the data using obspy modules (which also rely on scipy and numpy).

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
from obspy.signal.filter import bandpass

def _check_daylong(tr):
    """
    Function to check the data quality of the daylong file - check to see that
    the day isn't just zeros, with large steps, if it is then the resampling will
    hate it.

    :type tr: obspy.Trace
    :param tr: Trace to check if the data are daylong.

    :return qual: bool
    """
    import numpy as np
    if len(tr.data)-len(np.nonzero(tr.data)) < 0.5*len(tr.data):
        qual=False
    else:
        qual=True
    return qual

def shortproc(st, lowcut, highcut, filt_order, samp_rate, debug=0):
    """
    Basic function to bandpass, downsample.  Works in place
    on data.  This is employed to ensure all parts of the data are processed
    in the same way.

    :type st: obspy.Stream
    :param st: Stream to process
    :type highcut: float
    :param highcut: High cut for bandpass in Hz
    :type lowcut: float
    :param lowcut: Low cut for bandpass in Hz
    :type filt_order: int
    :param filt_order: Number of corners for bandpass filter
    :type samp_rate: float
    :param samp_rate: Sampling rate desired in Hz
    :type debug: int
    :param debug: Debug flag from 0-5, higher numbers = more output

    :return: obspy.Stream
    """

    for tr in st:
        if debug > 4:
            tr.plot()
        # Check data quality first
        qual=_check_daylong(tr)
        if not qual:
            msg="Data have more zeros than actual data, please check the raw" +\
                    "data set-up and manually sort it"
            raise ValueError(msg)
        # Check sampling rate and resample
        if tr.stats.sampling_rate != samp_rate:
            tr.resample(samp_rate)

        # Filtering section
        tr=tr.detrend('simple')    # Detrend data before filtering
        tr.data=bandpass(tr.data, lowcut, highcut,
                    tr.stats.sampling_rate, filt_order, True)
        # Correct FRAN N,E channels to 1,2 as this is correct
        # if tr.stats.station=='FRAN' and tr.stats.channel=='SHN':
            # print 'Correcting FRAN SHN to SH1'
            # tr.stats.channel='SH1'
        # if tr.stats.station=='FRAN' and tr.stats.channel=='SHE':
            # print 'Correcting FRAN SHE to SH2'
            # tr.stats.channel='SH2'
        # Account for two letter channel names in s-files and therefore templates
        tr.stats.channel=tr.stats.channel[0]+tr.stats.channel[2]

        # Final visual check for debug
        if debug > 4:
            tr.plot()
    return st

def dayproc(tr, lowcut, highcut, filt_order, samp_rate, debug, starttime):
    """
    Basic function to bandpass, downsample and check headers and length of trace
    to ensure files start at the start of a day and are daylong.  Works in place
    on data.  This is employed to ensure all parts of the data are processed
    in the same way.

    :type tr: obspy.Trace
    :param tr: Trace to process
    :type highcut: float
    :param highcut: High cut in Hz for bandpass
    :type lowcut: float
    :type lowcut: Low cut in Hz for bandpass
    :type filt_order: int
    :param filt_order: Corners for bandpass
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type debug: int
    :param debug: Debug output level from 0-5, higher numbers = more output
    :type starttime: obspy.UTCDateTime
    :param starttime: Desired start of trace

    :return: obspy.Stream
    """

    day=str(starttime.year)+str(starttime.month).zfill(2)+\
        str(starttime.day).zfill(2)
    if debug>=2:
        print 'Working on: '+tr.stats.station+'.'+tr.stats.channel
    if debug >= 5:
        tr.plot()
    # Do a brute force quality check
    qual=_check_daylong(tr)
    if not qual:
        msg="Data have more zeros than actual data, please check the raw" +\
                "data set-up and manually sort it"
        raise ValueError(msg)
    tr=tr.detrend('simple')    # Detrend data before filtering

    # If there is one sample too many remove the first sample - this occurs
    # at station FOZ where the first sample is zero when it shouldn't be,
    # Not real sample generated during data download
    if len(tr.data)==(86400*tr.stats.sampling_rate)+1:
        tr.data=tr.data[1:len(tr.data)]
    print 'I have '+str(len(tr.data))+' data points for '+tr.stats.station+'.'\
            +tr.stats.channel+' before processing'

    # Sanity check to ensure files are daylong
    if float(tr.stats.npts/tr.stats.sampling_rate) != 86400.0:
        if debug >= 2:
            print 'Data for '+tr.stats.station+'.'+tr.stats.channel+\
                    ' is not of daylong length, will zero pad'
        # Work out when the trace thinks it is starting - Aaron's time headers
        # are often wrong
        traceday=UTCDateTime(str(tr.stats.starttime.year)+'-'+\
                            str(tr.stats.starttime.month)+'-'+\
                            str(tr.stats.starttime.day))
        # Use obspy's trim function with zero padding
        tr=tr.trim(traceday,traceday+86400,pad=True,fill_value=0,\
                    nearest_sample=True)
        # If there is one sample too many after this remove the last one
        # by convention
        if len(tr.data)==(86400*tr.stats.sampling_rate)+1:
            tr.data=tr.data[1:len(tr.data)]
        if not tr.stats.sampling_rate*86400 == tr.stats.npts:
                raise ValueError ('Data are not daylong for '+tr.stats.station+\
                                  '.'+tr.stats.channel)

    print 'I now have '+str(len(tr.data))+' data points after enforcing day length'

    # Check sampling rate and resample
    if tr.stats.sampling_rate != samp_rate:
        if debug>=2:
            print 'Resampling'
        tr.resample(samp_rate)


    # Filtering section
    tr=tr.detrend('simple')    # Detrend data before filtering
    if debug>=2:
        print 'Bandpassing'
    tr.data=bandpass(tr.data, lowcut, highcut,
                tr.stats.sampling_rate, filt_order, True)

    ################### CHANNEL renaming, does not need to be in the main codes!
    # Correct FRAN N,E channels to 1,2 as this is correct
    if tr.stats.station=='FRAN' and tr.stats.channel=='SHN':
        print 'Correcting FRAN SHN to SH1'
        tr.stats.channel='SH1'
    if tr.stats.station=='FRAN' and tr.stats.channel=='SHE':
        print 'Correcting FRAN SHE to SH2'
        tr.stats.channel='SH2'

    # Correct FOZ channels
    if tr.stats.station=='FOZ':
        if len(tr.stats.channel)==3:
            tr.stats.channel='HH'+tr.stats.channel[2]
        else:
            tr.stats.channel='HH'+tr.stats.channel[1]
    # Account for two letter channel names in s-files and therefore templates
    if len(tr.stats.channel)==3:
        tr.stats.channel=tr.stats.channel[0]+tr.stats.channel[2]

    # Sanity check the time header
    if str(tr.stats.starttime.year)+str(tr.stats.starttime.month).zfill(2)+\
        str(tr.stats.starttime.day).zfill(2) != day:
        if debug >= 2:
            print "Time headers are wrong: "+str(tr.stats.starttime)
            print "Correcting to: "+str(starttime)
            tr.stats.starttime=starttime

    # Sanity check to ensure files are daylong
    if float(tr.stats.npts/tr.stats.sampling_rate) != 86400.0:
        if debug >= 2:
            print 'Data for '+tr.stats.station+'.'+tr.stats.channel+\
                    ' is not of daylong length, will zero pad'
        # Work out when the trace thinks it is starting - Aaron's time headers
        # are often wrong
        traceday=UTCDateTime(str(tr.stats.starttime.year)+'-'+\
                            str(tr.stats.starttime.month)+'-'+\
                            str(tr.stats.starttime.day))
        # Use obspy's trim function with zero padding
        tr=tr.trim(traceday,traceday+86400,pad=True,fill_value=0,\
                    nearest_sample=True)
        # If there is one sample too many after this remove the last one
        # by convention
        if len(tr.data)==(86400*tr.stats.sampling_rate)+1:
            tr.data=tr.data[1:len(tr.data)]
        if not tr.stats.sampling_rate*86400 == tr.stats.npts:
                raise ValueError ('Data are not daylong for '+tr.stats.station+\
                                  '.'+tr.stats.channel)
    # Final visual check for debug
    if debug >= 4:
        tr.plot()
    return tr
