#!/usr/bin/python
"""
Utilities module for the EQcorrscan package written by Calum Chamberlain of
Victoria University Wlelington.  These functions are designed to do the basic
processing of the data using obspy modules (which also rely on scipy and numpy).
"""

# There is a problem with Aaron's processing of FRAN, whereby channels
# are not named consistently, horizontal channels are intermitently
# referred to as N,E and 1,2 after the borehole sensor has been deployed
# and these come from the borehole.  This was worked around in matlab,
# and while Laura-May is working on the matlab codes I will impliment
# a work aroudn here, but this should be corrected by a bulk process
# of the raw data to miniseed, which will then serve as the go-to
# resource.


from obspy import read, UTCDateTime
from obspy.signal.filter import bandpass

def shortproc(st, lowcut, highcut, filt_order, samp_rate, debug):
    """
    Basic function to bandpass, downsample.  Works in place
    on data.  This is employed to ensure all parts of the data are processed
    in the same way.

    :type st: obspy.Stream
    :type highcut: float
    :type lowcut: float
    :type filt_order: int
    :type samp_rate: float
    :type debug: int
    :type starttime: obspy.UTCDateTime

    :return: obspy.Stream
    """

    for tr in st:
        if debug > 4:
            tr.plot()
        # Check sampling rate and resample
        if tr.stats.sampling_rate != samp_rate:
            tr.resample(samp_rate)

        # Filtering section
        tr=tr.detrend('simple')    # Detrend data before filtering
        tr.data=bandpass(tr.data, lowcut, highcut,
                    tr.stats.sampling_rate, filt_order, True)
        # Correct FRAN N,E channels to 1,2 as this is correct
        if tr.stats.station=='FRAN' and tr.stats.channel=='SHN':
            print 'Correcting FRAN SHN to SH1'
            tr.stats.channel='SH1'
        if tr.stats.station=='FRAN' and tr.stats.channel=='SHE':
            print 'Correcting FRAN SHE to SH2'
            tr.stats.channel='SH2'
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
    :type highcut: float
    :type lowcut: float
    :type filt_order: int
    :type samp_rate: float
    :type debug: int
    :type starttime: obspy.UTCDateTime

    :return: obspy.Stream
    """

    day=str(starttime.year)+str(starttime.month).zfill(2)+\
        str(starttime.day).zfill(2)
    if debug>=2:
        print 'Working on: '+tr.stats.station+'.'+tr.stats.channel
    if debug >= 4:
        tr.plot()
    # Check sampling rate and resample
    if tr.stats.sampling_rate != samp_rate:
        if debug>=2:
            print 'Resampling'
        tr.resample(samp_rate)

    # If there is one sample too many remove the first sample - this occurs
    # at station FOZ where the first sample is zero when it shouldn't be,
    # Not real sample generated during data download
    if len(tr.data)==(86400*tr.stats.sampling_rate)+1:
        tr.data=tr.data[1:len(tr.data)]

    # Filtering section
    tr=tr.detrend('simple')    # Detrend data before filtering
    if debug>=2:
        print 'Bandpassing'
    tr.data=bandpass(tr.data, lowcut, highcut,
                tr.stats.sampling_rate, filt_order, True)
    # Correct FRAN N,E channels to 1,2 as this is correct
    if tr.stats.station=='FRAN' and tr.stats.channel=='SHN':
        print 'Correcting FRAN SHN to SH1'
        tr.stats.channel='SH1'
    if tr.stats.station=='FRAN' and tr.stats.channel=='SHE':
        print 'Correcting FRAN SHE to SH2'
        tr.stats.channel='SH2'

    # Correct FOZ channels
    if tr.stats.station=='FOZ':
        tr.stats.channel='HH'+tr.stats.channel[2]
    # Account for two letter channel names in s-files and therefore templates
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


