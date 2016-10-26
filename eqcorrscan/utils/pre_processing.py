"""
Utilities module whose functions are designed to do the basic \
processing of the data using obspy modules (which also rely on scipy and \
numpy).

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
import warnings
import datetime as dt

from multiprocessing import Pool, cpu_count

from obspy import Stream, Trace, UTCDateTime
from obspy.signal.filter import bandpass, lowpass, highpass


def _check_daylong(tr):
    """
    Check the data quality of the daylong file.

    Check to see that the day isn't just zeros, with large steps, if it is
    then the resampling will hate it.

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to check if the data are daylong.

    :return quality (simply good or bad)
    :rtype: bool

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import _check_daylong
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> _check_daylong(st[0])
    True
    """
    if len(tr.data) - len(np.nonzero(tr.data)) < 0.5 * len(tr.data):
        qual = False
    else:
        qual = True
    return qual


def shortproc(st, lowcut, highcut, filt_order, samp_rate, debug=0,
              parallel=False, num_cores=False, starttime=None, endtime=None):
    """
    Basic function to bandpass and downsample.

    Works in place on data.  This is employed to ensure all parts of the
    data are processed in the same way.

    :type st: obspy.core.stream.Stream
    :param st: Stream to process
    :type lowcut: float
    :param lowcut: Low cut for bandpass in Hz
    :type highcut: float
    :param highcut: High cut for bandpass in Hz
    :type filt_order: int
    :param filt_order: Number of corners for bandpass filter
    :type samp_rate: float
    :param samp_rate: Sampling rate desired in Hz
    :type debug: int
    :param debug: Debug flag from 0-5, higher numbers = more output
    :type parallel: bool
    :param parallel: Set to True to process traces in parallel, for small \
        numbers of traces this is often slower than serial processing, \
        defaults to False
    :type num_cores: int
    :param num_cores: Control the number of cores for parallel processing, \
        if set to False then this will use all the cores.
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime:
        Desired data start time, will trim to this before processing
    :type endtime: obspy.core.utcdatetime.UTCDateTime
    :param endtime:
        Desired data end time, will trim to this before processing

    :return: Processed stream
    :rtype: :class:`obspy.core.stream.Stream`

    .. note:: Will convert channel names to two characters long.

    .. warning::
        If you intend to use this for processing templates you should consider
        how resampling will impact your cross-correlations. Minor differences
        in resampling between day-long files (which you are likely to use for
        continuous detection) and shorter files will reduce your
        cross-correlations!

    .. rubric:: Example, bandpass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=2, highcut=9, filt_order=3, samp_rate=20,
    ...                debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    AF.LABE..SZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z |\
 20.0 Hz, 1800 samples

    .. rubric:: Example, low-pass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=None, highcut=9, filt_order=3,
    ...                samp_rate=20, debug=0)
    >>> print(st[0])
    AF.LABE..SZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z |\
 20.0 Hz, 1800 samples

    .. rubric:: Example, high-pass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=2, highcut=None, filt_order=3,
    ...                samp_rate=20, debug=0)
    >>> print(st[0])
    AF.LABE..SZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z |\
 20.0 Hz, 1800 samples
    """
    if isinstance(st, Trace):
        tracein = True
        st = Stream(st)
    else:
        tracein = False
    # Add sanity check for filter
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the nyquist')
    if debug > 4:
        parallel = False
    if starttime and endtime:
        for tr in st:
            tr.trim(starttime, endtime)
            print(len(tr))
            if len(tr.data) == ((endtime - starttime) *
                                tr.stats.sampling_rate) + 1:
                tr.data = tr.data[1:len(tr.data)]
    elif starttime:
        for tr in st:
            tr.trim(starttime=starttime)
    elif endtime:
        for tr in st:
            tr.trim(endtime=endtime)
    if parallel:
        if not num_cores:
            num_cores = cpu_count()
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(process, (tr,), {'lowcut': lowcut,
                                                     'highcut': highcut,
                                                     'filt_order': filt_order,
                                                     'samp_rate': samp_rate,
                                                     'debug': debug,
                                                     'starttime': False,
                                                     'clip': False})
                   for tr in st]
        pool.close()
        stream_list = [p.get() for p in results]
        pool.join()
        st = Stream(stream_list)
    else:
        for tr in st:
            process(tr=tr, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, samp_rate=samp_rate, debug=debug,
                    starttime=False, clip=False)
    if tracein:
        st.merge()
        return st[0]
    return st


def dayproc(st, lowcut, highcut, filt_order, samp_rate, starttime, debug=0,
            parallel=True, num_cores=False, ignore_length=False):
    """
    Wrapper for dayproc to parallel multiple traces in a stream.

    Works in place on data.  This is employed to ensure all parts of the data \
    are processed in the same way.

    :type st: obspy.core.stream.Stream
    :param st: Stream to process (can be trace).
    :type lowcut: float
    :param lowcut: Low cut in Hz for bandpass.
    :type highcut: float
    :param highcut: High cut in Hz for bandpass.
    :type filt_order: int
    :param filt_order: Corners for bandpass.
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz.
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Desired start-date of trace.
    :type debug: int
    :param debug: Debug output level from 0-5, higher numbers = more output.
    :type parallel: bool
    :param parallel:
        Set to True to process traces in parallel, this is often faster than
        serial processing of traces: defaults to True.
    :type num_cores: int
    :param num_cores:
        Control the number of cores for parallel processing, if set to False
        then this will use all the cores.
    :type ignore_length: bool
    :param ignore_length: See warning below.

    :return: Processed stream.
    :rtype: :class:`obspy.core.stream.Stream`

    .. note:: Will convert channel names to two characters long.

    .. warning::
        Will fail if data are less than 19.2 hours long - this number is
        arbitrary and is chosen to alert the user to the dangers of padding
        to day-long, if you don't care you can ignore this error by setting
        `ignore_length=True`. Use this option at your own risk!  It will also
        warn any-time it has to pad data - if you see strange artifacts in your
        detections, check whether the data have gaps.

    .. rubric:: Example, bandpass

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('GEONET')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('NZ', 'FOZ', '10', 'HHE', t1, t2),
    ...              ('NZ', 'FOZ', '10', 'HHE', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st = dayproc(st=st, lowcut=2, highcut=9, filt_order=3, samp_rate=20,
    ...              starttime=t1, debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    NZ.FOZ.10.HE | 2012-03-25T23:59:59.998393Z - 2012-03-26T23:59:59.948393Z |\
 20.0 Hz, 1728000 samples


    .. rubric:: Example, low-pass

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('GEONET')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('NZ', 'FOZ', '10', 'HHE', t1, t2),
    ...              ('NZ', 'FOZ', '10', 'HHE', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st = dayproc(st=st, lowcut=None, highcut=9, filt_order=3, samp_rate=20,
    ...              starttime=t1, debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    NZ.FOZ.10.HE | 2012-03-25T23:59:59.998393Z - 2012-03-26T23:59:59.948393Z |\
 20.0 Hz, 1728000 samples

    .. rubric:: Example, high-pass

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('GEONET')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('NZ', 'FOZ', '10', 'HHE', t1, t2),
    ...              ('NZ', 'FOZ', '10', 'HHE', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st = dayproc(st=st, lowcut=2, highcut=None, filt_order=3, samp_rate=20,
    ...              starttime=t1, debug=0, parallel=True, num_cores=2)
    >>> print(st[0])
    NZ.FOZ.10.HE | 2012-03-25T23:59:59.998393Z - 2012-03-26T23:59:59.948393Z |\
 20.0 Hz, 1728000 samples
    """
    # Add sanity check for filter
    if isinstance(st, Trace):
        st = Stream(st)
        tracein = True
    else:
        tracein = False
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the nyquist')
    if debug > 4:
        parallel = False
    if parallel:
        if not num_cores:
            num_cores = cpu_count()
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(process, (tr,), {'lowcut': lowcut,
                                                     'highcut': highcut,
                                                     'filt_order': filt_order,
                                                     'samp_rate': samp_rate,
                                                     'debug': debug,
                                                     'starttime': starttime,
                                                     'clip': True,
                                                     'ignore_length':
                                                         ignore_length,
                                                     'length': 86400})
                   for tr in st]
        pool.close()
        stream_list = [p.get() for p in results]
        pool.join()
        st = Stream(stream_list)
    else:
        for tr in st:
            process(tr=tr, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, samp_rate=samp_rate, debug=debug,
                    starttime=starttime, clip=True,
                    ignore_length=ignore_length, length=86400)
    if tracein:
        st.merge()
        return st[0]
    return st


def process(tr, lowcut, highcut, filt_order, samp_rate, debug,
            starttime=False, clip=False, length=86400,
            seisan_chan_names=True, ignore_length=False):
    """
    Basic function to process data, usually called by dayproc or shortproc.

    Functionally, this will bandpass, downsample and check headers and length
    of trace to ensure files start at the start of a day and are daylong.
    This is a simple wrapper on obspy functions, we include it here to provide
    a system to ensure all parts of the dataset are processed in the same way.

    .. note:: Usually this function is called via dayproc or shortproc.

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to process
    :type lowcut: float
    :param lowcut: Low cut in Hz, if set to None and highcut is set, will use \
        a lowpass filter.
    :type highcut: float
    :param highcut: High cut in Hz, if set to None and lowcut is set, will \
        use a highpass filter.
    :type filt_order: int
    :param filt_order: Number of corners for filter.
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz.
    :type debug: int
    :param debug: Debug output level from 0-5, higher numbers = more output.
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Desired start of trace
    :type clip: bool
    :param clip: Whether to expect, and enforce a set length of data or not.
    :type length: float
    :param length: Use to set a fixed length for data from the given starttime.
    :type seisan_chan_names: bool
    :param seisan_chan_names:
        Whether channels are named like seisan channels (which are two letters
        rather than SEED convention of three) - defaults to True.
    :type ignore_length: bool
    :param ignore_length: See warning in dayproc.

    :return: Processed stream.
    :type: :class:`obspy.core.stream.Stream`
    """
    # Add sanity check
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the nyquist')

    # Define the start-time
    if starttime:
        # Be nice and allow a datetime object.
        if isinstance(starttime, dt.date) or isinstance(starttime,
                                                        dt.datetime):
            starttime = UTCDateTime(starttime)
        day = starttime.date
    else:
        day = tr.stats.starttime.date

    if debug >= 2:
        print('Working on: ' + tr.stats.station + '.' + tr.stats.channel)
    if debug >= 5:
        tr.plot()
    # Do a brute force quality check
    qual = _check_daylong(tr)
    if not qual:
        msg = ("Data have more zeros than actual data, please check the raw",
               " data set-up and manually sort it: " + tr.stats.station + "." +
               tr.stats.channel)
        raise ValueError(msg)
    tr = tr.detrend('simple')
    # Detrend data before filtering
    if debug > 0:
        print('I have ' + str(len(tr.data)) + ' data points for ' +
              tr.stats.station + '.' + tr.stats.channel +
              ' before processing')

    # Sanity check to ensure files are daylong
    if float(tr.stats.npts / tr.stats.sampling_rate) != length and clip:
        if debug >= 2:
            print('Data for ' + tr.stats.station + '.' + tr.stats.channel +
                  ' are not of daylong length, will zero pad')
        if tr.stats.endtime - tr.stats.starttime < 0.8 * length\
           and not ignore_length:
            msg = ('Data for %s.%s is %i hours long, which is less than 0.8 '
                   'of the desired length, will not pad' %
                   (tr.stats.station, tr.stats.channel,
                    (tr.stats.endtime - tr.stats.starttime) / 3600))
            raise NotImplementedError(msg)
        # Use obspy's trim function with zero padding
        tr = tr.trim(starttime, starttime + length, pad=True, fill_value=0,
                     nearest_sample=True)
        # If there is one sample too many after this remove the first one
        # by convention
        if len(tr.data) == (length * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        if not tr.stats.sampling_rate * length == tr.stats.npts:
                raise ValueError('Data are not daylong for ' +
                                 tr.stats.station + '.' + tr.stats.channel)

        print('I now have %i data points after enforcing length'
              % len(tr.data))
    # Check sampling rate and resample
    if tr.stats.sampling_rate != samp_rate:
        if debug >= 2:
            print('Resampling')
        tr.resample(samp_rate)
    # Filtering section
    tr = tr.detrend('simple')    # Detrend data again before filtering
    if highcut and lowcut:
        if debug >= 2:
            print('Bandpassing')
        tr.data = bandpass(tr.data, lowcut, highcut,
                           tr.stats.sampling_rate, filt_order, True)
    elif highcut:
        if debug >= 2:
            print('Lowpassing')
        tr.data = lowpass(tr.data, highcut, tr.stats.sampling_rate,
                          filt_order, True)
    elif lowcut:
        if debug >= 2:
            print('Highpassing')
        tr.data = highpass(tr.data, lowcut, tr.stats.sampling_rate,
                           filt_order, True)
    else:
        warnings.warn('No filters applied')
    # Account for two letter channel names in s-files and therefore templates
    if seisan_chan_names:
        tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]

    # Sanity check the time header
    if tr.stats.starttime.day != day and clip:
        warnings.warn("Time headers do not match expected date: " +
                      str(tr.stats.starttime))

    # Sanity check to ensure files are daylong
    if float(tr.stats.npts / tr.stats.sampling_rate) != length and clip:
        if debug >= 2:
            print('Data for ' + tr.stats.station + '.' + tr.stats.channel +
                  ' is not of daylong length, will zero pad')
        # Use obspy's trim function with zero padding
        tr = tr.trim(starttime, starttime + length, pad=True, fill_value=0,
                     nearest_sample=True)
        # If there is one sample too many after this remove the last one
        # by convention
        if len(tr.data) == (length * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        if not tr.stats.sampling_rate * length == tr.stats.npts:
                raise ValueError('Data are not daylong for ' +
                                 tr.stats.station + '.' + tr.stats.channel)
    # Final visual check for debug
    if debug > 4:
        tr.plot()
    return tr


if __name__ == "__main__":
    import doctest
    doctest.testmod()
