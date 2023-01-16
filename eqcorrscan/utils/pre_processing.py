"""
Utilities module whose functions are designed to do the basic processing of
the data using obspy modules (which also rely on scipy and numpy).

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import logging
import datetime as dt
import copy

from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count

from obspy import Stream, Trace, UTCDateTime
from obspy.core.trace import Stats
from obspy.signal.filter import bandpass, lowpass, highpass


Logger = logging.getLogger(__name__)


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
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> st = read(TEST_PATH + '/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> _check_daylong(st[0])
    True
    """
    if len(np.nonzero(tr.data)[0]) < 0.5 * len(tr.data):
        qual = False
    else:
        qual = True
    return qual


def shortproc(st, lowcut, highcut, filt_order, samp_rate, parallel=False,
              num_cores=False, starttime=None, endtime=None,
              seisan_chan_names=False, fill_gaps=True, ignore_length=False,
              ignore_bad_data=False, fft_threads=1):
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
    :type parallel: bool
    :param parallel:
        Set to True to process traces in parallel, for small numbers of traces
        this is often slower than serial processing, defaults to False
    :type num_cores: int
    :param num_cores:
        Control the number of cores for parallel processing, if set to False
        then this will use all the cores available.
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime:
        Desired data start time, will trim to this before processing
    :type endtime: obspy.core.utcdatetime.UTCDateTime
    :param endtime:
        Desired data end time, will trim to this before processing
    :type seisan_chan_names: bool
    :param seisan_chan_names:
        Whether channels are named like seisan channels (which are two letters
        rather than SEED convention of three) - defaults to True.
    :type fill_gaps: bool
    :param fill_gaps: Whether to pad any gaps found with zeros or not.
    :type ignore_length: bool
    :param ignore_length:
        Whether to allow data that are less than 80% of the requested length.
        Defaults to False which will error if short data are found.
    :type ignore_bad_data: bool
    :param ignore_bad_data:
        If False (default), errors will be raised if data are excessively
        gappy or are mostly zeros. If True then no error will be raised, but
        an empty trace will be returned.
    :type fft_threads: int
    :param fft_threads:
        Number of threads to use for pyFFTW FFT in resampling. Note that it
        is not recommended to use fft_threads > 1 and num_cores > 1.


    :return: Processed stream
    :rtype: :class:`obspy.core.stream.Stream`

    .. note::
        If your data contain gaps you should *NOT* fill those gaps before
        using the pre-process functions. The pre-process functions will fill
        the gaps internally prior to processing, process the data, then re-fill
        the gaps with zeros to ensure correlations are not incorrectly
        calculated within gaps. If your data have gaps you should pass a merged
        stream without the `fill_value` argument (e.g.: `st = st.merge()`).

    .. warning::
        If you intend to use this for processing templates you should consider
        how resampling will impact your cross-correlations. Minor differences
        in resampling between day-long files (which you are likely to use for
        continuous detection) and shorter files will reduce your
        cross-correlations!

    .. rubric:: Example, bandpass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> st = read(TEST_PATH + '/WAV/TEST_/2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=2, highcut=9, filt_order=3, samp_rate=20,
    ...                parallel=True, num_cores=2)
    >>> print(st[0])
    AF.LABE..SHZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z \
| 20.0 Hz, 1800 samples

    .. rubric:: Example, low-pass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> st = read(TEST_PATH + '/WAV/TEST_/2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=None, highcut=9, filt_order=3,
    ...                samp_rate=20)
    >>> print(st[0])
    AF.LABE..SHZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z \
| 20.0 Hz, 1800 samples

    .. rubric:: Example, high-pass

    >>> from obspy import read
    >>> from eqcorrscan.utils.pre_processing import shortproc
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> st = read(TEST_PATH + '/WAV/TEST_/2013-09-01-0410-35.DFDPC_024_00')
    >>> st = shortproc(st=st, lowcut=2, highcut=None, filt_order=3,
    ...                samp_rate=20)
    >>> print(st[0])
    AF.LABE..SHZ | 2013-09-01T04:10:35.700000Z - 2013-09-01T04:12:05.650000Z \
| 20.0 Hz, 1800 samples
    """
    if isinstance(st, Trace):
        tracein = True
        st = Stream(st)
    else:
        tracein = False
    # Add sanity check for filter
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the nyquist')
    length = None
    clip = False
    if starttime is not None and endtime is not None:
        for tr in st:
            tr.trim(starttime, endtime)
            if len(tr.data) == ((endtime - starttime) *
                                tr.stats.sampling_rate) + 1:
                tr.data = tr.data[1:len(tr.data)]
        length = endtime - starttime
        clip = True
    elif starttime:
        for tr in st:
            tr.trim(starttime=starttime)
    elif endtime:
        for tr in st:
            tr.trim(endtime=endtime)
    for tr in st:
        if len(tr.data) == 0:
            st.remove(tr)
            Logger.warning('No data for {0} after trim'.format(tr.id))
    if parallel:
        if not num_cores:
            num_cores = cpu_count()
        if num_cores > len(st):
            num_cores = len(st)
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(process, (tr,), {
            'lowcut': lowcut, 'highcut': highcut, 'filt_order': filt_order,
            'samp_rate': samp_rate, 'starttime': starttime,
            'clip': clip, 'seisan_chan_names': seisan_chan_names,
            'fill_gaps': fill_gaps, 'length': length,
            'ignore_length': ignore_length, 'fft_threads': fft_threads,
            'ignore_bad_data': ignore_bad_data})
                   for tr in st]
        pool.close()
        try:
            stream_list = [p.get() for p in results]
        except KeyboardInterrupt as e:  # pragma: no cover
            pool.terminate()
            raise e
        pool.join()
        st = Stream(stream_list)
    else:
        for i, tr in enumerate(st):
            st[i] = process(
                tr=tr, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
                samp_rate=samp_rate, starttime=starttime,
                clip=clip, seisan_chan_names=seisan_chan_names,
                fill_gaps=fill_gaps, length=length,
                ignore_length=ignore_length, ignore_bad_data=ignore_bad_data,
                fft_threads=fft_threads)
    if tracein:
        st.merge()
        return st[0]
    return st


def dayproc(st, lowcut, highcut, filt_order, samp_rate, starttime,
            parallel=True, num_cores=False, ignore_length=False,
            seisan_chan_names=False, fill_gaps=True, ignore_bad_data=False,
            fft_threads=1):
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
    :type seisan_chan_names: bool
    :param seisan_chan_names:
        Whether channels are named like seisan channels (which are two letters
        rather than SEED convention of three) - defaults to True.
    :type fill_gaps: bool
    :param fill_gaps: Whether to pad any gaps found with zeros or not.
    :type ignore_bad_data: bool
    :param ignore_bad_data:
        If False (default), errors will be raised if data are excessively
        gappy or are mostly zeros. If True then no error will be raised, but
        an empty trace will be returned.
    :type fft_threads: int
    :param fft_threads:
        Number of threads to use for pyFFTW FFT in resampling. Note that it
        is not recommended to use fft_threads > 1 and num_cores > 1.

    :return: Processed stream.
    :rtype: :class:`obspy.core.stream.Stream`

    .. note::
        If your data contain gaps you should *NOT* fill those gaps before
        using the pre-process functions. The pre-process functions will fill
        the gaps internally prior to processing, process the data, then re-fill
        the gaps with zeros to ensure correlations are not incorrectly
        calculated within gaps. If your data have gaps you should pass a merged
        stream without the `fill_value` argument (e.g.: `st = st.merge()`).

    .. warning::
        Will fail if data are less than 19.2 hours long - this number is
        arbitrary and is chosen to alert the user to the dangers of padding
        to day-long, if you don't care you can ignore this error by setting
        `ignore_length=True`. Use this option at your own risk!  It will also
        warn any-time it has to pad data - if you see strange artifacts in your
        detections, check whether the data have gaps.

    .. rubric:: Example

    >>> import obspy
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...     from obspy.clients.fdsn import Client
    ... else:
    ...     from obspy.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.pre_processing import dayproc
    >>> client = Client('NCEDC')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> bulk_info = [('BP', 'JCNB', '40', 'SP1', t1, t2)]
    >>> st = client.get_waveforms_bulk(bulk_info)
    >>> st_keep = st.copy()  # Copy the stream for later examples
    >>> # Example of bandpass filtering
    >>> st = dayproc(st=st, lowcut=2, highcut=9, filt_order=3, samp_rate=20,
    ...              starttime=t1, parallel=True, num_cores=2)
    >>> print(st[0])
    BP.JCNB.40.SP1 | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.\
950000Z | 20.0 Hz, 1728000 samples
    >>> # Example of lowpass filtering
    >>> st = dayproc(st=st, lowcut=None, highcut=9, filt_order=3, samp_rate=20,
    ...              starttime=t1, parallel=True, num_cores=2)
    >>> print(st[0])
    BP.JCNB.40.SP1 | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.\
950000Z | 20.0 Hz, 1728000 samples
    >>> # Example of highpass filtering
    >>> st = dayproc(st=st, lowcut=2, highcut=None, filt_order=3, samp_rate=20,
    ...              starttime=t1, parallel=True, num_cores=2)
    >>> print(st[0])
    BP.JCNB.40.SP1 | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.\
950000Z | 20.0 Hz, 1728000 samples
    """
    # Add sanity check for filter
    if isinstance(st, Trace):
        st = Stream(st)
        tracein = True
    else:
        tracein = False
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the nyquist')
    # Set the start-time to a day start - cope with
    if starttime is None:
        startdates = []
        for tr in st:
            if abs(tr.stats.starttime - (UTCDateTime(
                    tr.stats.starttime.date) + 86400)) < tr.stats.delta:
                # If the trace starts within 1 sample of the next day, use the
                # next day as the startdate
                startdates.append((tr.stats.starttime + 86400).date)
                Logger.warning(
                    '{0} starts within 1 sample of the next day, using this '
                    'time {1}'.format(
                        tr.id, (tr.stats.starttime + 86400).date))
            else:
                startdates.append(tr.stats.starttime.date)
        # Check that all traces start on the same date...
        if not len(set(startdates)) == 1:
            raise NotImplementedError('Traces start on different days')
        starttime = UTCDateTime(startdates[0])
    if parallel:
        if not num_cores:
            num_cores = cpu_count()
        if num_cores > len(st):
            num_cores = len(st)
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(process, (tr,), {
            'lowcut': lowcut, 'highcut': highcut, 'filt_order': filt_order,
            'samp_rate': samp_rate, 'starttime': starttime, 'clip': True,
            'ignore_length': ignore_length, 'length': 86400,
            'seisan_chan_names': seisan_chan_names, 'fill_gaps': fill_gaps,
            'ignore_bad_data': ignore_bad_data, 'fft_threads': fft_threads})
                   for tr in st]
        pool.close()
        try:
            stream_list = [p.get() for p in results]
        except KeyboardInterrupt as e:  # pragma: no cover
            pool.terminate()
            raise e
        pool.join()
        st = Stream(stream_list)
    else:
        for i, tr in enumerate(st):
            st[i] = process(
                tr=tr, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
                samp_rate=samp_rate, starttime=starttime, clip=True,
                length=86400, ignore_length=ignore_length,
                seisan_chan_names=seisan_chan_names, fill_gaps=fill_gaps,
                ignore_bad_data=ignore_bad_data, fft_threads=fft_threads)
    for tr in st:
        if len(tr.data) == 0:
            st.remove(tr)
    if tracein:
        st.merge()
        return st[0]
    return st


def process(tr, lowcut, highcut, filt_order, samp_rate,
            starttime=False, clip=False, length=86400,
            seisan_chan_names=False, ignore_length=False, fill_gaps=True,
            ignore_bad_data=False, fft_threads=1):
    """
    Basic function to process data, usually called by dayproc or shortproc.

    Functionally, this will bandpass, downsample and check headers and length
    of trace to ensure files start when they should and are the correct length.
    This is a simple wrapper on obspy functions, we include it here to provide
    a system to ensure all parts of the dataset are processed in the same way.

    .. note:: Usually this function is called via dayproc or shortproc.

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to process
    :type lowcut: float
    :param lowcut:
        Low cut in Hz, if set to None and highcut is set, will use
        a lowpass filter.
    :type highcut: float
    :param highcut:
        High cut in Hz, if set to None and lowcut is set, will use
        a highpass filter.
    :type filt_order: int
    :param filt_order: Number of corners for filter.
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz.
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
    :type fill_gaps: bool
    :param fill_gaps: Whether to pad any gaps found with zeros or not.
    :type ignore_bad_data: bool
    :param ignore_bad_data:
        If False (default), errors will be raised if data are excessively
        gappy or are mostly zeros. If True then no error will be raised, but
        an empty trace will be returned.
    :type fft_threads: int
    :param fft_threads: Number of threads to use for pyFFTW FFT in resampling

    :return: Processed trace.
    :type: :class:`obspy.core.stream.Trace`

    .. note::
        If your data contain gaps you should *NOT* fill those gaps before
        using the pre-process functions. The pre-process functions will fill
        the gaps internally prior to processing, process the data, then re-fill
        the gaps with zeros to ensure correlations are not incorrectly
        calculated within gaps. If your data have gaps you should pass a merged
        stream without the `fill_value` argument (e.g.: `tr = tr.merge()`).
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

    Logger.debug('Working on: {0}'.format(tr.id))
    # Check if the trace is gappy and pad if it is.
    gappy = False
    if isinstance(tr.data, np.ma.MaskedArray):
        gappy = True
        gaps, tr = _fill_gaps(tr)
    # Do a brute force quality check
    qual = _check_daylong(tr)
    if not qual:
        msg = ("Data have more zeros than actual data, please check the raw",
               " data set-up and manually sort it: " + tr.stats.station + "." +
               tr.stats.channel)
        if not ignore_bad_data:
            raise ValueError(msg)
        else:
            Logger.warning(msg)
            return Trace(data=np.array([]), header={
                "station": tr.stats.station, "channel": tr.stats.channel,
                "network": tr.stats.network, "location": tr.stats.location,
                "starttime": tr.stats.starttime,
                "sampling_rate": tr.stats.sampling_rate})
    tr = tr.detrend('simple')
    # Detrend data before filtering
    Logger.debug('I have {0} data points for {1} before processing'.format(
        tr.stats.npts, tr.id))

    # Sanity check to ensure files are daylong
    padded = False
    if clip:
        tr = tr.trim(starttime, starttime + length, nearest_sample=True)
    if float(tr.stats.npts / tr.stats.sampling_rate) != length and clip:
        Logger.info(
            'Data for {0} are not long-enough, will zero pad'.format(
                tr.id))
        if tr.stats.endtime - tr.stats.starttime < 0.8 * length\
           and not ignore_length:
            msg = (
                "Data for {0}.{1} is {2:.2f} seconds long, which is less than "
                "80 percent of the desired length ({3} seconds), will not "
                "pad".format(
                    tr.stats.station, tr.stats.channel,
                    tr.stats.endtime - tr.stats.starttime, length))
            if not ignore_bad_data:
                raise NotImplementedError(msg)
            else:
                Logger.warning(msg)
                return Trace(data=np.array([]), header={
                    "station": tr.stats.station, "channel": tr.stats.channel,
                    "network": tr.stats.network, "location": tr.stats.location,
                    "starttime": tr.stats.starttime,
                    "sampling_rate": tr.stats.sampling_rate})
        # trim, then calculate length of any pads required
        pre_pad_secs = tr.stats.starttime - starttime
        post_pad_secs = (starttime + length) - tr.stats.endtime
        if pre_pad_secs > 0 or post_pad_secs > 0:
            padded = True
            pre_pad = np.zeros(int(pre_pad_secs * tr.stats.sampling_rate))
            post_pad = np.zeros(int(post_pad_secs * tr.stats.sampling_rate))
            Logger.debug(str(tr))
            Logger.info("Padding to length with {0} s before and {1} s "
                        "at end".format(pre_pad_secs, post_pad_secs))
            tr.data = np.concatenate([pre_pad, tr.data, post_pad])
            # Use this rather than the expected pad because of rounding samples
            tr.stats.starttime -= len(pre_pad) * tr.stats.delta
            Logger.debug(str(tr))
        # If there is one sample too many after this remove the first one
        # by convention
        if tr.stats.npts == (length * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        # Cope with time precision.
        if abs((tr.stats.sampling_rate * length) -
               tr.stats.npts) > tr.stats.delta:
            msg = ("Data sampling-rate * length ({0} * {1} = {2}) does not "
                   "match number of samples ({3}) for {4}".format(
                    tr.stats.sampling_rate, length,
                    tr.stats.sampling_rate * length, tr.stats.npts, tr.id))
            if not ignore_bad_data:
                raise ValueError(msg)
            else:
                Logger.warning(msg)
                return Trace(data=np.array([]), header={
                    "station": tr.stats.station, "channel": tr.stats.channel,
                    "network": tr.stats.network, "location": tr.stats.location,
                    "starttime": tr.stats.starttime,
                    "sampling_rate": tr.stats.sampling_rate})
        Logger.debug(
            'I now have {0} data points after enforcing length'.format(
                tr.stats.npts))
    # Check sampling rate and resample
    if tr.stats.sampling_rate != samp_rate:
        Logger.debug('Resampling')
        tr = _resample(tr, samp_rate, threads=fft_threads)
    # Filtering section
    tr = tr.detrend('simple')    # Detrend data again before filtering
    if highcut and lowcut:
        Logger.debug('Bandpassing')
        tr.data = bandpass(tr.data, lowcut, highcut,
                           tr.stats.sampling_rate, filt_order, True)
    elif highcut:
        Logger.debug('Lowpassing')
        tr.data = lowpass(tr.data, highcut, tr.stats.sampling_rate,
                          filt_order, True)
    elif lowcut:
        Logger.debug('Highpassing')
        tr.data = highpass(tr.data, lowcut, tr.stats.sampling_rate,
                           filt_order, True)
    else:
        Logger.warning('No filters applied')
    # Account for two letter channel names in s-files and therefore templates
    if seisan_chan_names:
        tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]

    if padded:
        Logger.debug("Reapplying zero pads post processing")
        Logger.debug(str(tr))
        pre_pad = np.zeros(int(pre_pad_secs * tr.stats.sampling_rate))
        post_pad = np.zeros(int(post_pad_secs * tr.stats.sampling_rate))
        pre_pad_len = len(pre_pad)
        post_pad_len = len(post_pad)
        Logger.debug(
            "Taking only valid data between {0} and {1} samples".format(
                pre_pad_len, tr.stats.npts - post_pad_len))
        # Re-apply the pads, taking only the data section that was valid
        tr.data = np.concatenate(
            [pre_pad, tr.data[pre_pad_len: len(tr.data) - post_pad_len],
             post_pad])
        Logger.debug(str(tr))
    # Sanity check to ensure files are correct length
    if float(tr.stats.npts * tr.stats.delta) != length and clip:
        Logger.info(
            'Data for {0} are not of required length, will zero pad'.format(
                tr.id))
        # Use obspy's trim function with zero padding
        tr = tr.trim(starttime, starttime + length, pad=True, fill_value=0,
                     nearest_sample=True)
        # If there is one sample too many after this remove the last one
        # by convention
        if len(tr.data) == (length * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        if abs((tr.stats.sampling_rate * length) -
               tr.stats.npts) > tr.stats.delta:
            raise ValueError('Data are not required length for ' +
                             tr.stats.station + '.' + tr.stats.channel)
    # Replace the gaps with zeros
    if gappy:
        tr = _zero_pad_gaps(tr, gaps, fill_gaps=fill_gaps)
    return tr


def _resample(tr, sampling_rate, threads=1):
    """
    Provide a pyfftw version of obspy's trace resampling.  This code is
    modified from obspy's Trace.resample method.
    """
    from scipy.signal import get_window
    from pyfftw.interfaces.scipy_fftpack import rfft, irfft

    factor = tr.stats.sampling_rate / float(sampling_rate)
    # resample in the frequency domain. Make sure the byteorder is native.
    x = rfft(tr.data.newbyteorder("="), threads=threads)
    # Cast the value to be inserted to the same dtype as the array to avoid
    # issues with numpy rule 'safe'.
    x = np.insert(x, 1, x.dtype.type(0))
    if tr.stats.npts % 2 == 0:
        x = np.append(x, [0])
    x_r = x[::2]
    x_i = x[1::2]

    large_w = np.fft.ifftshift(
        get_window("hann", tr.stats.npts))
    x_r *= large_w[:tr.stats.npts // 2 + 1]
    x_i *= large_w[:tr.stats.npts // 2 + 1]

    # interpolate
    num = int(tr.stats.npts / factor)
    df = 1.0 / (tr.stats.npts * tr.stats.delta)
    d_large_f = 1.0 / num * sampling_rate
    f = df * np.arange(0, tr.stats.npts // 2 + 1, dtype=np.int32)
    n_large_f = num // 2 + 1
    large_f = d_large_f * np.arange(0, n_large_f, dtype=np.int32)
    large_y = np.zeros((2 * n_large_f))
    large_y[::2] = np.interp(large_f, f, x_r)
    large_y[1::2] = np.interp(large_f, f, x_i)

    large_y = np.delete(large_y, 1)
    if num % 2 == 0:
        large_y = np.delete(large_y, -1)
    tr.data = irfft(large_y, threads=threads) * (
            float(num) / float(tr.stats.npts))
    tr.stats.sampling_rate = sampling_rate

    return tr


def _zero_pad_gaps(tr, gaps, fill_gaps=True):
    """
    Replace padded parts of trace with zeros.

    Will cut around gaps, detrend, then pad the gaps with zeros.

    :type tr: :class:`osbpy.core.stream.Trace`
    :param tr: A trace that has had the gaps padded
    :param gaps: List of dict of start-time and end-time as UTCDateTime objects
    :type gaps: list

    :return: :class:`obspy.core.stream.Trace`
    """
    start_in, end_in = (tr.stats.starttime, tr.stats.endtime)
    for gap in gaps:
        stream = Stream()
        if gap['starttime'] > tr.stats.starttime:
            stream += tr.slice(tr.stats.starttime, gap['starttime']).copy()
        if gap['endtime'] < tr.stats.endtime:
            # Note this can happen when gaps are calculated for a trace that
            # is longer than `length`, e.g. gaps are calculated pre-trim.
            stream += tr.slice(gap['endtime'], tr.stats.endtime).copy()
        tr = stream.merge()[0]
    if fill_gaps:
        tr = tr.split()
        tr = tr.detrend()
        tr = tr.merge(fill_value=0)[0]
        # Need to check length - if a gap happened overlapping the end or start
        #  of the trace this will be lost.
        if tr.stats.starttime != start_in:
            # pad with zeros
            tr.data = np.concatenate(
                [np.zeros(int(tr.stats.starttime - start_in)), tr.data])
            tr.stats.starttime = start_in
        if tr.stats.endtime != end_in:
            tr.data = np.concatenate(
                [tr.data, np.zeros(int(end_in - tr.stats.endtime))])
    return tr


def _fill_gaps(tr):
    """
    Interpolate through gaps and work-out where gaps are.

    :param tr: Gappy trace (e.g. tr.data is np.ma.MaskedArray)
    :type tr: `obspy.core.stream.Trace`

    :return: gaps, trace, where gaps is a list of dict
    """
    tr = tr.split()
    gaps = tr.get_gaps()
    tr = tr.detrend().merge(fill_value=0)[0]
    gaps = [{'starttime': gap[4], 'endtime': gap[5]} for gap in gaps]
    return gaps, tr


def _quick_copy_trace(trace, deepcopy_data=True):
    """
    Function to quickly copy a trace. Sets values in the traces' and trace
    header's dict directly, circumventing obspy's init functions.
    Speedup: from 37 us to 12 us per trace - 3x faster

    :type trace: :class:`obspy.core.trace.Trace`
    :param trace: Stream to quickly copy
    :type deepcopy_data: bool
    :param deepcopy_data:
        Whether to deepcopy trace data (with `deepcopy_data=False` expect up to
        20 % speedup, but use only when you know that data trace contents will
        not change or affect results). Warning: do not use this option to copy
        traces with processing history or response information.
    :rtype: :class:`obspy.core.trace.Trace`
    return: trace
    """
    new_trace = Trace()
    for key, value in trace.__dict__.items():
        if key == 'stats':
            new_stats = new_trace.stats
            for key_2, value_2 in value.__dict__.items():
                if isinstance(value_2, UTCDateTime):
                    new_stats.__dict__[key_2] = UTCDateTime(
                        ns=value_2.__dict__['_UTCDateTime__ns'])
                else:
                    new_stats.__dict__[key_2] = value_2
        elif deepcopy_data:
            # data needs to be deepcopied (and anything else, to be safe)
            new_trace.__dict__[key] = copy.deepcopy(value)
        else:  # No deepcopy, e.g. for NaN-traces with no effect on results
            new_trace.__dict__[key] = value
    return new_trace


def _quick_copy_stream(stream, deepcopy_data=True):
    """
    Function to quickly copy a stream.
    Speedup for simple trace:
        from 112 us to 44 (35) us per 3-trace stream - 2.8x (3.2x) faster

    Warning: use `deepcopy_data=False` (saves extra ~20 % time) only when the
             changing the data in the stream later does not change results
             (e.g., for NaN-trace or when data array will not be changed).

    This is what takes longest (1 empty trace, total time to copy 27 us):
    copy header: 18 us (vs create new empty header: 683 ns)
    Two points that can speed up copying / creation:
        1. circumvent trace.__init__ and trace.__set_attr__ by setting value
           directly in trace's __dict__
        2. when setting trace header, circumvent that Stats(header) is called
           when header is already a Stats instance

    :type stream: :class:`obspy.core.stream.Stream`
    :param stream: Stream to quickly copy
    :type deepcopy_data: bool
    :param deepcopy_data:
        Whether to deepcopy data (with `deepcopy_data=False` expect up to 20 %
        speedup, but use only when you know that data trace contents will not
        change or affect results).

    :rtype: :class:`obspy.core.stream.Stream`
    return: stream
    """
    new_traces = list()
    for trace in stream:
        new_traces.append(
            _quick_copy_trace(trace, deepcopy_data=deepcopy_data))
    return Stream(new_traces)


def _stream_quick_select(stream, seed_id):
    """
    4x quicker selection of traces in stream by full Seed-ID. Does not support
    wildcards or selection by network/station/location/channel alone.
    """
    net, sta, loc, chan = seed_id.split('.')
    stream = Stream(
        [tr for tr in stream
         if (tr.stats.network == net and
             tr.stats.station == sta and
             tr.stats.location == loc and
             tr.stats.channel == chan)])
    return stream


def _prep_data_for_correlation(stream, templates, template_names=None,
                               force_stream_epoch=True):
    """
    Check that all channels are the same length and that all channels have data
    for both template and stream.

    Works in place on data - will cut to shortest length

    :param stream: Stream to compare data to
    :param templates:
        List of streams that will be forced to have the same channels as stream
    :param template_names:
        List of strings same length as templates
    :type force_stream_epoch: bool
    :param force_stream_epoch:
        Whether to force all channels in stream to cover the same time period

    :return: stream, templates, template_names (if template_names given)
    """
    n_templates = len(templates)
    template_samp_rates = {
        tr.stats.sampling_rate for template in templates for tr in template}
    stream_samp_rates = {tr.stats.sampling_rate for tr in stream}
    samp_rates = template_samp_rates.union(stream_samp_rates)
    assert len(samp_rates) == 1, "Sampling rates differ"
    samp_rate = samp_rates.pop()

    out_stream = Stream()

    named = True
    if template_names is None:
        named = False
        template_names = range(n_templates)

    # Work out shapes.
    stream_start = min([tr.stats.starttime for tr in stream])
    stream_end = max([tr.stats.endtime for tr in stream])
    if force_stream_epoch:
        stream_length = int(samp_rate * (stream_end - stream_start)) + 1
    else:
        stream_length = max([tr.stats.npts for tr in stream])

    template_length = {
        tr.stats.npts for template in templates for tr in template}
    assert len(template_length) == 1, "Template traces not all the same length"
    template_length = template_length.pop()

    stream_ids = {tr.id for tr in stream}

    # Need to ensure that a channel can be in the template multiple times.
    all_template_ids = [
        Counter([tr.id for tr in template]) for template in templates]
    template_ids = {
        stream_id: max(tid.get(stream_id, 0) for tid in all_template_ids)
        for stream_id in stream_ids}
    template_ids = {_id: value for _id, value in template_ids.items() if value}

    seed_ids = sorted(
        [key.split('.') + [i] for key, value in template_ids.items()
         for i in range(value)])
    seed_ids = [('.'.join(seed_id[0:-1]), seed_id[-1]) for seed_id in seed_ids]

    for channel_number, seed_id in enumerate(template_ids.keys()):
        stream_data = np.zeros(stream_length, dtype=np.float32)
        stream_channel = stream.select(id=seed_id)
        if len(stream_channel) > 1:
            raise NotImplementedError(
                "Multiple channels in continuous data for {0}".format(seed_id))
        stream_channel = stream_channel[0]
        if stream_channel.stats.npts == stream_length:
            stream_data = stream_channel.data
        else:
            Logger.info('Data for {0} is not as long as needed, '
                        'padding'.format(stream_channel.id))
            if force_stream_epoch:
                start_pad = int(samp_rate * (
                        stream_channel.stats.starttime - stream_start))
                end_pad = stream_length - (
                        start_pad + stream_channel.stats.npts)
                # In some cases there will be one sample missing when sampling
                # time-stamps are not set consistently between channels, this
                # results in start_pad and end_pad being len==0
                if start_pad == 0 and end_pad == 0:
                    Logger.debug("Start and end pad are both zero, padding "
                                 "at one end")
                    if (stream_channel.stats.starttime - stream_start) > (
                       stream_end - stream_channel.stats.endtime):
                        start_pad = int(
                            stream_length - stream_channel.stats.npts)
                    else:
                        end_pad = int(
                            stream_length - stream_channel.stats.npts)
                stream_channel.stats.starttime -= (start_pad / samp_rate)
            else:
                start_pad = 0
                end_pad = stream_length - stream_channel.stats.npts
            if end_pad == 0:
                stream_data[start_pad:] = stream_channel.data
            else:
                stream_data[start_pad:-end_pad] = stream_channel.data
        header = stream_channel.stats.copy()
        header.npts = stream_length
        out_stream += Trace(data=stream_data, header=header)

    # Initialize nan template for speed.
    nan_channel = np.full(template_length, np.nan, dtype=np.float32)
    nan_channel = np.require(nan_channel, requirements=['C_CONTIGUOUS'])
    nan_template = Stream()
    for _seed_id in seed_ids:
        net, sta, loc, chan = _seed_id[0].split('.')
        nan_template += Trace(header=Stats({
            'network': net, 'station': sta, 'location': loc,
            'channel': chan, 'starttime': UTCDateTime(ns=0),
            'npts': template_length, 'sampling_rate': samp_rate}))

    # Remove templates with no matching channels
    filt = np.ones(len(template_names)).astype(bool)
    for i, template in enumerate(templates):
        trace_ids = {tr.id for tr in template}
        if len(trace_ids.intersection(stream_ids)) == 0:
            filt[i] = 0

    _out = dict(zip(
        [_tn for _tn, _filt in zip(template_names, filt) if _filt],
        [_t for _t, _filt in zip(templates, filt) if _filt]))
    flt_templates = list(_out.values())

    if len(_out) != len(templates):
        Logger.debug("Some templates not used due to no matching channels")

    # Ensure that the templates' earliest traces are kept, even if there is no
    # continuous data for them. If this happens, we need to add a NaN-stream to
    # the continuous data to avoid inconsistent detection times.
    n_template_traces = np.array([len(temp) for temp in flt_templates])
    n_stream_traces = sum([n+1 for s, n in seed_ids])
    # These checks are not necessary if all templates will get NaN-traces,
    # because the NaN-traces will save the right starttime for the template.
    nan_stream_ids = list()
    if any(n_template_traces > n_stream_traces):
        earliest_templ_trace_ids = set(
            [template.sort(['starttime'])[0].id for template in flt_templates])
        for earliest_templ_trace_id in earliest_templ_trace_ids:
            if earliest_templ_trace_id not in template_ids:
                nan_stream_ids.append(earliest_templ_trace_id)
                net, sta, loc, chan = earliest_templ_trace_id.split('.')
                nan_template += Trace(header=Stats({
                    'network': net, 'station': sta, 'location': loc,
                    'channel': chan, 'starttime': UTCDateTime(ns=0),
                    'sampling_rate': samp_rate}))
                stream_nan_data = np.full(
                    stream_length, np.nan, dtype=np.float32)
                out_stream += Trace(
                    data=np.ma.masked_array(stream_nan_data, stream_nan_data),
                    header=Stats({
                        'network': net, 'station': sta, 'location': loc,
                        'channel': chan, 'starttime': stream_start,
                        'npts': stream_length, 'sampling_rate': samp_rate}))
                seed_ids.append((earliest_templ_trace_id, 0))

    incomplete_templates = {
        template_name for template_name, template in _out.items() if
        sorted([tr.id for tr in template]) != [tr.id for tr in nan_template]}

    # Fill out the templates with nan channels
    for template_name in incomplete_templates:
        template = _out[template_name]
        template_starttime = min(tr.stats.starttime for tr in template)
        out_template = _quick_copy_stream(nan_template, deepcopy_data=False)

        # Select traces very quickly: assume that trace order does not change,
        # make dict of trace-ids and list of indices and use indices to select
        stream_trace_id_dict = defaultdict(list)
        for n, tr in enumerate(template.traces):
            stream_trace_id_dict[tr.id].append(n)

        for channel_number, _seed_id in enumerate(seed_ids):
            seed_id, channel_index = _seed_id
            # Select all traces with same seed_id, based on indices for
            # corresponding traces stored in stream_trace_id_dict
            # Much quicker than: template_channel = template.select(id=seed_id)
            template_channel = Stream([
                template.traces[idx] for idx in stream_trace_id_dict[seed_id]])
            if len(template_channel) <= channel_index:
                # out_template[channel_number].data = nan_channel  # quicker:
                out_template[channel_number].__dict__['data'] = copy.deepcopy(
                    nan_channel)
                out_template[channel_number].stats.__dict__['npts'] = \
                    template_length
                out_template[channel_number].stats.__dict__['starttime'] = \
                    template_starttime
                out_template[channel_number].stats.__dict__['endtime'] = \
                    UTCDateTime(ns=int(
                        round(template_starttime.ns
                              + (template_length * samp_rate) * 1e9)))
            else:
                out_template[channel_number] = template_channel[channel_index]
        # If a template-trace matches a NaN-trace in the stream , then set
        # template-trace to NaN so that this trace does not appear in channel-
        # list of detections.
        if len(nan_stream_ids) > 0:
            for tr in out_template:
                if tr.id in nan_stream_ids:
                    tr.data = nan_channel
        _out.update({template_name: out_template})

    out_templates = list(_out.values())
    out_template_names = list(_out.keys())

    if named:
        return out_stream, out_templates, out_template_names
    return out_stream, out_templates


if __name__ == "__main__":
    import doctest
    doctest.testmod()
