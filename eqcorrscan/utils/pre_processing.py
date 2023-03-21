"""
Utilities module whose functions are designed to do the basic processing of
the data using obspy modules (which also rely on scipy and numpy).

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import os

import numpy as np
import logging
import copy

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache, partial
from scipy.signal import iirfilter, sosfilt, zpk2sos

from obspy import Stream, Trace, UTCDateTime
from obspy.core.trace import Stats


Logger = logging.getLogger(__name__)


def _check_daylong(data, threshold=0.5):
    """
    Check data continuity.

    Check to see that the day is more than threshold of zeros, if it is
    then the resampling will hate it.

    :type data: np.array
    :param data: Data from Trace to check if the data are okay.
    :type threshold: float
    :param threshold: Fraction of data to accept as zeros.

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
    >>> _check_daylong(st[0].data)
    True
    >>> zeroed_data = st[0].copy().data
    >>> zeroed_data[0:9100] = np.zeros(9100)
    >>> _check_daylong(zeroed_data)
    False
    """
    return np.nonzero(data)[0].shape[0] >= threshold * data.shape[0]


def _simple_qc(st, max_workers=None, chunksize=1):
    qual = dict()
    with ThreadPoolExecutor(max_workers) as executor:
        for tr, _qual in zip(st, executor.map(
                _check_daylong, (tr.data for tr in st), chunksize=chunksize)):
            qual[tr.id] = _qual
    return qual


def _sanitize_length(st, starttime=None, endtime=None, daylong=False):
    """

    :param st:
    :param starttime:
    :param endtime:
    :param daylong:
    :return:
    """
    length, clip = None, False

    if daylong:
        length, clip = 86400, True
        # Set the start-time to a day start - cope with
        if starttime is None:
            startdates = []
            for tr in st:
                if abs(tr.stats.starttime - (UTCDateTime(
                        tr.stats.starttime.date) + 86400)) < tr.stats.delta:
                    # If the trace starts within 1 sample of the next day,
                    # use the next day as the startdate
                    startdates.append((tr.stats.starttime + 86400).date)
                    Logger.warning(
                        f'{tr.id} starts within 1 sample of the next day, '
                        f'using this time {(tr.stats.starttime + 86400).date}')
                else:
                    startdates.append(tr.stats.starttime.date)
            # Check that all traces start on the same date...
            if not len(set(startdates)) == 1:
                raise NotImplementedError('Traces start on different days')
            starttime = UTCDateTime(startdates[0])
    else:
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
    return st, length, clip, starttime


@lru_cache
def _get_window(window, npts):
    from scipy.signal import get_window
    return np.fft.ifftshift(get_window(window, npts))


def multi_process(st, lowcut, highcut, filt_order, samp_rate, parallel=False,
                  num_cores=False, starttime=None, endtime=None,
                  daylong=False, seisan_chan_names=False, fill_gaps=True,
                  ignore_length=False, ignore_bad_data=False):
    """

    :param st:
    :param lowcut:
    :param highcut:
    :param filt_order:
    :param samp_rate:
    :param parallel:
    :param num_cores:
    :param starttime:
    :param endtime:
    :param daylong:
    :param seisan_chan_names:
    :param fill_gaps:
    :param ignore_length:
    :param ignore_bad_data:
    :return:
    """
    if isinstance(st, Trace):
        tracein = True
        st = Stream(st)
    else:
        tracein = False
    # Add sanity check for filter
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the Nyquist')
    if highcut and lowcut:
        assert lowcut < highcut, f"Lowcut: {lowcut} above highcut: {highcut}"

    # Allow datetimes for starttime and endtime
    if starttime and not isinstance(starttime, UTCDateTime):
        starttime = UTCDateTime(starttime)
    if starttime is False:
        starttime = None
    if endtime and not isinstance(endtime, UTCDateTime):
        endtime = UTCDateTime(endtime)
    if endtime is False:
        endtime = None

    # Make sensible choices about workers and chunk sizes
    if parallel:
        if not num_cores:
            # We don't want to over-specify threads, we don't have IO
            # bound tasks
            max_workers = min(len(st), os.cpu_count())
        else:
            max_workers = min(len(st), num_cores)
    else:
        max_workers = 1
    chunksize = len(st) // max_workers

    st, length, clip, starttime = _sanitize_length(
        st=st, starttime=starttime, endtime=endtime, daylong=daylong)

    for tr in st:
        if len(tr.data) == 0:
            st.remove(tr)
            Logger.warning('No data for {0} after trim'.format(tr.id))

    # Do work
    # 1. Fill gaps and keep track of them
    gappy = {tr.id: False for tr in st}
    gaps = dict()
    for i, tr in enumerate(st):
        if isinstance(tr.data, np.ma.MaskedArray):
            gappy[tr.id] = True
            gaps[tr.id], tr = _fill_gaps(tr)
            st[i] = tr

    # 2. Check for zeros and cope with bad data
    ## ~ 4x speedup for 50 100 Hz daylong traces on 12 threads
    qual = _simple_qc(st, max_workers=max_workers, chunksize=chunksize)
    for trace_id, _qual in qual.items():
        if not _qual:
            msg = ("Data have more zeros than actual data, please check the "
                   f"raw data set-up and manually sort it: {tr.id}")
            if not ignore_bad_data:
                raise ValueError(msg)
            else:
                # Remove bad traces from the stream
                st.remove(st.select(id=trace_id))

    # 3. Detrend
    ## ~ 2x speedup for 50 100 Hz daylong traces on 12 threads
    st = _multi_detrend(st, max_workers=max_workers, chunksize=chunksize)

    # 4. Check length and pad to length
    padded = {tr.id: (0., 0.) for tr in st}
    if clip:
        st.trim(starttime, starttime + length, nearest_sample=True)
        # Indexing because we are going to overwrite traces
        for i, _ in enumerate(st):
            if float(st[i].stats.npts / st[i].stats.sampling_rate) != length:
                Logger.info(
                    'Data for {0} are not long-enough, will zero pad'.format(
                        st[i].id))
                st[i], padded[st[i].id] = _length_check(
                    st[i], starttime=starttime, length=length,
                    ignore_length=ignore_length,
                    ignore_bad_data=ignore_bad_data)
        # Remove None traces that might be returned from length checking
        st.traces = [tr for tr in st if tr is not None]

    # Check that we actually still have some data
    if not _stream_has_data(st):
        if tracein:
            return st[0]
        return st

    # 5. Resample
    ## ~ 3.25x speedup for 50 100 Hz daylong traces on 12 threads
    st = _multi_resample(
        st, sampling_rate=samp_rate, max_workers=max_workers,
        chunksize=chunksize)

    # Detrend again before filtering
    st = _multi_detrend(st, max_workers=max_workers, chunksize=chunksize)

    # 6. Filter
    ## ~3.25x speedup for 50 100 Hz daylong traces on 12 threads
    st = _multi_filter(
        st, highcut=highcut, lowcut=lowcut, filt_order=filt_order,
        max_workers=max_workers, chunksize=chunksize)

    # 7. Reapply zeros after processing from 4
    for tr in st:
        # Pads default to (0., 0.), pads should only ever be positive.
        if sum(padded[tr.id]) == 0:
            continue
        Logger.debug("Reapplying zero pads post processing")
        Logger.debug(str(tr))
        pre_pad = np.zeros(int(padded[tr.id][0] * tr.stats.sampling_rate))
        post_pad = np.zeros(int(padded[tr.id][1] * tr.stats.sampling_rate))
        pre_pad_len = len(pre_pad)
        post_pad_len = len(post_pad)
        Logger.debug(
            f"Taking only valid data between {pre_pad_len} and "
            f"{tr.stats.npts - post_pad_len} samples")
        # Re-apply the pads, taking only the data section that was valid
        tr.data = np.concatenate(
            [pre_pad, tr.data[pre_pad_len: len(tr.data) - post_pad_len],
             post_pad])
        Logger.debug(str(tr))

    # 8. Recheck length
    for tr in st:
        if float(tr.stats.npts * tr.stats.delta) != length and clip:
            Logger.info(f'Data for {tr.id} are not of required length, will '
                        f'zero pad')
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

    # 9. Re-insert gaps from 1
    # TODO: If this is slow it could probably be threaded
    for i, tr in enumerate(st):
        if gappy[tr.id]:
            st[i] = _zero_pad_gaps(tr, gaps[tr.id], fill_gaps=fill_gaps)

    # 10. Clean up
    for tr in st:
        if len(tr.data) == 0:
            st.remove(tr)

    # 11. Account for seisan channel naming
    if seisan_chan_names:
        for tr in st:
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]

    if tracein:
        st.merge()
        return st[0]
    return st


def _empty_trace(tr):
    """
    Generate an empty trace with a basic header matching the input trace
    """
    bad_trace = Trace(
        data=np.array([]), header={
            "station": tr.stats.station, "channel": tr.stats.channel,
            "network": tr.stats.network, "location": tr.stats.location,
            "starttime": tr.stats.starttime,
            "sampling_rate": tr.stats.sampling_rate})
    return bad_trace


def _stream_has_data(st):
    return sum(tr.stats.npts for tr in st) > 0


def _length_check(tr, starttime, length, ignore_length, ignore_bad_data):
    trace_length = tr.stats.endtime - tr.stats.starttime
    if trace_length < 0.8 * length and not ignore_length:
        msg = f"Data for {tr.id} is {trace_length:.2f} seconds "\
              f"long, which is less than 80 percent of the desired "\
              f"length ({length} seconds), will not pad"
        if not ignore_bad_data:
            raise NotImplementedError(msg)
        else:
            Logger.warning(msg)
            return _empty_trace(tr), (0., 0.)
    # trim, then calculate length of any pads required
    pre_pad_secs = tr.stats.starttime - starttime
    post_pad_secs = (starttime + length) - tr.stats.endtime
    if pre_pad_secs > 0 or post_pad_secs > 0:
        pre_pad = np.zeros(int(pre_pad_secs * tr.stats.sampling_rate))
        post_pad = np.zeros(int(post_pad_secs * tr.stats.sampling_rate))
        Logger.debug(str(tr))
        Logger.info(
            f"Padding to length with {pre_pad_secs} s before "
            f"and {post_pad} s at end")
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
        msg = (f"Data sampling-rate * length ({tr.stats.sampling_rate} *"
               f" {length} = {tr.stats.sampling_rate * length}) does not "
               f"match number of samples ({tr.stats.npts}) for {tr.id}")
        if not ignore_bad_data:
            raise ValueError(msg)
        else:
            Logger.warning(msg)
            return _empty_trace(tr), (0., 0.)
    Logger.debug(
        f'I now have {tr.stats.npts} data points after enforcing length')
    return tr, (pre_pad_secs, post_pad_secs)


def _multi_filter(st, highcut, lowcut, filt_order, max_workers=None,
                  chunksize=1):
    """ Zero-phase filtering of multi-channel data. """
    if not highcut and not lowcut:
        Logger.warning("No filters applied")
        return st
    # Require that all channels are the same sampling frequency
    samp_rate = set(tr.stats.sampling_rate for tr in st)
    assert len(samp_rate) == 1, "Different sampling rates found"
    samp_rate = samp_rate.pop()
    # Sanity check filter bounds
    if highcut:
        assert highcut * 2 < samp_rate, "Highcut must be below Nyquist"
    if highcut and lowcut:
        assert lowcut < highcut, "Lowcut must be below highcut"

    fe = 0.5 * samp_rate
    if lowcut:
        low = lowcut / fe
    if highcut:
        high = highcut / fe

    # Design filter
    if highcut and lowcut:
        z, p, k = iirfilter(
            filt_order, [low, high], btype='band',
            ftype='butter', output='zpk')
    elif highcut:
        z, p, k = iirfilter(
            filt_order, high, btype='lowpass', ftype='butter',
            output='zpk')
    elif lowcut:
        z, p, k = iirfilter(
            filt_order, low, btype='highpass', ftype='butter',
            output='zpk')

    sos = zpk2sos(z, p, k)

    _filter = partial(_zerophase_filter, sos)

    with ThreadPoolExecutor(max_workers) as executor:
        results = executor.map(
            _filter, (tr.data for tr in st), chunksize=chunksize)

    for r, tr in zip(results, st):
        tr.data = r

    return st


def _zerophase_filter(sos, data):
    firstpass = sosfilt(sos, data)
    return sosfilt(sos, firstpass[::-1])[::-1]


def _multi_detrend(st, max_workers=None, chunksize=1):
    """
    Multithreaded detrending using simple linear detrend between
    first and last values. Follows obspy "simple" detrend.

    :param st:
    :return:
    """
    for tr in st:
        tr.data = np.require(tr.data, np.float64)
    with ThreadPoolExecutor(max_workers) as executor:
        results = executor.map(_detrend, (tr.data for tr in st),
                               chunksize=chunksize)
    # Ensure tasks complete
    _ = (r for r in results)
    return st


def _detrend(data):
    """
    Detrend signal simply by subtracting a line through the first and last
    point of the trace

    :param data: Data to detrend, type numpy.ndarray.
    :return: Nothing - works in place
    """
    ndat = data.shape[0]
    x1, x2 = data[0], data[-1]
    data -= x1 + np.arange(ndat) * (x2 - x1) / np.float64(ndat - 1)
    return


def _multi_resample(st, sampling_rate, max_workers=None, chunksize=1):
    """

    :param st:
    :param sampling_rate:
    :param max_workers:
    :return:
    """
    # Get the windows, and downsampling factors ahead of time
    to_resample = (
        (tr.data, tr.stats.delta,
         tr.stats.sampling_rate / float(sampling_rate),
         sampling_rate, _get_window("hann", tr.stats.npts))
        for tr in st)
    with ThreadPoolExecutor(max_workers) as executor:
        # Unpack tuple using lamba
        results = executor.map(lambda args: _resample(*args), to_resample,
                               chunksize=chunksize)
    for r, tr in zip(results, st):
        tr.data = r
        tr.stats.sampling_rate = sampling_rate
    return st


def _resample(data, delta, factor, sampling_rate, large_w):
    """
    Resample data in the frequency domain - adapted from obspy resample

    """
    if factor == 1:
        # No resampling needed, don't waste time.
        return data
    # Need to work with numpy objects to release the GIL
    npts = data.shape[0]
    df = np.float32(1.0) / (npts * delta)
    num = np.int32(npts / factor)
    d_large_f = np.float32(1.0) / num * sampling_rate

    # Forward fft
    x = np.fft.rfft(data)
    # Window
    x *= large_w[:npts // 2 + 1]

    # interpolate
    f = df * np.arange(0, npts // 2 + 1, dtype=np.int32)
    n_large_f = num // 2 + 1
    large_f = d_large_f * np.arange(0, n_large_f, dtype=np.int32)

    # Have to split into real and imaginary parts for interpolation.
    y = np.interp(large_f, f, np.real(x)) + (1j * np.interp(
        large_f, f, np.imag(x)))
    # Try to reduce memory before doing the ifft
    del large_f, f, x

    return np.fft.irfft(y, n=num)[0:num] * (np.float32(num) / np.float32(npts))


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
                out_template.traces[channel_number].__dict__[
                    'data'] = np.copy(nan_channel)
                out_template.traces[channel_number].stats.__dict__[
                    'npts'] = template_length
                out_template.traces[channel_number].stats.__dict__[
                    'starttime'] = template_starttime
                out_template.traces[channel_number].stats.__dict__[
                    'endtime'] = UTCDateTime(ns=int(
                        round(template_starttime.ns
                              + (template_length / samp_rate) * 1e9)))
            else:
                out_template.traces[channel_number] = template_channel.traces[
                    channel_index]

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


def shortproc(st, lowcut, highcut, filt_order, samp_rate, parallel=False,
              num_cores=False, starttime=None, endtime=None,
              seisan_chan_names=False, fill_gaps=True, ignore_length=False,
              ignore_bad_data=False, fft_threads=1):
    """
    Deprecated
    """
    Logger.warning("Shortproc is depreciated after 0.4.4 and will "
                   "be removed in a future version. Use multi_process"
                   " instead")
    st = multi_process(
        st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
        samp_rate=samp_rate, parallel=parallel, num_cores=num_cores,
        starttime=starttime, endtime=endtime, daylong=False,
        seisan_chan_names=seisan_chan_names, fill_gaps=fill_gaps,
        ignore_length=ignore_length, ignore_bad_data=ignore_bad_data)
    return st


def dayproc(st, lowcut, highcut, filt_order, samp_rate, starttime,
            parallel=True, num_cores=False, ignore_length=False,
            seisan_chan_names=False, fill_gaps=True, ignore_bad_data=False,
            fft_threads=1):
    """
    Deprecated
    """
    Logger.warning("dayproc is depreciated after 0.4.4 and will be "
                   "removed in a future version. Use multi_process instead")
    st = multi_process(
        st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
        samp_rate=samp_rate, parallel=parallel, num_cores=num_cores,
        starttime=starttime, endtime=None, daylong=True,
        seisan_chan_names=seisan_chan_names, fill_gaps=fill_gaps,
        ignore_length=ignore_length, ignore_bad_data=ignore_bad_data)
    return st


def process(tr, lowcut, highcut, filt_order, samp_rate,
            starttime=False, clip=False, length=86400,
            seisan_chan_names=False, ignore_length=False, fill_gaps=True,
            ignore_bad_data=False, fft_threads=1):
    """
    Deprecated
    """
    Logger.warning("process is depreciated after 0.4.4 and will be removed "
                   "in a future version. Use multi_process instead")
    if length == 86400:
        daylong = True
    else:
        daylong = False

    endtime = None
    if clip:
        if not starttime:
            starttime = tr.stats.starttime
        elif not isinstance(starttime, UTCDateTime):
            starttime = UTCDateTime(starttime)
        endtime = starttime + length
    st = multi_process(
        st=tr, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
        samp_rate=samp_rate, parallel=False, num_cores=1,
        starttime=starttime, endtime=endtime, daylong=daylong,
        seisan_chan_names=seisan_chan_names, fill_gaps=fill_gaps,
        ignore_length=ignore_length, ignore_bad_data=ignore_bad_data)
    return st


if __name__ == "__main__":
    import doctest
    doctest.testmod()
