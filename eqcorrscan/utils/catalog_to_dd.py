"""
Functions to generate hypoDD input files from catalogs.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import os
import numpy as np
import logging
from collections import namedtuple, defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

from obspy import UTCDateTime, Stream
from obspy.core.event import (
    Catalog, Event, Origin, Magnitude, Pick, WaveformStreamID, Arrival,
    OriginQuality)
from eqcorrscan.utils.clustering import dist_mat_km
from eqcorrscan.core.lag_calc import _concatenate_and_correlate, _xcorr_interp


Logger = logging.getLogger(__name__)

SeedPickID = namedtuple("SeedPickID", "seed_id phase_hint")


# Generic helpers

class _DTObs(object):
    """ Holder for phase observations """
    def __init__(self, station, tt1, tt2, weight, phase):
        self.station = station
        assert len(self.station) <= 7, "Station must be five characters or less"
        self.tt1 = tt1
        self.tt2 = tt2
        self.weight = weight
        self.phase = phase
        assert self.phase in "PS", "Only P or S phases accepted"

    @property
    def ct_string(self):
        return "{sta:<7s} {tt1:7.3f} {tt2:7.3f} {weight:6.4f} {phase}".format(
            sta=self.station, tt1=self.tt1, tt2=self.tt2, weight=self.weight,
            phase=self.phase)

    @property
    def cc_string(self):
        return "{sta:<7s} {dt:7.3f} {weight:6.4f} {phase}".format(
            sta=self.station, dt=self.tt1 - self.tt2, weight=self.weight,
            phase=self.phase)


class _EventPair(object):
    """ Holder for event paid observations. """
    def __init__(self, event_id_1, event_id_2, obs=None):
        self.event_id_1 = event_id_1
        self.event_id_2 = event_id_2
        self.obs = obs or list()

    @property
    def _header(self):
        return "# {event_id_1:>9d} {event_id_2:>9d}".format(
            event_id_1=self.event_id_1, event_id_2=self.event_id_2)

    @property
    def ct_string(self):
        parts = [self._header]
        parts.extend([obs.ct_string for obs in self.obs])
        return '\n'.join(parts)

    @property
    def cc_string(self):
        parts = [self._header + " 0.0"]  # No origin time correction
        parts.extend([obs.cc_string for obs in self.obs])
        return '\n'.join(parts)


def _generate_event_id_mapper(catalog, event_id_mapper=None):
    """
    Generate or fill an event_id_mapper mapping event resource id to integer id
    """
    event_id_mapper = event_id_mapper or dict()
    try:
        largest_event_id = max(event_id_mapper.values())
    except ValueError:
        largest_event_id = 0
    for event in catalog:
        if event.resource_id.id not in event_id_mapper.keys():
            event_id = largest_event_id + 1
            largest_event_id = event_id
            event_id_mapper.update({event.resource_id.id: event_id})
    return event_id_mapper


def _get_arrival_for_pick(event, pick):
    """
    Get a matching arrival for a given pick.
    """
    matched_arrivals = [
        arr for origin in event.origins for arr in origin.arrivals
        if arr.pick_id.get_referred_object() == pick]
    if len(matched_arrivals) == 0:  # pragma: no cover
        return None
    if len(matched_arrivals) > 1:  # pragma: no cover
        Logger.warning("Multiple arrivals found for pick: {0} on {1}".format(
            pick.phase_hint, pick.waveform_id.get_seed_string()))
    return matched_arrivals[0]


def _combined_weight(arr1, arr2):
    """
    Combine the weights between two arrivals to give a weight between 0-1.

    Uses the mean of the time_weight attributes of the arrivals
    """
    if arr1:
        w1 = arr1.time_weight or 1.0
    else:
        w1 = 1.0
    if arr2:
        w2 = arr2.time_weight or 1.0
    else:
        w2 = 1.0
    return (w1 + w2) / 2.0


def _prepare_stream(stream, event, extract_len, pre_pick, seed_pick_ids=None):
    """
    Slice stream around picks

    returns a dictionary of traces keyed by phase_hint.
    """
    seed_pick_ids = seed_pick_ids or {
        SeedPickID(pick.waveform_id.get_seed_string(), pick.phase_hint[0])
        for pick in event.picks}
    stream_sliced = defaultdict(lambda: Stream())
    for seed_pick_id in seed_pick_ids:
        pick = [pick for pick in event.picks
                if pick.waveform_id.get_seed_string() == seed_pick_id.seed_id
                and pick.phase_hint[0] == seed_pick_id.phase_hint]
        if len(pick) > 1:
            Logger.warning(
                "Multiple picks for {seed_id}, phase-hint {phase_hint}, using "
                "the earliest".format(
                    seed_id=seed_pick_id.seed_id,
                    phase_hint=seed_pick_id.phase_hint))
            pick = sorted(pick, key=lambda p: p.time)
        elif len(pick) == 0:
            continue
        pick = pick[0]
        tr = stream.select(id=seed_pick_id.seed_id).slice(
            starttime=pick.time - pre_pick,
            endtime=(pick.time - pre_pick) + extract_len).merge()
        if len(tr) == 0:
            continue
        if len(tr) > 1:
            Logger.error("Multiple traces for {seed_id}".format(
                seed_id=seed_pick_id.seed_id))
            continue
        stream_sliced.update(
            {seed_pick_id.phase_hint:
                 stream_sliced[seed_pick_id.phase_hint] + tr[0]})
    return stream_sliced


# Time calculators

def _compute_dt_correlations(catalog, master, min_link, event_id_mapper,
                             stream_dict, min_cc, extract_len, pre_pick,
                             shift_len, interpolate):
    differential_times_dict = dict()
    master_stream = _prepare_stream(
        stream=stream_dict[master.resource_id.id], event=master,
        extract_len=extract_len, pre_pick=pre_pick)
    master_seed_ids = {
        SeedPickID(pick.waveform_id.get_seed_string(), pick.phase_hint[0])
        for pick in master.picks if pick.phase_hint[0] in "PS"}
    # Dictionary of travel-times for master keyed by {station}_{phase_hint}
    master_tts = dict()
    master_origin_time = (master.preferred_origin() or master.origins[0]).time
    for pick in master.picks:
        if pick.phase_hint[0] not in "PS":
            continue
        tt1 = pick.time - master_origin_time
        master_tts.update({
            "{0}_{1}".format(
                pick.waveform_id.station_code, pick.phase_hint[0]): tt1})

    matched_length = extract_len + (2 * shift_len)
    matched_pre_pick = pre_pick + shift_len
    # We will use this to maintain order
    event_dict = {event.resource_id.id: event for event in catalog}
    event_ids = list(event_dict.keys())
    matched_streams = {
        event_id: _prepare_stream(
            stream=stream_dict[event_id], event=event_dict[event_id],
            extract_len=matched_length, pre_pick=matched_pre_pick,
            seed_pick_ids=master_seed_ids)
        for event_id in event_ids}

    sampling_rates = {tr.stats.sampling_rate for st in master_stream.values()
                      for tr in st}
    for phase_hint in master_stream.keys():  # Loop over P and S separately
        for sampling_rate in sampling_rates:   # Loop over separate samp rates
            delta = 1.0 / sampling_rate
            _master_stream = master_stream[phase_hint].select(
                sampling_rate=sampling_rate)
            _matched_streams = [
                matched_streams[event.resource_id.id][phase_hint].select(
                    sampling_rate=sampling_rate) for event in catalog]
            # Check lengths
            master_length = Counter(
                (tr.stats.npts for tr in _master_stream)).most_common(1)[0][0]
            _master_stream = _master_stream.select(npts=master_length)
            matched_length = Counter(
                (tr.stats.npts for st in _matched_streams
                 for tr in st)).most_common(1)[0][0]
            # Remove empty streams and generate an ordered list of event_ids
            used_event_ids, used_matched_streams = [], []
            for event_id, _matched_stream in zip(event_ids, _matched_streams):
                _matched_stream = _matched_stream.select(npts=matched_length)
                if len(_matched_stream) > 0:
                    used_event_ids.append(event_id)
                    used_matched_streams.append(_matched_stream)
            ccc_out, used_chans = _concatenate_and_correlate(
                template=_master_stream, streams=used_matched_streams, cores=1)
            # Convert ccc_out to pick-time
            for i, used_event_id in enumerate(used_event_ids):
                for j, chan in used_chans[i]:
                    if not chan.used:
                        continue
                    correlation = ccc_out[i][j]
                    if interpolate:
                        shift, cc_max = _xcorr_interp(correlation, dt=delta)
                    else:
                        cc_max = np.amax(correlation)
                        shift = np.argmax(correlation) * delta
                    if cc_max < min_cc:
                        continue
                    pick = [p for p in event_dict[used_event_id].picks
                            if p.phase_hint == phase_hint
                            and p.waveform_id.station_code == chan.channel[0]
                            and p.waveform_id.channel_code == chan.channel[1]]
                    pick = sorted(pick, key=lambda p: p.time)[0]
                    tt2 = pick.time - (
                        event_dict[used_event_id].preferred_origin() or
                        event_dict[used_event_id].origins[0]).time
                    tt2 += shift
                    diff_time = differential_times_dict.get(
                        used_event_id, None)
                    if diff_time is None:
                        diff_time = _EventPair(
                            event_id_1=event_id_mapper[master.resource_id.id],
                            event_id_2=event_id_mapper[used_event_id])
                    diff_time.obs.append(
                        _DTObs(station=used_chans[0],
                               tt1=master_tts["{0}_{1}".format(
                                   used_chans[0], phase_hint)],
                               tt2=tt2, weight=cc_max ** 2, phase=phase_hint))
    # Threshold on min_link
    differential_times = [dt for dt in differential_times_dict.values()
                          if len(dt.obs) >= min_link]
    return differential_times


def _compute_dt(catalog, master, min_link, event_id_mapper):
    """
    Inner function to compute differential times between a catalog and a
    master event.
    """
    if len(catalog) == 0:
        return []
    differential_times = [
        _make_event_pair(event=event, master=master,
                         event_id_mapper=event_id_mapper) for event in catalog
        if event.resource_id != master.resource_id]
    return [d for d in differential_times if len(d.obs) >= min_link]


def _make_event_pair(event, master, event_id_mapper):
    """
    Make an event pair for a given event and master event.
    """
    differential_times = _EventPair(
        event_id_1=event_id_mapper[master.resource_id.id],
        event_id_2=event_id_mapper[event.resource_id.id])
    for master_pick in master.picks:
        if master_pick.phase_hint not in "PS":
            continue
        master_arr = _get_arrival_for_pick(master, master_pick)
        matched_picks = [p for p in event.picks
                         if p.waveform_id == master_pick.waveform_id
                         and p.phase_hint == master_pick.phase_hint]
        master_moveout = master_pick.time - (
                master.preferred_origin() or master.origins[0]).time
        for matched_pick in matched_picks:
            matched_arr = _get_arrival_for_pick(event, matched_pick)
            matched_moveout = matched_pick.time - (
                event.preferred_origin() or event.origins[0]).time
            differential_times.obs.append(
                _DTObs(station=master_pick.waveform_id.station_code,
                       tt1=master_moveout, tt2=matched_moveout,
                       weight=_combined_weight(master_arr, matched_arr),
                       phase=master_pick.phase_hint))
    return differential_times


def compute_differential_times(catalog, correlation, stream_dict=None,
                               event_id_mapper=None, max_sep=8., min_link=8,
                               min_cc=None, extract_len=None, pre_pick=None,
                               shift_len=None, interpolate=False,
                               max_workers=None):
    """
    Generate groups of differential times for a catalog.

    :type catalog: `obspy.core.event.Catalog`
    :param catalog: Catalog of events to get differential times for
    :type correlation: bool
    :param correlation:
        If True will generate cross-correlation derived differential-times for
        a dt.cc file. If false, will generate catalog times for a dt.ct file.
    :type stream_dict: dict
    :param stream_dict:
        Dictionary of streams keyed by event-id (the event.resource_id.id,
        NOT the hypoDD event-id)
    :type event_id_mapper: dict
    :param event_id_mapper:
        Dictionary mapping event resource id to an integer event id for hypoDD.
        If this is None, or missing events then the dictionary will be updated
        to include appropriate event-ids. This should be of the form
        {event.resource_id.id: integer_id}
    :type max_sep: float
    :param max_sep: Maximum hypocentral separation in km to link events
    :type min_link: int
    :param min_link: Minimum shared phase observations to link events
    :type min_cc: float
    :param min_cc: Threshold to include cross-correlation results.
    :type extract_len: float
    :param extract_len: Length in seconds to extract around the pick
    :type pre_pick: float
    :param pre_pick: Time before the pick to start the correlation window
    :type shift_len: float
    :param shift_len: Time to allow pick to vary in seconds
    :type interpolate: bool
    :param interpolate:
        Whether to interpolate correlations or not. Allows subsample accuracy
    :type max_workers: int
    :param max_workers:
        Maximum number of workers for parallel processing. If None then all
        threads will be used.

    :rtype: dict
    :return: Dictionary of differential times keyed by event id.
    :rtype: dict
    :return: Dictionary of event_id_mapper

    .. note::
        The arguments min_cc, stream_dict, extract_len, pre_pick, shift_len
        and interpolate are only required if correlation=True.
    """
    correlation_kwargs = dict(
        min_cc=min_cc, stream_dict=stream_dict, extract_len=extract_len,
        pre_pick=pre_pick, shift_len=shift_len, interpolate=interpolate)
    if correlation:
        for arg, name in correlation_kwargs.items():
            assert arg is not None, "{0} is required for correlation".format(
                name)
    max_workers = max_workers or cpu_count()
    # Ensure all events have locations and picks.
    event_id_mapper = _generate_event_id_mapper(
        catalog=catalog, event_id_mapper=event_id_mapper)
    distances = dist_mat_km(catalog)
    distance_filter = distances <= max_sep

    differential_times = {}
    sub_catalogs = [[ev for i, ev in enumerate(catalog)
                     if master_filter[i]] for master_filter in distance_filter]

    additional_args = dict(min_link=min_link, event_id_mapper=event_id_mapper)
    if correlation:
        process_func = _compute_dt_correlations
        additional_args.update(correlation_kwargs)
    else:
        process_func = _compute_dt
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = [executor.submit(
            process_func, sub_catalog, master, **additional_args)
                   for sub_catalog, master in zip(sub_catalogs, catalog)]
    for master, result in zip(catalog, results):
        differential_times.update(
            {master.resource_id.id: result.result()})
    return differential_times, event_id_mapper


# dt.ct functions

def write_catalog(catalog, event_id_mapper=None, max_sep=8, min_link=8,
                  max_workers=None):
    """
    Generate a dt.ct file for hypoDD for a series of events.

    :type catalog: `obspy.core.event.Catalog`
    :param catalog: Catalog of events
    :type event_id_mapper: dict
    :param event_id_mapper:
        Dictionary mapping event resource id to an integer event id for hypoDD.
        If this is None, or missing events then the dictionary will be updated
        to include appropriate event-ids. This should be of the form
        {event.resource_id.id: integer_id}
    :type max_sep: float
    :param max_sep: Maximum separation between event pairs in km
    :type min_link: int
    :param min_link:
        Minimum links for an event to be paired, e.g. minimum number of picks
        from the same station and channel (and phase) that are shared between
        two events for them to be paired.
    :type max_workers: int
    :param max_workers:
        Maximum number of workers for parallel processing. If None, then all
        threads will be used.

    :returns: event_id_mapper
    """
    differential_times, event_id_mapper = compute_differential_times(
        catalog=catalog, correlation=False, event_id_mapper=event_id_mapper,
        max_sep=max_sep, min_link=min_link, max_workers=max_workers)
    with open("dt.ct", "w") as f:
        for master_id, linked_events in differential_times.items():
            for linked_event in linked_events:
                f.write(linked_event.ct_string)
                f.write("\n")
    return event_id_mapper


# dt.cc functions

def _filter_stream(event_id, st, lowcut, highcut):
    if lowcut is not None and highcut is not None:
        st_out = st.copy().detrend().filter(
            "bandpass", freqmin=lowcut, freqmax=highcut, corners=4,
            zerophase=True)
    elif lowcut is None and highcut is not None:
        st_out = st.copy().detrend().filter(
           "lowpass", freq=highcut, corners=4, zerophase=True)
    elif lowcut is not None and highcut is None:
        st_out = st.copy().detrend().filter(
           "highpass", freq=lowcut, corners=4, zerophase=True)
    else:
        st_out = st  # Don't need to copy if we aren't doing anything.
    return {event_id: st_out}


def write_correlations(catalog, stream_dict, extract_len, pre_pick,
                       shift_len, event_id_mapper=None, lowcut=1.0,
                       highcut=10.0, max_sep=8, min_link=8, cc_thresh=0.0,
                       interpolate=False, max_workers=None):
    """
    Write a dt.cc file for hypoDD input for a given list of events.

    Takes an input list of events and computes pick refinements by correlation.
    Outputs a single files, dt.cc file with weights as the square of the
    cross-correlation.

    :type catalog: `obspy.core.event.Catalog`
    :param catalog: Catalog of events to get differential times for
    :type event_id_mapper: dict
    :param event_id_mapper:
        Dictionary mapping event resource id to an integer event id for hypoDD.
        If this is None, or missing events then the dictionary will be updated
        to include appropriate event-ids. This should be of the form
        {event.resource_id.id: integer_id}
    :type stream_dict: dict
    :param stream_dict:
        Dictionary of streams keyed by event-id (the event.resource_id.id,
        NOT the hypoDD event-id)
    :type extract_len: float
    :param extract_len: Length in seconds to extract around the pick
    :type pre_pick: float
    :param pre_pick: Time before the pick to start the correlation window
    :type shift_len: float
    :param shift_len: Time to allow pick to vary in seconds
    :type lowcut: float
    :param lowcut: Lowcut in Hz - set to None to apply no lowcut
    :type highcut: float
    :param highcut: Highcut in Hz - set to None to apply no highcut
    :type max_sep: float
    :param max_sep: Maximum separation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired
    :type cc_thresh: float
    :param cc_thresh: Threshold to include cross-correlation results.
    :type interpolate: bool
    :param interpolate:
        Whether to interpolate correlations or not. Allows subsample accuracy
    :type max_workers: int
    :param max_workers:
        Maximum number of workers for parallel processing. If None then all
        threads will be used.

    :rtype: dict
    :returns: event_id_mapper

    .. note::
        You can provide processed waveforms, or let this function filter your
        data for you.  Filtering is undertaken by detrending and bandpassing
        with a 8th order zerophase butterworth filter.
    """
    max_workers = max_workers or cpu_count()
    # Process the streams
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = [executor.submit(
            _filter_stream, event_id, stream, lowcut, highcut)
            for event_id, stream in stream_dict.items()]
    processed_stream_dict = dict()
    for result in results:
        processed_stream_dict.update(result.result())

    correlation_times, event_id_mapper = compute_differential_times(
        catalog=catalog, correlation=True, event_id_mapper=event_id_mapper,
        max_sep=max_sep, min_link=min_link, max_workers=max_workers,
        stream_dict=processed_stream_dict, min_cc=cc_thresh,
        extract_len=extract_len, pre_pick=pre_pick, shift_len=shift_len,
        interpolate=interpolate)
    with open("dt.cc", "w") as f:
        for master_id, linked_events in correlation_times.items():
            for linked_event in linked_events:
                f.write(linked_event.cc_string)
                f.write("\n")
    return event_id_mapper


# Phase-file functions

def _hypodd_phase_pick_str(pick, event):
    """ Make a hypodd phase.dat style pick string. """
    arr = _get_arrival_for_pick(event=event, pick=pick)
    tt = pick.time - (event.preferred_origin() or event.origins[0]).time
    if arr:
        weight = arr.time_weight or 1.0
    else:
        weight = 1.0
    pick_str = "{station:5s} {tt:7.2f} {weight:5.3f} {phase:1s}".format(
        station=pick.waveform_id.station_code,
        tt=tt, weight=weight, phase=pick.phase_hint[0].upper())
    return pick_str


def _hypodd_phase_str(event, event_id_mapper):
    """ Form a phase-string for hypoDD from an obspy event. """
    try:
        origin = event.preferred_origin() or event.origins[0]
    except (IndexError, AttributeError):
        Logger.error("No origin for event {0}".format(event.resource_id.id))
        return
    try:
        magnitude = (event.preferred_magnitude() or event.magnitudes[0]).mag
    except (IndexError, AttributeError):
        Logger.warning("No magnitude")
        magnitude = 0.0
    try:
        time_error = origin.quality['standard_error']
    except AttributeError:
        Logger.warning('No time residual in header')
        time_error = 0.0

    z_err = (origin.depth_errors.uncertainty or 0.0) / 1000.
    # Note that these should be in degrees, but GeoNet uses meters.
    x_err = (origin.longitude_errors.uncertainty or 0.0) / 1000.
    y_err = (origin.latitude_errors.uncertainty or 0.0) / 1000.
    x_err = max(x_err, y_err)

    event_str = [(
        "# {year:4d} {month:02d} {day:02d} {hour:02d} {minute:02d} "
        "{second:02d} {latitude:10.6f} {longitude: 10.6f} {depth:10.6f} "
        "{magnitude:4.2f} {err_h:4.2f} {err_z:4.2f} {t_rms:4.2f} "
        "{event_id:9d}".format(
            year=origin.time.year, month=origin.time.month,
            day=origin.time.day, hour=origin.time.hour,
            minute=origin.time.minute, second=origin.time.second,
            latitude=origin.latitude, longitude=origin.longitude,
            depth=origin.depth / 1000., magnitude=magnitude,
            err_h=x_err, err_z=z_err, t_rms=time_error,
            event_id=event_id_mapper[event.resource_id.id]))]
    event_str.extend([_hypodd_phase_pick_str(pick, event)
                      for pick in event.picks if pick.phase_hint[0] in "PS"])
    return "\n".join(event_str)


def write_phase(catalog, event_id_mapper=None):
    """
    Write a phase.dat formatted file from an obspy catalog.

    :type catalog: `obspy.core.event.Catalog`
    :param catalog: Catalog of events
    :type event_id_mapper: dict
    :param event_id_mapper:
        Dictionary mapping event resource id to an integer event id for hypoDD.
        If this is None, or missing events then the dictionary will be updated
        to include appropriate event-ids. This should be of the form
        {event.resource_id.id: integer_id}

    :return: event_id_mapper
    """
    event_id_mapper = _generate_event_id_mapper(catalog, event_id_mapper)
    phase_strings = "\n".join(
        [_hypodd_phase_str(event, event_id_mapper) for event in catalog])
    with open("phase.dat", "w") as f:
        f.write(phase_strings)
    return event_id_mapper


def read_phase(ph_file):
    """
    Read hypoDD phase files into Obspy catalog class.

    :type ph_file: str
    :param ph_file: Phase file to read event info from.

    :returns: Catalog of events from file.
    :rtype: :class:`obspy.core.event.Catalog`

    >>> from obspy.core.event.catalog import Catalog
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> catalog = read_phase(TEST_PATH + '/tunnel.phase')
    >>> isinstance(catalog, Catalog)
    True
    """
    ph_catalog = Catalog()
    f = open(ph_file, 'r')
    # Topline of each event is marked by # in position 0
    for line in f:
        if line[0] == '#':
            if 'event_text' not in locals():
                event_text = {'header': line.rstrip(),
                              'picks': []}
            else:
                ph_catalog.append(_phase_to_event(event_text))
                event_text = {'header': line.rstrip(),
                              'picks': []}
        else:
            event_text['picks'].append(line.rstrip())
    ph_catalog.append(_phase_to_event(event_text))
    f.close()
    return ph_catalog


def _phase_to_event(event_text):
    """
    Function to convert the text for one event in hypoDD phase format to \
    event object.

    :type event_text: dict
    :param event_text: dict of two elements, header and picks, header is a \
        str, picks is a list of str.

    :returns: obspy.core.event.Event
    """
    ph_event = Event()
    # Extract info from header line
    # YR, MO, DY, HR, MN, SC, LAT, LON, DEP, MAG, EH, EZ, RMS, ID
    header = event_text['header'].split()
    ph_event.origins.append(Origin())
    ph_event.origins[0].time =\
        UTCDateTime(year=int(header[1]), month=int(header[2]),
                    day=int(header[3]), hour=int(header[4]),
                    minute=int(header[5]), second=int(header[6].split('.')[0]),
                    microsecond=int(float(('0.' + header[6].split('.')[1])) *
                                    1000000))
    ph_event.origins[0].latitude = float(header[7])
    ph_event.origins[0].longitude = float(header[8])
    ph_event.origins[0].depth = float(header[9]) * 1000
    ph_event.origins[0].quality = OriginQuality(
        standard_error=float(header[13]))
    ph_event.magnitudes.append(Magnitude())
    ph_event.magnitudes[0].mag = float(header[10])
    ph_event.magnitudes[0].magnitude_type = 'M'
    # Extract arrival info from picks!
    for i, pick_line in enumerate(event_text['picks']):
        pick = pick_line.split()
        _waveform_id = WaveformStreamID(station_code=pick[0])
        pick_time = ph_event.origins[0].time + float(pick[1])
        ph_event.picks.append(Pick(waveform_id=_waveform_id,
                                   phase_hint=pick[3],
                                   time=pick_time))
        ph_event.origins[0].arrivals.append(Arrival(phase=ph_event.picks[i],
                                                    pick_id=ph_event.picks[i].
                                                    resource_id))
        ph_event.origins[0].arrivals[i].time_weight = float(pick[2])
    return ph_event


# Event file functions

def _hypodd_event_str(event, event_id):
    """
    Make an event.dat style string for an event.

    :type event_id: int
    :param event_id: Must be an integer.
    """
    assert isinstance(event_id, int)
    try:
        origin = event.preferred_origin() or event.origins[0]
    except (IndexError, AttributeError):
        Logger.error("No origin for event {0}".format(event.resource_id.id))
        return
    try:
        magnitude = (event.preferred_magnitude() or event.magnitudes[0]).mag
    except (IndexError, AttributeError):
        Logger.warning("No magnitude")
        magnitude = 0.0
    try:
        time_error = origin.quality['standard_error']
    except AttributeError:
        Logger.warning('No time residual in header')
        time_error = 0.0

    z_err = (origin.depth_errors.uncertainty or 0.0) / 1000.
    # Note that these should be in degrees, but GeoNet uses meters.
    x_err = (origin.longitude_errors.uncertainty or 0.0) / 1000.
    y_err = (origin.latitude_errors.uncertainty or 0.0) / 1000.
    x_err = max(x_err, y_err)

    event_str = (
        "{year:4d}{month:02d}{day:02d}  {hour:2d}{minute:02d}"
        "{seconds:02d}{microseconds:02d}  {latitude:8.4f}  {longitude:9.4f}"
        "   {depth:8.4f}  {magnitude:5.2f}  {x_err:5.2f}  {z_err:5.2f}  "
        "{time_err:5.2f}  {event_id:9d}".format(
            year=origin.time.year, month=origin.time.month,
            day=origin.time.day, hour=origin.time.hour,
            minute=origin.time.minute, seconds=origin.time.second,
            microseconds=round(origin.time.microsecond / 1e4),
            latitude=origin.latitude, longitude=origin.longitude,
            depth=origin.depth / 1000., magnitude=magnitude,
            x_err=x_err, z_err=z_err, time_err=time_error, event_id=event_id))
    return event_str


def write_event(catalog, event_id_mapper=None):
    """
    Write obspy.core.event.Catalog to a hypoDD format event.dat file.

    :type catalog: obspy.core.event.Catalog
    :param catalog: A catalog of obspy events.
    :type event_id_mapper: dict
    :param event_id_mapper:
        Dictionary mapping event resource id to an integer event id for hypoDD.
        If this is None, or missing events then the dictionary will be updated
        to include appropriate event-ids. This should be of the form
        {event.resource_id.id: integer_id}

    :returns: dictionary of event-id mappings.
    """
    event_id_mapper = _generate_event_id_mapper(
        catalog=catalog, event_id_mapper=event_id_mapper)
    event_strings = [
        _hypodd_event_str(event, event_id_mapper[event.resource_id.id])
        for event in catalog]
    event_strings = "\n".join(event_strings)
    with open("event.dat", "w") as f:
        f.write(event_strings)
    return event_id_mapper


# Station.dat functions

def write_station(inventory):
    station_strings = []
    for network in inventory:
        for station in network:
            station_strings.append(
                "{station:<7s} {latitude:6.3f} {longitude:6.3f}".format(
                    station=station.code,
                    latitude=station.latitude,
                    longitude=station.longitude))
    with open("station.dat", "w") as f:
        f.write("\n".join(station_strings))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
