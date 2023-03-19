"""
Functions to generate hypoDD input files from catalogs.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import logging
from collections import namedtuple, defaultdict, Counter
from multiprocessing import cpu_count, Pool, shared_memory

from obspy import UTCDateTime, Stream
from obspy.core.event import (
    Catalog, Event, Origin, Magnitude, Pick, WaveformStreamID, Arrival,
    OriginQuality)
from eqcorrscan.utils.clustering import dist_mat_km
from eqcorrscan.core.lag_calc import _concatenate_and_correlate, _xcorr_interp
from eqcorrscan.utils.correlate import pool_boy

Logger = logging.getLogger(__name__)

SeedPickID = namedtuple("SeedPickID", ["seed_id", "phase_hint"])


# Some hypoDD specific event holders - classes were faster than named-tuples

class SparseEvent(object):
    def __init__(self, resource_id, picks, origin_time):
        self.resource_id = resource_id
        self.picks = picks
        self.origin_time = origin_time

    def __repr__(self):
        return ("SparseEvent(resource_id={0}, origin_time={1}, picks=[{2} "
                "picks])".format(
                    self.resource_id, self.origin_time, len(self.picks)))


class SparsePick(object):
    def __init__(self, tt, time, time_weight, seed_id, phase_hint,
                 waveform_id):
        self.tt = tt
        self.time = time
        self.time_weight = time_weight
        self.seed_id = seed_id
        self.phase_hint = phase_hint
        self.waveform_id = waveform_id

    def __repr__(self):
        return ("SparsePick(seed_id={0}, phase_hint={1}, tt={2:.2f}, "
                "time_weight={3})".format(
                    self.seed_id, self.phase_hint, self.tt, self.time_weight))

    @property
    def station(self):
        return self.seed_id.split('.')[1]

    @property
    def channel(self):
        return self.seed_id.split('.')[-1]


# Generic helpers

class _DTObs(object):
    """ Holder for phase observations """

    def __init__(self, station, tt1, tt2, weight, phase):
        self.station = station
        assert len(self.station) <= 7, "Station must be <= five characters"
        self.tt1 = tt1
        self.tt2 = tt2
        self.weight = weight
        self.phase = phase
        assert self.phase in "PS", "Only P or S phases accepted"

    def __repr__(self):
        return ("_DTObs(station={0}, tt1={1:.2f}, tt2={2:.2f}, weight={3:.2f},"
                " phase={4})".format(self.station, self.tt1, self.tt2,
                                     self.weight, self.phase))

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

    def __repr__(self):
        return ("_EventPair(event_id_1={0}, event_id_2={1}, obs=[{2} "
                "observations])".format(self.event_id_1, self.event_id_2,
                                        len(self.obs)))

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
        if str(event.resource_id) not in event_id_mapper.keys():
            event_id = largest_event_id + 1
            largest_event_id = event_id
            event_id_mapper.update({str(event.resource_id): event_id})
    return event_id_mapper


def _make_sparse_event(event):
    """ Make a sparse event with just the info hypoDD needs. """
    origin_time = (event.preferred_origin() or event.origins[0]).time
    time_weight_dict = {
        arr.pick_id: arr.time_weight or 1.0 for origin in event.origins
        for arr in origin.arrivals}
    sparse_event = SparseEvent(
        resource_id=event.resource_id.id,
        origin_time=origin_time,
        picks=[SparsePick(
            tt=pick.time - origin_time,
            time=pick.time,
            seed_id=pick.waveform_id.get_seed_string(),
            phase_hint=pick.phase_hint[0],  # Only use P or S hints.
            time_weight=time_weight_dict.get(pick.resource_id, 1.0),
            waveform_id=pick.waveform_id)
            for pick in event.picks])
    return sparse_event


def _prepare_stream(stream, event, extract_len, pre_pick, seed_pick_ids=None):
    """
    Slice stream around picks

    returns a dictionary of traces keyed by phase_hint.
    """
    seed_pick_ids = seed_pick_ids or {
        SeedPickID(pick.waveform_id.get_seed_string(), pick.phase_hint[0])
        for pick in event.picks if pick.phase_hint.startswith(("P", "S"))}
    stream_sliced = defaultdict(Stream)
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
        tr = stream.select(id=seed_pick_id.seed_id).merge()
        if len(tr) == 0:
            continue
        else:
            tr = tr[0]
        Logger.debug(
            f"Trimming trace on {tr.id} between {tr.stats.starttime} - "
            f"{tr.stats.endtime} to {pick.time - pre_pick} - "
            f"{(pick.time - pre_pick) + extract_len}")
        tr = stream.select(id=seed_pick_id.seed_id).slice(
            starttime=pick.time - pre_pick,
            endtime=(pick.time - pre_pick) + extract_len).merge()
        if len(tr) == 0:
            continue
        if len(tr) > 1:
            Logger.error("Multiple traces for {seed_id}".format(
                seed_id=seed_pick_id.seed_id))
            continue
        tr = tr[0]

        # If there is one sample too many after this remove the first one
        # by convention
        n_samples_intended = extract_len * tr.stats.sampling_rate
        if len(tr.data) == n_samples_intended + 1:
            tr.data = tr.data[1:len(tr.data)]
        # if tr.stats.endtime - tr.stats.starttime != extract_len:
        if tr.stats.npts < n_samples_intended:
            Logger.warning(
                "Insufficient data ({rlen} s) for {tr_id}, discarding. Check "
                "that your traces are at least of length {length} s, with a "
                "pre_pick time of at least {prepick} s!".format(
                    rlen=tr.stats.endtime - tr.stats.starttime,
                    tr_id=tr.id, length=extract_len, prepick=pre_pick))
            continue
        stream_sliced.update(
            {seed_pick_id.phase_hint:
             stream_sliced[seed_pick_id.phase_hint] + tr})
    return stream_sliced


# Time calculators

def _compute_dt_correlations(catalog, master, min_link, event_id_mapper,
                             stream_dict, min_cc, extract_len, pre_pick,
                             shift_len, interpolate, max_workers=1,
                             shm_data_shape=None, shm_dtype=None,
                             weight_by_square=True, **kwargs):
    """ Compute cross-correlation delay times. """
    max_workers = max_workers or 1
    Logger.info(
        f"Correlating {str(master.resource_id)} with {len(catalog)} events")
    differential_times_dict = dict()
    # Assign trace data from shared memory
    for (key, stream) in stream_dict.items():
        for tr in stream:
            if len(tr.data) == 0 and hasattr(tr, 'shared_memory_name'):
                shm = shared_memory.SharedMemory(name=tr.shared_memory_name)
                # Reconstructing numpy data array
                sm_data = np.ndarray(
                    shm_data_shape, dtype=shm_dtype, buffer=shm.buf)
                tr.data = np.zeros_like(sm_data)
                # Copy data into process memory
                tr.data[:] = sm_data[:]

    master_stream = _prepare_stream(
        stream=stream_dict[str(master.resource_id)], event=master,
        extract_len=extract_len, pre_pick=pre_pick)
    available_seed_ids = {tr.id for st in master_stream.values() for tr in st}
    Logger.debug(f"The channels provided are: {available_seed_ids}")
    master_seed_ids = {
        SeedPickID(pick.waveform_id.get_seed_string(), pick.phase_hint[0])
        for pick in master.picks if
        pick.phase_hint[0] in "PS" and
        pick.waveform_id.get_seed_string() in available_seed_ids}
    Logger.debug(f"Using channels: {master_seed_ids}")
    # Dictionary of travel-times for master keyed by {station}_{phase_hint}
    master_tts = dict()
    try:
        master_origin_time = (
            master.preferred_origin() or master.origins[0]).time
    except AttributeError:  # In case it's a SparseEvent
        master_origin_time = master.origin_time
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
    event_dict = {str(event.resource_id): event for event in catalog}
    event_ids = set(event_dict.keys())
    # Check for overlap
    _stream_event_ids = set(stream_dict.keys())
    if len(event_ids.difference(_stream_event_ids)):
        Logger.warning(
            f"Missing streams for {event_ids.difference(_stream_event_ids)}")
        # Just use the event ids that we actually have streams for!
        event_ids = event_ids.intersection(_stream_event_ids)
    # Reorder event_ids according to original order
    event_ids = [key for key in event_dict.keys() if key in event_ids]

    if max_workers > 1:
        with pool_boy(Pool, len(event_ids), cores=max_workers) as pool:
            results = [pool.apply_async(
                _prepare_stream,
                args=(stream_dict[event_id], event_dict[event_id],
                      matched_length, matched_pre_pick),
                kwds=dict(seed_pick_ids=master_seed_ids))
                        for event_id in event_ids]
        matched_streams = {id_res[0]: id_res[1].get()
                           for id_res in zip(event_ids, results)}
    else:
        matched_streams = {
            event_id: _prepare_stream(
                stream=stream_dict[event_id], event=event_dict[event_id],
                extract_len=matched_length, pre_pick=matched_pre_pick,
                seed_pick_ids=master_seed_ids)
            for event_id in event_ids}

    sampling_rates = {tr.stats.sampling_rate for st in master_stream.values()
                      for tr in st}
    for phase_hint in master_stream.keys():  # Loop over P and S separately
        for sampling_rate in sampling_rates:  # Loop over separate samp rates
            delta = 1.0 / sampling_rate
            _master_stream = master_stream[phase_hint].select(
                sampling_rate=sampling_rate)
            if len(_master_stream) == 0:
                continue
            _matched_streams = dict()
            for key, value in matched_streams.items():
                _st = value[phase_hint].select(sampling_rate=sampling_rate)
                if len(_st) > 0:
                    _matched_streams.update({key: _st})
            if len(_matched_streams) == 0:
                Logger.info("No matching data for {0}, {1} phase".format(
                    str(master.resource_id), phase_hint))
                continue
            # Check lengths
            master_length = [tr.stats.npts for tr in _master_stream]
            if len(set(master_length)) > 1:
                Logger.warning("Multiple lengths found - check that you "
                               "are providing sufficient data")
            master_length = Counter(master_length).most_common(1)[0][0]
            _master_stream = _master_stream.select(npts=master_length)
            matched_length = Counter(
                (tr.stats.npts for st in _matched_streams.values()
                 for tr in st))
            if len(matched_length) > 1:
                Logger.warning("Multiple lengths of stream found - taking "
                               "the most common. Check that you are "
                               "providing sufficient data")
            matched_length = matched_length.most_common(1)[0][0]
            if matched_length < master_length:
                Logger.error("Matched streams are shorter than the master, "
                             "will not correlate")
                continue
            # Remove empty streams and generate an ordered list of event_ids
            used_event_ids, used_matched_streams = [], []
            for event_id, _matched_stream in _matched_streams.items():
                _matched_stream = _matched_stream.select(npts=matched_length)
                if len(_matched_stream) > 0:
                    used_event_ids.append(event_id)
                    used_matched_streams.append(_matched_stream)
            # Check that there are matching seed ids.
            master_seed_ids = set(tr.id for tr in _master_stream)
            matched_seed_ids = set(
                tr.id for st in used_matched_streams for tr in st)
            if not matched_seed_ids.issubset(master_seed_ids):
                Logger.warning(
                    "After checking length there are no matched traces: "
                    f"master: {master_seed_ids}, matched: {matched_seed_ids}")
                continue
            # Do the correlations
            Logger.debug(
                f"Correlating channels: {[tr.id for tr in _master_stream]}")
            ccc_out, used_chans = _concatenate_and_correlate(
                template=_master_stream, streams=used_matched_streams,
                cores=max_workers)
            # Convert ccc_out to pick-time
            for i, used_event_id in enumerate(used_event_ids):
                for j, chan in enumerate(used_chans[i]):
                    if not chan.used:
                        continue
                    correlation = ccc_out[i][j]
                    if interpolate:
                        shift, cc_max = _xcorr_interp(correlation, dt=delta,
                                                      **kwargs)
                    else:
                        cc_max = np.amax(correlation)
                        shift = np.argmax(correlation) * delta
                    if cc_max < min_cc:
                        continue
                    shift -= shift_len
                    pick = [p for p in event_dict[used_event_id].picks
                            if p.phase_hint[0] == phase_hint
                            and p.waveform_id.station_code == chan.channel[0]
                            and p.waveform_id.channel_code == chan.channel[1]]
                    pick = sorted(pick, key=lambda p: p.time)[0]
                    try:
                        tt2 = pick.time - (
                                event_dict[used_event_id].preferred_origin() or
                                event_dict[used_event_id].origins[0]).time
                    except AttributeError:
                        tt2 = pick.time - event_dict[used_event_id].origin_time
                    tt2 += shift
                    diff_time = differential_times_dict.get(
                        used_event_id, None)
                    if diff_time is None:
                        diff_time = _EventPair(
                            event_id_1=event_id_mapper[
                                str(master.resource_id)],
                            event_id_2=event_id_mapper[used_event_id])
                    weight = cc_max
                    if weight_by_square:
                        weight **= 2
                    diff_time.obs.append(
                        _DTObs(station=chan.channel[0],
                               tt1=master_tts["{0}_{1}".format(
                                   chan.channel[0], phase_hint)],
                               tt2=tt2, weight=weight,
                               phase=phase_hint[0]))
                    differential_times_dict.update({used_event_id: diff_time})
    # Threshold on min_link
    differential_times = [dt for dt in differential_times_dict.values()
                          if len(dt.obs) >= min_link]
    return differential_times


def _compute_dt(sparse_catalog, master, min_link, event_id_mapper):
    """
    Inner function to compute differential times between a catalog and a
    master event.
    """
    return [_make_event_pair(
        sparse_event=event, master=master, event_id_mapper=event_id_mapper,
        min_link=min_link) for event in sparse_catalog]


def _make_event_pair(sparse_event, master, event_id_mapper, min_link):
    """
    Make an event pair for a given event and master event.
    """
    differential_times = _EventPair(
        event_id_1=event_id_mapper[master.resource_id],
        event_id_2=event_id_mapper[sparse_event.resource_id])
    for master_pick in master.picks:
        if master_pick.phase_hint and \
                master_pick.phase_hint not in "PS":  # pragma: no cover
            continue
        matched_picks = [p for p in sparse_event.picks
                         if p.station == master_pick.station
                         and p.phase_hint == master_pick.phase_hint]
        for matched_pick in matched_picks:
            differential_times.obs.append(
                _DTObs(station=master_pick.station,
                       tt1=master_pick.tt, tt2=matched_pick.tt,
                       weight=(master_pick.time_weight +
                               matched_pick.time_weight) / 2.0,
                       phase=master_pick.phase_hint))
    if len(differential_times.obs) >= min_link:
        return differential_times
    return


def _prep_horiz_picks(catalog, stream_dict, event_id_mapper):
    """
    Fill in horizontal picks for the alternate horizontal channel for events in
    catalog.
    """
    # keep user input safe
    catalog = catalog.copy()
    for event in catalog:
        event_S_picks = [
            pick for pick in event.picks if pick.phase_hint.upper().startswith(
                'S') and pick.waveform_id.get_seed_string()[-1] in 'EN12XY']
        st = stream_dict[str(event.resource_id)]
        st = Stream([tr for tr in st if tr.stats.channel[-1] in 'EN12XY'])
        for tr in st:
            tr_picks = [
                pick for pick in event_S_picks
                if tr.id == pick.waveform_id.get_seed_string()]
            if len(tr_picks) > 0:
                continue
            else:
                tr_picks = [
                    pick for pick in event_S_picks
                    if tr.id[0:-1] == pick.waveform_id.get_seed_string()[0:-1]]
                new_wav_id = WaveformStreamID(network_code=tr.stats.network,
                                              station_code=tr.stats.station,
                                              location_code=tr.stats.location,
                                              channel_code=tr.stats.channel)
                for pick in tr_picks:
                    new_pick = SparsePick(tt=pick.tt, time=pick.time,
                                          time_weight=pick.time_weight,
                                          seed_id=new_wav_id.get_seed_string(),
                                          phase_hint=pick.phase_hint,
                                          waveform_id=new_wav_id)
                    event.picks.append(new_pick)
    return catalog


def stream_dict_to_shared_mem(stream_dict):
    """
    Move the data of streams from a dict of (key, obspy.stream) into shared
    memory so that the data can be retrieved by multiple processes in parallel.
    This can help speed up parallel execution because the initiation of each
    worker process becomes cheaper (less data to transfer). For now this only
    puts the numpy array in trace.data into shared memory (because it's easy).

    :type stream_dict: dict of (key, `obspy.stream`)
    :param stream_dict: dict of streams that should be moved to shared memory

    :returns: stream_dict, shm_name_list, shm_data_shapes, shm_data_dtypes

    :rtype: dict
    :return: Dictionary streams that were moved to shared memory
    :rtype: list
    :return: List of names to the shared memory address for each trace.
    :rtype: list
    :return:
        List of numpy-array shaped for each trace-data array in shared memory.
    :rtype: list
    :return: List of data types for each trace-data-array in shared memory.

    """
    shm_name_list = []
    shm_data_shapes = []
    shm_data_dtypes = []
    shm_references = []
    for (key, stream) in stream_dict.items():
        for tr in stream:
            data_array = tr.data
            # Let SharedMemory create suitable filename itself:
            shm = shared_memory.SharedMemory(
                create=True, size=data_array.nbytes)
            shm_name_list.append(shm.name)
            shm_references.append(shm)
            # Now create a NumPy array backed by shared memory
            shm_data_shape = data_array.shape
            shm_data_dtype = data_array.dtype
            shared_data_array = np.ndarray(
                shm_data_shape, dtype=shm_data_dtype, buffer=shm.buf)
            # Copy the original data into shared memory
            shared_data_array[:] = data_array[:]
            # tr.data = shared_data_array
            tr.data = np.array([])
            tr.shared_memory_name = shm.name
            shm_data_shapes.append(shm_data_shape)
            shm_data_dtypes.append(shm_data_dtype)
    shm_data_shapes = list(set(shm_data_shapes))
    shm_data_dtypes = list(set(shm_data_dtypes))
    return (stream_dict, shm_name_list, shm_references, shm_data_shapes,
            shm_data_dtypes)


def compute_differential_times(catalog, correlation, stream_dict=None,
                               event_id_mapper=None, max_sep=8., min_link=8,
                               min_cc=None, extract_len=None, pre_pick=None,
                               shift_len=None, interpolate=False,
                               all_horiz=False, max_workers=None,
                               max_trace_workers=1, use_shared_memory=False,
                               weight_by_square=True, *args, **kwargs):
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
    :param shift_len:
        Time (+/-) to allow pick to vary in seconds. e.g. if shift_len
        is set to 1s, the pick will be allowed to shift between
        pick_time - 1 and pick_time + 1.
    :type interpolate: bool
    :param interpolate:
        Whether to interpolate correlations or not. Allows subsample accuracy
    :type max_workers: int
    :param max_workers:
        Maximum number of workers for parallel correlation of events. If None
        then all threads will be used.
    :type max_trace_workers: int
    :param max_trace_workers:
        Maximum number of workers for parallel correlation of traces insted of
        events. If None then all threads will be used (but can only be used
        when max_workers = 1).
    :type use_shared_memory: bool
    :param use_shared_memory:
        Whether to move trace data arrays into shared memory for computing
        trace correlations. Can speed up total execution time by ~20 % for
        hypodd-correlations with a lot of clustered seismicity.
    :type weight_by_square: bool
    :param weight_by_square:
        Whether to compute correlation weights as the square of the maximum
        correlation (True), or the maximum correlation (False).

    :rtype: dict
    :return: Dictionary of differential times keyed by event id.
    :rtype: dict
    :return: Dictionary of event_id_mapper

    .. note::
        The arguments min_cc, stream_dict, extract_len, pre_pick, shift_len
        and interpolate are only required if correlation=True.

        Two parallelization strategies are available for correlating waveforms:
        parallelization across events (default) or across each event's traces
        (when max_workers = 1 and max_traces_workers > 1). The former is often
        quicker for short traces because it generally loads the CPU better for
        multiple events and may require more memory, but the latter can be
        quicker for few events with many or very long traces and requires less
        memory.

    .. note::
        Differential times are computed as travel-time for event 1 minus
        travel-time for event 2 (tt1 - tt2).
    """
    include_master = kwargs.get("include_master", False)
    correlation_kwargs = dict(
        min_cc=min_cc, stream_dict=stream_dict, extract_len=extract_len,
        pre_pick=pre_pick, shift_len=shift_len, interpolate=interpolate,
        max_workers=max_workers, weight_by_square=weight_by_square)
    for key, value in kwargs.items():
        correlation_kwargs.update({key: value})
    if correlation:
        for arg, name in correlation_kwargs.items():
            assert arg is not None, "{0} is required for correlation".format(
                name)
    # Ensure all events have locations and picks.
    event_id_mapper = _generate_event_id_mapper(
        catalog=catalog, event_id_mapper=event_id_mapper)
    distances = dist_mat_km(catalog)
    distance_filter = distances <= max_sep
    if not include_master:
        np.fill_diagonal(distance_filter, 0)
        # Do not match events to themselves - this is the default,
        # only included for testing
    # Reformat catalog to sparse catalog
    sparse_catalog = [_make_sparse_event(ev) for ev in catalog]
    if all_horiz:
        sparse_catalog = _prep_horiz_picks(sparse_catalog, stream_dict,
                                           event_id_mapper)

    additional_args = dict(min_link=min_link, event_id_mapper=event_id_mapper)
    if correlation:
        differential_times = {}
        additional_args.update(correlation_kwargs)
        n = len(sparse_catalog)
        if max_workers == 1:
            # If desired, parallelize over traces instead of events:
            max_trace_workers = max_trace_workers or cpu_count()
            additional_args.update(dict(max_workers=max_trace_workers))
            for i, master in enumerate(sparse_catalog):
                master_id = str(master.resource_id)
                sub_catalog = [ev for j, ev in enumerate(sparse_catalog)
                               if distance_filter[i][j]]
                if master_id not in additional_args["stream_dict"].keys():
                    Logger.warning(
                        f"{master_id} not in waveforms, skipping")
                    continue
                differential_times.update({
                    master_id: _compute_dt_correlations(
                        sub_catalog, master, **additional_args)})
                Logger.info(
                    f"Completed correlations for core event {i} of {n}")
        else:
            sub_catalogs = ([ev for i, ev in enumerate(sparse_catalog)
                             if master_filter[i]]
                            for master_filter in distance_filter)
            # Move trace data into shared memory
            if use_shared_memory:
                (shm_stream_dict, shm_name_list, shm_references,
                 shm_data_shapes, shm_dtypes) = (
                    stream_dict_to_shared_mem(stream_dict))
                if len(shm_data_shapes) == 1 and len(shm_dtypes) == 1:
                    shm_data_shape = shm_data_shapes[0]
                    shm_dtype = shm_dtypes[0]
                    additional_args.update({'stream_dict': shm_stream_dict})
                    additional_args.update({'shm_data_shape': shm_data_shape})
                    additional_args.update({'shm_dtype': shm_dtype})
                else:
                    use_shared_memory = False
            with pool_boy(Pool, n, cores=max_workers) as pool:
                # Parallelize over events instead of traces
                additional_args.update(dict(max_workers=1))
                results = [
                    pool.apply_async(
                        _compute_dt_correlations,
                        args=(sub_catalog, master), kwds=additional_args)
                    for sub_catalog, master in zip(sub_catalogs,
                                                   sparse_catalog)
                    if str(master.resource_id) in additional_args[
                        "stream_dict"].keys()]
                Logger.info('Submitted asynchronous jobs to workers.')
                differential_times = {
                    master.resource_id: result.get()
                    for master, result in zip(sparse_catalog, results)
                    if str(master.resource_id) in additional_args[
                        "stream_dict"].keys()}
                Logger.debug('Got results from workers.')
                # Destroy shared memory
                if use_shared_memory:
                    for shm_name in shm_name_list:
                        shm = shared_memory.SharedMemory(name=shm_name)
                        shm.close()
                        shm.unlink()
    else:
        sub_catalogs = ([ev for i, ev in enumerate(sparse_catalog)
                         if master_filter[i]]
                        for master_filter in distance_filter)
        max_workers = max_workers or cpu_count()
        if max_workers > 1:
            with pool_boy(
                    Pool, len(sparse_catalog), cores=max_workers) as pool:
                results = [pool.apply_async(
                    _compute_dt,
                    args=(sub_catalog, master), kwds=additional_args)
                           for master, sub_catalog in zip(
                               sparse_catalog, sub_catalogs)]
                differential_times = {
                    master.resource_id: result.get()
                    for master, result in zip(sparse_catalog, results)}
        else:
            differential_times = {
                master.resource_id: _compute_dt(
                    sub_catalog, master, **additional_args)
                for master, sub_catalog in zip(sparse_catalog, sub_catalogs)}

    # Remove Nones
    for key, value in differential_times.items():
        differential_times.update({key: [v for v in value if v is not None]})
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
        Maximum number of workers for parallel processing. If None then all
        threads will be used.

    :returns: event_id_mapper

    .. note::
        Differential times are computed as travel-time for event 1 minus
        travel-time for event 2 (tt1 - tt2).
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

def _meta_filter_stream(event_id, stream_dict, lowcut, highcut):
    return _filter_stream(
        event_id=event_id, st=stream_dict[event_id], lowcut=lowcut,
        highcut=highcut)


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
                       highcut=10.0, max_sep=8, min_link=8,  min_cc=0.0,
                       interpolate=False, all_horiz=False, max_workers=None,
                       parallel_process=False, weight_by_square=True,
                       *args, **kwargs):
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
    :type min_cc: float
    :param min_cc: Threshold to include cross-correlation results.
    :type interpolate: bool
    :param interpolate:
        Whether to interpolate correlations or not. Allows subsample accuracy
    :type max_workers: int
    :param max_workers:
        Maximum number of workers for parallel processing. If None then all
        threads will be used.
    :type parallel_process: bool
    :param parallel_process:
        Whether to process streams in parallel or not. Experimental, may use
        too much memory.
    :type weight_by_square: bool
    :param weight_by_square:
        Whether to compute correlation weights as the square of the maximum
        correlation (True), or the maximum correlation (False).

    :rtype: dict
    :returns: event_id_mapper

    .. note::
        You can provide processed waveforms, or let this function filter your
        data for you.  Filtering is undertaken by detrending and bandpassing
        with a 8th order zerophase butterworth filter.

    .. note::
        Differential times are computed as travel-time for event 1 minus
        travel-time for event 2 (tt1 - tt2).
    """
    # Depreciated argument
    cc_thresh = kwargs.get("cc_thresh", None)
    if cc_thresh:
        min_cc = cc_thresh
        Logger.warning("cc_thresh is depreciated, use min_cc instead")
    max_workers = max_workers or cpu_count()
    processed_stream_dict = stream_dict
    # Process the streams
    if not (lowcut is None and highcut is None):
        processed_stream_dict = dict()
        if parallel_process:
            max_process_workers = int(max(np.array(
                [max_workers, kwargs.get('max_trace_workers')],
                dtype=np.float64)))
            with pool_boy(
                    Pool, len(stream_dict), cores=max_process_workers) as pool:
                results = [pool.apply_async(
                    _meta_filter_stream,
                    (key, stream_dict, lowcut, highcut))
                           for key in stream_dict.keys()]
            for result in results:
                processed_stream_dict.update(result.get())
        else:
            for key in stream_dict.keys():
                processed_stream_dict.update(_meta_filter_stream(
                    stream_dict=stream_dict, lowcut=lowcut, highcut=highcut,
                    event_id=key))
    correlation_times, event_id_mapper = compute_differential_times(
        catalog=catalog, correlation=True, event_id_mapper=event_id_mapper,
        max_sep=max_sep, min_link=min_link, max_workers=max_workers,
        stream_dict=processed_stream_dict, min_cc=min_cc,
        extract_len=extract_len, pre_pick=pre_pick, shift_len=shift_len,
        interpolate=interpolate, all_horiz=all_horiz,
        weight_by_square=weight_by_square, **kwargs)
    with open("dt.cc", "w") as f:
        for master_id, linked_events in correlation_times.items():
            for linked_event in linked_events:
                f.write(linked_event.cc_string)
                f.write("\n")
    return event_id_mapper


# Phase-file functions

def _hypodd_phase_pick_str(pick, sparse_event):
    """ Make a hypodd phase.dat style pick string. """
    pick_str = "{station:5s} {tt:7.4f} {weight:5.3f} {phase:1s}".format(
        station=pick.waveform_id.station_code,
        tt=pick.tt, weight=pick.weight, phase_hint=pick.phase_hint[0].upper())
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
        time_error = origin.quality['standard_error'] or 0.0
    except (TypeError, AttributeError):
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
    sparse_event = _make_sparse_event(event)
    for pick in sparse_event.picks:
        if pick.phase_hint[0] not in "PS":
            continue
        event_str.append(
            "{station:5s} {tt:7.2f} {weight:5.3f} {phase_hint:1s}".format(
                station=pick.station, tt=pick.tt, weight=pick.time_weight,
                phase_hint=pick.phase_hint[0].upper()))
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
    ph_event.origins[0].time = \
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
        time_error = origin.quality['standard_error'] or 0.0
    except (TypeError, AttributeError):
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

def write_station(inventory, use_elevation=False, filename="station.dat"):
    """
    Write a hypoDD formatted station file.

    :type inventory: obspy.core.Inventory
    :param inventory:
        Inventory of stations to write - should include channels if
        use_elevation=True to incorporate channel depths.
    :type use_elevation: bool
    :param use_elevation: Whether to write elevations (requires hypoDD >= 2)
    :type filename: str
    :param filename: File to write stations to.
    """
    station_strings = []
    formatter = "{sta:<7s} {lat:>9.5f} {lon:>10.5f}"
    if use_elevation:
        formatter = " ".join([formatter, "{elev:>5.0f}"])

    for network in inventory:
        for station in network:
            parts = dict(sta=station.code, lat=station.latitude,
                         lon=station.longitude)
            if use_elevation:
                channel_depths = {chan.depth for chan in station}
                if len(channel_depths) == 0:
                    Logger.warning("No channels provided, using 0 depth.")
                    depth = 0.0
                else:
                    depth = channel_depths.pop()
                if len(channel_depths) > 1:
                    Logger.warning(
                        f"Multiple depths for {station.code}, using {depth}")
                parts.update(dict(elev=station.elevation - depth))
            station_strings.append(formatter.format(**parts))
    with open(filename, "w") as f:
        f.write("\n".join(station_strings))


if __name__ == '__main__':
    import doctest

    doctest.testmod()
