"""
Functions to generate hypoDD input files from catalogs.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import os
import glob
import logging
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

from obspy import read, UTCDateTime
from obspy.core.event import (
    Catalog, Event, Origin, Magnitude, Pick, WaveformStreamID, Arrival,
    OriginQuality)
from obspy.signal.cross_correlation import xcorr_pick_correction
from eqcorrscan.utils.mag_calc import dist_calc
from eqcorrscan.utils.clustering import dist_mat_km


Logger = logging.getLogger(__name__)


class _DTObs(object):
    """ Holder for phase observations """
    def __init__(self, station, tt1, tt2, weight, phase):
        self.station = station
        self.tt1 = tt1
        self.tt2 = tt2
        self.weight = weight
        self.phase = phase

    def __repr__(self):
        return "{sta:<5s} {tt1:7.3f} {tt2:7.3f} {weight:6.4f} {phase}".format(
            sta=self.station, tt1=self.tt1, tt2=self.tt2, weight=self.weight,
            phase=self.phase)


class _EventPair(object):
    """ Holder for event paid observations. """
    def __init__(self, event_id_1, event_id_2, obs=None):
        self.event_id_1 = event_id_1
        self.event_id_2 = event_id_2
        self.obs = obs or list()

    def __repr__(self):
        header = "# {event_id_1:>9d} {event_id_2:>9d}".format(
            event_id_1=self.event_id_1, event_id_2=self.event_id_2)
        parts = [header]
        parts.extend([str(obs) for obs in self.obs])
        return '\n'.join(parts)


def _cc_round(num, dp):
    """
    Convenience function to take a float and round it to dp padding with zeros
    to return a string

    :type num: float
    :param num: Number to round
    :type dp: int
    :param dp: Number of decimal places to round to.

    :returns: str

    >>> print(_cc_round(0.25364, 2))
    0.25
    """
    num = round(num, dp)
    num = '{0:.{1}f}'.format(num, dp)
    return num


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


def _get_arrival_for_pick(event, pick):
    """
    Get a matching arrival for a given pick.
    """
    matched_arrivals = [
        arr for origin in event.origins for arr in origin.arrivals
        if arr.pick_id.get_referred_object() == pick]
    if len(matched_arrivals) == 0:
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


def compute_differential_times(catalog, event_id_mapper=None, max_sep=8.,
                               min_link=8, max_workers=None):
    """
    Generate groups of differential times for a catalog.

    :type catalog: `obspy.core.event.Catalog`
    :param catalog: Catalog of events to get differential times for
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
    :type max_workers: int
    :param max_workers:
        Maximum number of workers for parallel processing. If None then all
        threads will be used.

    :rtype: dict
    :return: Dictionary of differential times keyed by event id.
    :rtype: dict
    :return: Dictionary of event_id_mapper
    """
    max_workers = max_workers or cpu_count()
    # Ensure all events have locations and picks.
    event_id_mapper = _generate_event_id_mapper(
        catalog=catalog, event_id_mapper=event_id_mapper)
    distances = dist_mat_km(catalog)
    distance_filter = distances <= max_sep

    differential_times = {}
    sub_catalogs = [[ev for i, ev in enumerate(catalog)
                     if master_filter[i]] for master_filter in distance_filter]
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = [executor.submit(
            _compute_dt, sub_catalog, master, min_link, event_id_mapper)
                   for sub_catalog, master in zip(sub_catalogs, catalog)]
    for master, result in zip(catalog, results):
        differential_times.update(
            {master.resource_id.id: result.result()})
    return differential_times, event_id_mapper


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
                      for pick in event.picks])
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
        catalog=catalog, event_id_mapper=event_id_mapper, max_sep=max_sep,
        min_link=min_link, max_workers=max_workers)
    with open("dt.ct", "w") as f:
        for master_id, linked_events in differential_times.items():
            for linked_event in linked_events:
                f.write(str(linked_event))
                f.write("\n")
    return event_id_mapper

# TODO: Use new correlation functions - wrap lag-calc
def write_correlations(catalog, event_id_mapper, wavbase, extract_len,
                       pre_pick, shift_len, lowcut=1.0, highcut=10.0,
                       max_sep=8, min_link=8, cc_thresh=0.0, plotvar=False):
    """
    Write a dt.cc file for hypoDD input for a given list of events.

    Takes an input list of events and computes pick refinements by correlation.
    Outputs two files, dt.cc and dt.cc2, each provides a different weight,
    dt.cc uses weights of the cross-correlation, and dt.cc2 provides weights
    as the square of the cross-correlation.

    :type wavbase: str
    :param wavbase: Path to the seisan wave directory that the wavefiles in the
                    S-files are stored
    :type extract_len: float
    :param extract_len: Length in seconds to extract around the pick
    :type pre_pick: float
    :param pre_pick: Time before the pick to start the correlation window
    :type shift_len: float
    :param shift_len: Time to allow pick to vary
    :type lowcut: float
    :param lowcut: Lowcut in Hz - default=1.0
    :type highcut: float
    :param highcut: Highcut in Hz - default=10.0
    :type max_sep: float
    :param max_sep: Maximum separation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired
    :type cc_thresh: float
    :param cc_thresh: Threshold to include cross-correlation results.
    :type plotvar: bool
    :param plotvar: To show the pick-correction plots, defualts to False.

    .. warning:: This is not a fast routine!

    .. warning::
        In contrast to seisan's corr routine, but in accordance with the
        hypoDD manual, this outputs corrected differential time.

    .. note::
        Currently we have not implemented a method for taking
        unassociated event objects and wavefiles.  As such if you have events \
        with associated wavefiles you are advised to generate Sfiles for each \
        event using the sfile_util module prior to this step.

    .. note::
        There is no provision to taper waveforms within these functions, if you
        desire this functionality, you should apply the taper before calling
        this.  Note the :func:`obspy.Trace.taper` functions.
    """
    return


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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
