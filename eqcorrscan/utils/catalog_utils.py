"""
Helper functions for common handling tasks for catalog objects.

.. note:: These functions are tools to aid simplification of general scripts,
    they do not cover all use cases, however if you have a use case you want
    to see here, then let the authors know, or implement it yourself and
    contribute it back to the project, or, if its really good, give it to the
    obspy guys!

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import logging

from collections import Counter
from obspy.core.event import Catalog


Logger = logging.getLogger(__name__)


def filter_picks(catalog, stations=None, channels=None, networks=None,
                 locations=None, top_n_picks=None, evaluation_mode='all',
                 phase_hints=None, enforce_single_pick=False):
    """
    Filter events in the catalog based on a number of parameters.

    :param catalog: Catalog to filter.
    :type catalog: obspy.core.event.Catalog
    :param stations: List for stations to keep picks from.
    :type stations: list
    :param channels: List of channels to keep picks from.
    :type channels: list
    :param networks: List of networks to keep picks from.
    :type networks: list
    :param locations: List of location codes to use
    :type locations: list
    :param top_n_picks: Filter only the top N most used station-channel pairs.
    :type top_n_picks: int
    :param evaluation_mode:
        To select only manual or automatic picks, or use all (default).
    :type evaluation_mode: str
    :param phase_hints: List of retained phase hints, or None to use all
    :type phase_hints: list
    :param enforce_single_pick:
        Method to enforce using only one pick of each phase-hint per
        station or False to leave all. Can be {False, "earliest"}
    :type enforce_single_pick: str


    :return:
        Filtered Catalog - if events are left with no picks, they are removed
        from the catalog.
    :rtype: obspy.core.event.Catalog

    .. note::
        Will filter first by station, then by channel, then by network, if
        using top_n_picks, this will be done last, after the other filters
        have been applied.

    .. note::
        Doesn't work in place on the catalog, your input catalog will be safe
        unless you overwrite it.

    .. note:: Doesn't expand wildcard characters.

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from eqcorrscan.utils.catalog_utils import filter_picks
    >>> from obspy import UTCDateTime
    >>> client = Client('NCEDC')
    >>> t1 = UTCDateTime(2004, 9, 28)
    >>> t2 = t1 + 86400
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3,
    ...                             minlatitude=35.7, maxlatitude=36.1,
    ...                             minlongitude=-120.6, maxlongitude=-120.2,
    ...                             includearrivals=True)
    >>> print(len(catalog))
    12
    >>> filtered_catalog = filter_picks(catalog, stations=['BMS', 'BAP',
    ...                                                    'PAG', 'PAN',
    ...                                                    'PBI', 'PKY',
    ...                                                    'YEG', 'WOF'])
    >>> print(len(filtered_catalog))
    12
    >>> stations = []
    >>> for event in filtered_catalog:
    ...     for pick in event.picks:
    ...         stations.append(pick.waveform_id.station_code)
    >>> print(sorted(list(set(stations))))
    ['BAP', 'BMS', 'PAG', 'PAN', 'PBI', 'PKY', 'WOF', 'YEG']
    """
    assert enforce_single_pick in {False, "earliest"}, \
        f"enforce_single_pick={enforce_single_pick} unknown"
    # Don't work in place on the catalog
    filtered_catalog = catalog.copy()

    if stations:
        for event in filtered_catalog:
            if len(event.picks) == 0:
                continue
            event.picks = [pick for pick in event.picks
                           if pick.waveform_id.station_code in stations]
    if channels:
        for event in filtered_catalog:
            if len(event.picks) == 0:
                continue
            event.picks = [pick for pick in event.picks
                           if pick.waveform_id.channel_code in channels]
    if networks:
        for event in filtered_catalog:
            if len(event.picks) == 0:
                continue
            event.picks = [pick for pick in event.picks
                           if pick.waveform_id.network_code in networks]
    if locations:
        for event in filtered_catalog:
            if len(event.picks) == 0:
                continue
            event.picks = [pick for pick in event.picks
                           if pick.waveform_id.location_code in locations]
    if phase_hints:
        for event in filtered_catalog:
            if len(event.picks) == 0:
                continue
            event.picks = [pick for pick in event.picks
                           if pick.phase_hint in phase_hints]

    if evaluation_mode == 'manual':
        for event in filtered_catalog:
            event.picks = [pick for pick in event.picks
                           if pick.evaluation_mode == 'manual']
    elif evaluation_mode == 'automatic':
        for event in filtered_catalog:
            event.picks = [pick for pick in event.picks
                           if pick.evaluation_mode == 'automatic']
    elif evaluation_mode != 'all':
        Logger.warning(
            'Unrecognised evaluation_mode: {0}, using all picks'.format(
                evaluation_mode))
    if top_n_picks:
        all_picks = []
        for event in filtered_catalog:
            all_picks += [(pick.waveform_id.station_code,
                           pick.waveform_id.channel_code)
                          for pick in event.picks]
        counted = Counter(all_picks).most_common()
        all_picks = []
        # Hack around sorting the counter object: Py 2 does it differently to 3
        for i in range(counted[0][1]):
            highest = [item[0] for item in counted
                       if item[1] >= counted[0][1] - i]
            # Sort them by alphabetical order in station
            highest = sorted(highest, key=lambda tup: tup[0])
            for stachan in highest:
                if stachan not in all_picks:
                    all_picks.append(stachan)
            if len(all_picks) > top_n_picks:
                all_picks = all_picks[0:top_n_picks]
                break
        for event in filtered_catalog:
            if len(event.picks) == 0:
                continue
            event.picks = [pick for pick in event.picks
                           if (pick.waveform_id.station_code,
                               pick.waveform_id.channel_code) in all_picks]
    # Remove events without picks
    tmp_catalog = Catalog()
    for event in filtered_catalog:
        if len(event.picks) > 0:
            tmp_catalog.append(event)

    # Finally remove extra picks
    if enforce_single_pick:
        reverse = False
        # TODO: Allow other options
        for ev in tmp_catalog:
            retained_picks = []
            stations = {p.waveform_id.station_code for p in ev.picks}
            for station in stations:
                phase_hints = {p.phase_hint for p in ev.picks
                               if p.waveform_id.station_code == station}
                for phase_hint in phase_hints:
                    picks = [p for p in ev.picks
                             if p.waveform_id.station_code == station
                             and p.phase_hint == phase_hint]
                    picks.sort(key=lambda p: p.time, reverse=reverse)
                    retained_picks.append(picks[0])
            ev.picks = retained_picks

    return tmp_catalog


def spatial_clip(catalog, corners, mindepth=None, maxdepth=None):
    """
    Clip the catalog to a spatial box, can be irregular.

    Can only be irregular in 2D, depth must be between bounds.

    :type catalog: :class:`obspy.core.catalog.Catalog`
    :param catalog: Catalog to clip.
    :type corners: :class:`matplotlib.path.Path`
    :param corners: Corners to clip the catalog to
    :type mindepth: float
    :param mindepth: Minimum depth for earthquakes in km.
    :type maxdepth: float
    :param maxdepth: Maximum depth for earthquakes in km.

    .. Note::
        Corners is expected to be a :class:`matplotlib.path.Path` in the form
        of tuples of (lat, lon) in decimal degrees.
    """
    cat_out = catalog.copy()
    if mindepth is not None:
        for event in cat_out:
            try:
                origin = _get_origin(event)
            except IOError:
                continue
            if origin.depth < mindepth * 1000:
                cat_out.events.remove(event)
    if maxdepth is not None:
        for event in cat_out:
            try:
                origin = _get_origin(event)
            except IOError:
                continue
            if origin.depth > maxdepth * 1000:
                cat_out.events.remove(event)
    for event in cat_out:
        try:
            origin = _get_origin(event)
        except IOError:
            continue
        if not corners.contains_point((origin.latitude, origin.longitude)):
            cat_out.events.remove(event)
    return cat_out


def _get_origin(event):
    """
    Get the origin of an event.

    :type event: :class:`obspy.core.event.Event`
    :param event: Event to get the origin of.
    :return: :class:`obspy.core.event.Origin`
    """
    if event.preferred_origin() is not None:
        origin = event.preferred_origin()
    elif len(event.origins) > 0:
        origin = event.origins[0]
    else:
        raise IndexError('No origin set, cannot constrain')
    return origin


def get_ordered_trace_indices(stream, event, sort_by="distance"):
    """
    Sort the traces of a template by hypocentral distance.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream to sort
    :type event: obspy.core.event.event.Event
    :param event:
        Event object. Only if the event contains picks, and an origin and asso-
        ciated arrivals (for sort_by="distance"), then it can sort the traces.
    :type sort_by: string
    :param sort_by: "distance" (default) or "pick_time"

    :return:
        List of indices ordered according to the traces in the sorted Stream.
    :rtype: List

    .. note::
        Indices of traces without associated sortable information will appear
        last in returned order of indices. If event doesn't contain enough
        information, None will be returned.
        For sort_by="pick_time", the earliest pick time at the station is rele-
        vant, independant of channel. To sort a template-stream by picktime,
        use stream.sort(['starttime']) instead.
    """
    from operator import itemgetter
    from obspy import UTCDateTime

    if not event:
        Logger.warning('No event information found to sort stream.')
        return None
    if not event.picks:
        Logger.warning('No picks found to sort stream.')
        return None
    origin = _get_origin(event)
    value_list = [None for i in range(len(stream))]
    if sort_by == "distance":
        for j, tr in enumerate(stream):
            # Don't change stats.distance if it is set already
            if hasattr(tr.stats, 'distance'):
                if tr.stats.distance is not None:
                    value_list[j] = tr.stats.distance
                    continue
            # Set default distance to 999 degrees.
            value_list[j] = 999
            # find the arrival that matches the trace
            for arrival in origin.arrivals:
                pick = arrival.pick_id.get_referred_object()
                if tr.stats.station == pick.waveform_id.station_code:
                    if arrival.distance is not None:
                        value_list[j] = arrival.distance
                        break
    elif sort_by == "pick_time":
        for j, tr in enumerate(stream):
            value_list[j] = UTCDateTime(9999, 12, 31, 23, 59, 59)
            for pick in event.picks:
                if tr.stats.station == pick.waveform_id.station_code and\
                    (tr.stats.channel == pick.waveform_id.channel_code or
                     (tr.id[-1] in 'NE12' and pick.waveform_id.id in 'NE12' and
                      tr.id[0:-1] == pick.waveform_id.id[0:-1])):
                    if pick.time < value_list[j]:
                        value_list[j] = pick.time
    else:
        raise NotImplementedError('Sorting by ' + sort_by + 'not implemented.')
    trace_indices_ordered, values_sorted = zip(*sorted(enumerate(value_list),
                                                       key=itemgetter(1)))
    return trace_indices_ordered


if __name__ == "__main__":
    import doctest
    doctest.testmod()
