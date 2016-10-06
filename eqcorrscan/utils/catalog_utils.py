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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from obspy.core.event import Catalog


def filter_picks(catalog, stations=None, channels=None, networks=None,
                 locations=None, top_n_picks=None, evaluation_mode='all'):
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
    if evaluation_mode == 'manual':
        for event in filtered_catalog:
            event.picks = [pick for pick in event.picks
                           if pick.evaluation_mode == 'manual']
    elif evaluation_mode == 'automatic':
        for event in filtered_catalog:
            event.picks = [pick for pick in event.picks
                           if pick.evaluation_mode == 'automatic']
    elif evaluation_mode != 'all':
        warnings.warn('Unrecognised evaluation_mode: %s, using all picks' %
                      evaluation_mode)
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

    return tmp_catalog


if __name__ == "__main__":
    import doctest
    doctest.testmod()
