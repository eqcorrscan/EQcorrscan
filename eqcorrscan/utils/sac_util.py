"""
Part of the EQcorrscan package: tools to convert SAC headers to obspy event \
objects.

Authors: Calum Chamberlain and the EQcorrscan developers.

.. note:: This functionality is not supported for obspy versions below \
    1.0.0 as references times are not read in by SACIO, which are needed \
    for defining pick times.

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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import warnings


def _version_check():
    import obspy
    if int(obspy.__version__.split('.')[0]) < 1:
        msg = 'Obspy version: ' + obspy.__version__ + ' does not have ' +\
            'correct reference time handling, please upgrade to ' +\
            'version > 1.0.0'
        raise NotImplementedError(msg)


def sactoevent(st, debug=0):
    """
    Function to convert SAC headers to obspy event class.

    :type st: obspy.core.Stream
    :param st: Stream of waveforms including SAC headers.
    :type debug: int
    :pram debug: Debug level, larger number = more output.

    :returns: obspy.core.Event

    .. note:: This functionality is not supported for obspy versions below \
        1.0.0 as references times are not read in by SACIO, which are needed \
        for defining pick times.
    """
    from obspy.core.event import Event, Origin, WaveformStreamID, Pick
    from obspy import Stream, UTCDateTime
    import numpy as np

    # Check the version
    _version_check()
    # Set the default SAC nan values
    float_nan = -12345.0

    if not isinstance(st, Stream):
        msg = 'st must be a stream object, if you have just read in ' +\
            'multiple SAC files you may have a list of streams, convert ' +\
            'using: st = Stream([tr[0] for tr in st])'
        raise IOError(msg)
        # Note, don't do this internally as we need to ensure that we are not
        # taking traces from other events, the user should check this.
    for tr in st:
        if not tr.stats._format == 'SAC':
            msg = tr.stats.station + '.' + tr.stats.channel + ' is not SAC ' +\
                'formatted'
            raise IOError(msg)
    # Now we need to create an event!
    event = Event()
    event.origins.append(Origin())
    print(st[0].stats.sac.keys())
    event.origins[0].time = UTCDateTime(year=st[0].stats.sac.nzyear,
                                        julday=st[0].stats.sac.nzjday,
                                        hour=st[0].stats.sac.nzhour,
                                        minute=st[0].stats.sac.nzmin,
                                        second=st[0].stats.sac.nzsec,
                                        microsecond=st[0].stats.
                                        sac.nzmsec * 1000)
    try:
        event.origins[0].latitude = st[0].stats.sac.evla
        event.origins[0].longitude = st[0].stats.sac.evlo
        event.origins[0].depth = st[0].stats.sac.evdp
        # Catch filled with 12345.0 as nan
        if event.origins[0].latitude == float_nan:
            event.origins[0].latitude = np.nan
        if event.origins[0].longitude == float_nan:
            event.origins[0].longitude = np.nan
        if event.origins[0].depth == float_nan:
            event.origins[0].depth = np.nan
    except KeyError:
        event.origins[0].latitude = np.nan
        event.origins[0].longitude = np.nan
        event.origins[0].depth = np.nan

    # Add in the picks
    for tr in st:
        reference_time = UTCDateTime(year=tr.stats.sac.nzyear,
                                     julday=tr.stats.sac.nzjday,
                                     hour=tr.stats.sac.nzhour,
                                     minute=tr.stats.sac.nzmin,
                                     second=tr.stats.sac.nzsec,
                                     microsecond=tr.stats.sac.nzmsec * 1000)
        # Possible pick locations are in the t[0-9] slot
        for pick_number in range(10):
            pick_key = 't' + str(pick_number)
            phase_key = 'kt' + str(pick_number)
            try:
                if tr.stats.sac[pick_key] == float_nan:
                    # in version 0.10.2 and before. rather than not include
                    # the keys, the variables are filled with SAC nans.
                    if debug > 1:
                        msg = 'No pick in position ' + pick_key + \
                            ' for trace: ' + tr.stats.station + '.' + \
                            tr.stats.channel
                        warnings.warn(msg)
                    continue
                pick_time = reference_time + tr.stats.sac[pick_key]
                phase_hint = tr.stats.sac[phase_key].split()[0]
            except KeyError:
                if debug > 1:
                    msg = 'No pick in position ' + pick_key + ' for trace: ' +\
                        tr.stats.station + '.' + tr.stats.channel
                    warnings.warn(msg)
                continue
            waveform_id = WaveformStreamID(station_code=tr.stats.station,
                                           network_code=tr.stats.network,
                                           channel_code=tr.stats.channel)
            pick = Pick(waveform_id=waveform_id,
                        phase_hint=phase_hint,
                        time=pick_time)
            event.picks.append(pick)

    return event
