"""
Part of the EQcorrscan package: tools to convert SAC headers to obspy event \
objects.

.. note:: This functionality is not supported for obspy versions below \
    1.0.0 as references times are not read in by SACIO, which are needed \
    for defining pick times.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import logging

import obspy
from obspy import Stream, UTCDateTime
from obspy.core.event import Event, Origin, WaveformStreamID, Pick


PICK_KEYS = ['t0', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 'a']
PHASE_KEYS = ["k{0}".format(pk) for pk in PICK_KEYS]


Logger = logging.getLogger(__name__)


def _version_check():
    if int(obspy.__version__.split('.')[0]) < 1:
        msg = 'Obspy version: ' + obspy.__version__ + ' does not have ' +\
            'correct reference time handling, please upgrade to ' +\
            'version > 1.0.0'
        raise NotImplementedError(msg)


def sactoevent(st):
    """
    Convert SAC headers (picks only) to obspy event class.

    Picks are taken from header values a, t[0-9]. A phase-hint in the
    corresponding kt[0-9] slot is recommended.

    :type st: obspy.core.stream.Stream
    :param st: Stream of waveforms including SAC headers.

    :returns: Event with picks taken from SAC headers.
    :rtype: :class:`obspy.core.event.event.Event`

    .. note::
        This functionality is not supported for obspy versions below 1.0.0 as
        reference times are not read in by SACIO, which are needed for defining
        pick times.

    .. note::
        Takes the event origin information from the first trace in the
        stream - to ensure this works as you expect, please populate the
        evla, evlo, evdp and nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec
        for all traces with the same values.

    >>> from obspy import read
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> st = read(TEST_PATH + '/SAC/2014p611252/*')
    >>> event = sactoevent(st)
    >>> print(event.origins[0].time)
    2014-08-15T03:55:21.057000Z
    >>> print(event.picks[0].phase_hint)
    S
    """
    # Check the version
    _version_check()
    # Set the default SAC nan values
    float_nan = -12345.0

    if not isinstance(st, Stream):
        msg = ('st must be a stream object, if you have just read in ' +
               'multiple SAC files you may have a list of streams, convert ' +
               'using: st = Stream([tr[0] for tr in st])')
        raise ValueError(msg)
        # Note, don't do this internally as we need to ensure that we are not
        # taking traces from other events, the user should check this.
    for tr in st:
        if not tr.stats._format == 'SAC':
            msg = ('%s.%s is not SAC formatted.' % (tr.stats.station,
                                                    tr.stats.channel))
            raise ValueError(msg)
    # Now we need to create an event!
    event = Event()
    event.origins.append(Origin())
    # print(st[0].stats.sac.keys())
    event.origins[0].time = UTCDateTime(
        year=st[0].stats.sac.nzyear, julday=st[0].stats.sac.nzjday,
        hour=st[0].stats.sac.nzhour, minute=st[0].stats.sac.nzmin,
        second=st[0].stats.sac.nzsec,
        microsecond=st[0].stats.sac.nzmsec * 1000)
    try:
        event.origins[0].latitude = st[0].stats.sac.evla
        event.origins[0].longitude = st[0].stats.sac.evlo
        event.origins[0].depth = st[0].stats.sac.evdp
        # Catch filled with 12345.0 as nan
        if event.origins[0].latitude == float_nan:
            event.origins[0].latitude = None
        if event.origins[0].longitude == float_nan:
            event.origins[0].longitude = None
        if event.origins[0].depth == float_nan:
            event.origins[0].depth = None
    except KeyError:
        event.origins[0].latitude = None
        event.origins[0].longitude = None
        event.origins[0].depth = None
    except AttributeError:
        event.origins[0].latitude = None
        event.origins[0].longitude = None
        event.origins[0].depth = None

    # Add in the picks
    for tr in st:
        reference_time = UTCDateTime(
            year=tr.stats.sac.nzyear, julday=tr.stats.sac.nzjday,
            hour=tr.stats.sac.nzhour, minute=tr.stats.sac.nzmin,
            second=tr.stats.sac.nzsec,
            microsecond=tr.stats.sac.nzmsec * 1000)
        # Possible pick locations are in the t[0-9] slot
        for pick_key, phase_key in zip(PICK_KEYS, PHASE_KEYS):
            try:
                if tr.stats.sac[pick_key] == float_nan:
                    # in version 0.10.2 and before. rather than not include
                    # the keys, the variables are filled with SAC nans.
                    Logger.debug(
                        'No pick in position {0} for trace {1}'.format(
                            pick_key, tr.id))
                    continue
                pick_time = reference_time + tr.stats.sac[pick_key]
                if phase_key in tr.stats.sac.keys():
                    phase_hint = tr.stats.sac[phase_key].split()[0]
                else:
                    Logger.warning(
                        "No phase hint found for pick in {0}".format(pick_key))
                    phase_hint = None
            except KeyError:
                Logger.debug('No pick in position {0} for trace {1}'.format(
                    pick_key, tr.id))
                continue
            Logger.info(
                'Found pick in position {0} for {1}'.format(pick_key, tr.id))
            waveform_id = WaveformStreamID(station_code=tr.stats.station,
                                           network_code=tr.stats.network,
                                           channel_code=tr.stats.channel)
            pick = Pick(waveform_id=waveform_id,
                        phase_hint=phase_hint,
                        time=pick_time)
            event.picks.append(pick)

    return event


if __name__ == '__main__':
    import doctest
    doctest.testmod()
