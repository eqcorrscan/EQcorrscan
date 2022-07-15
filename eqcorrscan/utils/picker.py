"""
Functions to pick earthquakes detected by EQcorrscan.

Designed primarily locate stacks of detections to give family locations.
Extensions may later be written, not tested for accuracy, just simple wrappers.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import logging
import numpy as np

from obspy import UTCDateTime
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.signal.filter import envelope
from obspy.core.event import Event, Pick, WaveformStreamID
from obspy.core.event import CreationInfo, Comment, Origin
from obspy.signal.trigger import classic_sta_lta, trigger_onset

import eqcorrscan.utils.plotting as plotting


Logger = logging.getLogger(__name__)


def cross_net(stream, env=False, master=False):
    """
    Generate picks using a simple envelope cross-correlation.

    Picks are made for each channel based on optimal moveout defined by
    maximum cross-correlation with master trace.  Master trace will be the
    first trace in the stream if not set.  Requires good inter-station
    coherance.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream to pick
    :type env: bool
    :param env: To compute cross-correlations on the envelope or not.
    :type master: obspy.core.trace.Trace
    :param master:
        Trace to use as master, if False, will use the first trace in stream.

    :returns: :class:`obspy.core.event.event.Event`

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.picker import cross_net
    >>> st = read()
    >>> event = cross_net(st, env=True)
    >>> print(event.creation_info.author)
    EQcorrscan

    .. warning::
        This routine is not designed for accurate picking, rather it can be
        used for a first-pass at picks to obtain simple locations. Based on
        the waveform-envelope cross-correlation method.
    """
    event = Event()
    event.origins.append(Origin())
    event.creation_info = CreationInfo(author='EQcorrscan',
                                       creation_time=UTCDateTime())
    event.comments.append(Comment(text='cross_net'))
    samp_rate = stream[0].stats.sampling_rate
    if not env:
        Logger.info('Using the raw data')
        st = stream.copy()
        st.resample(samp_rate)
    else:
        st = stream.copy()
        Logger.info('Computing envelope')
        for tr in st:
            tr.resample(samp_rate)
            tr.data = envelope(tr.data)
    if not master:
        master = st[0]
    else:
        master = master
    master.data = np.nan_to_num(master.data)
    for i, tr in enumerate(st):
        tr.data = np.nan_to_num(tr.data)
        Logger.debug('Comparing {0} with the master'.format(tr.id))
        shift_len = int(0.3 * len(tr))
        Logger.debug('Shift length is set to ' + str(shift_len) + ' samples')
        corr_fun = correlate(master, tr, shift_len)
        index, cc = xcorr_max(corr_fun)
        wav_id = WaveformStreamID(station_code=tr.stats.station,
                                  channel_code=tr.stats.channel,
                                  network_code=tr.stats.network)
        event.picks.append(Pick(time=tr.stats.starttime +
                                (index / tr.stats.sampling_rate),
                                waveform_id=wav_id,
                                phase_hint='S',
                                onset='emergent'))
        Logger.debug(event.picks[i])
    event.origins[0].time = min([pick.time for pick in event.picks]) - 1
    # event.origins[0].latitude = float('nan')
    # event.origins[0].longitude = float('nan')
    # Set arbitrary origin time
    del st
    return event


def stalta_pick(stream, stalen, ltalen, trig_on, trig_off, freqmin=False,
                freqmax=False, show=False):
    """
    Basic sta/lta picker, suggest using alternative in obspy.

    Simple sta/lta (short-term average/long-term average) picker, using
    obspy's :func:`obspy.signal.trigger.classic_sta_lta` routine to generate
    the characteristic function.

    Currently very basic quick wrapper, there are many other (better) options
    in obspy in the :mod:`obspy.signal.trigger` module.

    :type stream: obspy.core.stream.Stream
    :param stream: The stream to pick on, can be any number of channels.
    :type stalen: float
    :param stalen: Length of the short-term average window in seconds.
    :type ltalen: float
    :param ltalen: Length of the long-term average window in seconds.
    :type trig_on: float
    :param trig_on: sta/lta ratio to trigger a detection/pick
    :type trig_off: float
    :param trig_off: sta/lta ratio to turn the trigger off - no further picks\
        will be made between exceeding trig_on until trig_off is reached.
    :type freqmin: float
    :param freqmin: Low-cut frequency in Hz for bandpass filter
    :type freqmax: float
    :param freqmax: High-cut frequency in Hz for bandpass filter
    :type show: bool
    :param show: Show picks on waveform.

    :returns: :class:`obspy.core.event.event.Event`

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.picker import stalta_pick
    >>> st = read()
    >>> event = stalta_pick(st, stalen=0.2, ltalen=4, trig_on=10,
    ...             trig_off=1, freqmin=3.0, freqmax=20.0)
    >>> print(event.creation_info.author)
    EQcorrscan

    .. warning::
        This function is not designed for accurate picking, rather it can give
        a first idea of whether picks may be possible.  Proceed with caution.
    """
    event = Event()
    event.origins.append(Origin())
    event.creation_info = CreationInfo(author='EQcorrscan',
                                       creation_time=UTCDateTime())
    event.comments.append(Comment(text='stalta'))
    picks = []
    for tr in stream:
        # We are going to assume, for now, that if the pick is made on the
        # horizontal channel then it is an S, otherwise we will assume it is
        # a P-phase: obviously a bad assumption...
        if tr.stats.channel[-1] == 'Z':
            phase = 'P'
        else:
            phase = 'S'
        if freqmin and freqmax:
            tr.detrend('simple')
            tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                      corners=3, zerophase=True)
        df = tr.stats.sampling_rate
        cft = classic_sta_lta(tr.data, int(stalen * df), int(ltalen * df))
        triggers = trigger_onset(cft, trig_on, trig_off)
        for trigger in triggers:
            on = tr.stats.starttime + (trigger[0] / df)
            # off = tr.stats.starttime + (trigger[1] / df)
            wav_id = WaveformStreamID(station_code=tr.stats.station,
                                      channel_code=tr.stats.channel,
                                      network_code=tr.stats.network)
            p = Pick(waveform_id=wav_id, phase_hint=phase, time=on)
            Logger.info('Pick made: {0}'.format(p))
            picks.append(p)
    # QC picks
    pick_stations = list(set([pick.waveform_id.station_code
                              for pick in picks]))
    for pick_station in pick_stations:
        station_picks = [pick for pick in picks if
                         pick.waveform_id.station_code == pick_station]
        # If P-pick is after S-picks, remove it.
        p_time = [pick.time for pick in station_picks
                  if pick.phase_hint == 'P']
        s_time = [pick.time for pick in station_picks
                  if pick.phase_hint == 'S']
        if p_time > s_time:
            p_pick = [pick for pick in station_picks if pick.phase_hint == 'P']
            for pick in p_pick:
                Logger.info('P pick after S pick, removing P pick')
                picks.remove(pick)
    event.picks = picks
    if show:
        plotting.pretty_template_plot(stream, event=event, title='Autopicks',
                                      size=(8, 9))
    if len(event.picks) > 0:
        event.origins[0].time = min([pick.time for pick in event.picks]) - 1
        # event.origins[0].latitude = float('nan')
        # event.origins[0].longitude = float('nan')
    # Set arbitrary origin time
    return event


if __name__ == "__main__":
    import doctest
    doctest.testmod()
