r"""Functions to locate earthquakes detected by EQcorrscan.
Designed primarily locate stacks of detections to give family locations.
Extensions may later be written.

Copyright 2015 Calum Chamberlain

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


def synth_compare(stream, stream_list, cores=4, debug=0):
    r"""Compare a specific stream to a list of synthetic templates, or \
    earthquakes of known source and find the best matching event.

    :type stream: :class:obspy.Stream
    :param stream: Stream to be compared to streams with known locations.
    :type stream_list: list
    :param stream_list: List of streams with known locations
    :type cores: int
    :param cores: Number of cores to parallel over
    :type debug: int
    :param debug: Debug level, high is more debug

    :returns: int, float: index of best match and cross-correlation sum
    """

    from eqcorrscan.core.match_filter import _channel_loop
    import numpy as np
    import copy
    from obspy import Trace

    stream_copy = stream.copy()
    templates = copy.deepcopy(stream_list)
    # Need to fill the stream_list - template - channels
    template_stachan = []
    for template in templates:
        for tr in template:
            template_stachan += [(tr.stats.station, tr.stats.channel)]
    template_stachan = list(set(template_stachan))

    for stachan in template_stachan:
        if not stream_copy.select(station=stachan[0], channel=stachan[1]):
            # Remove template traces rather than adding NaN data
            for template in templates:
                if template.select(station=stachan[0], channel=stachan[1]):
                    for tr in template.select(station=stachan[0],
                                              channel=stachan[1]):
                        template.remove(tr)
    # Remove un-needed channels
    for tr in stream_copy:
        if not (tr.stats.station, tr.stats.channel) in template_stachan:
            stream_copy.remove(tr)
    # Also pad out templates to have all channels
    for template in templates:
        for stachan in template_stachan:
            if not template.select(station=stachan[0], channel=stachan[1]):
                nulltrace = Trace()
                nulltrace.stats.station = stachan[0]
                nulltrace.stats.channel = stachan[1]
                nulltrace.stats.sampling_rate = template[0].stats.sampling_rate
                nulltrace.stats.starttime = template[0].stats.starttime
                nulltrace.data = np.array([np.NaN] * len(template[0].data),
                                          dtype=np.float32)
                template += nulltrace
    # Hand off  cross-correaltion to _channel_loop, which runs in parallel
    [cccsums, no_chans] = _channel_loop(templates, stream_copy, cores, debug)
    cccsums = [np.max(cccsum) for cccsum in cccsums]
    # Find the maximum cccsum and index thereof
    index = np.argmax(cccsums)
    cccsum = cccsums[index]
    return index, cccsum


def cross_net(stream, env=False, debug=0, master=False):
    r"""Function to generate picks for each channel based on optimal moveout \
    defined by maximum cross-correaltion with master trace.  Master trace \
    will be the first trace in the stream.

    :type stream: :class: obspy.Stream
    :param stream: Stream to pick
    :type envelope: bool
    :param envelope: To compute cross-correlations on the envelope or not.
    :type debug: int
    :param debug: Debug level from 0-5
    :type master: obspy.Trace
    :param master: Trace to use as master, if False, will use the first trace \
            in stream.

    :returns: list of pick class
    """
    from obspy.signal.cross_correlation import xcorr
    from obspy.signal.filter import envelope
    from eqcorrscan.utils.Sfile_util import PICK
    import matplotlib.pyplot as plt
    import numpy as np
    picks = []
    samp_rate = stream[0].stats.sampling_rate
    if not env:
        if debug > 2:
            print 'Using the raw data'
        st = stream.copy()
        st.resample(samp_rate)
    else:
        st = stream.copy()
        if debug > 2:
            print 'Computing envelope'
        for tr in st:
            tr.resample(samp_rate)
            tr.data = envelope(tr.data)
    if debug > 2:
        st.plot(equal_scale=False, size=(800, 600))
    if not master:
        master = st[0]
    else:
        master = master
    master.data = np.nan_to_num(master.data)
    for tr in st:
        tr.data = np.nan_to_num(tr.data)
        if debug > 2:
            msg = ' '.join(['Comparing', tr.stats.station, tr.stats.channel,
                            'with the master'])
            print msg
        shift_len = int(0.3 * len(tr))
        if debug > 2:
            print 'Shift length is set to ' + str(shift_len) + ' samples'
        if debug > 3:
            index, cc, cc_vec = xcorr(master, tr, shift_len, full_xcorr=True)
            cc_vec = np.nan_to_num(cc_vec)
            if debug > 4:
                print cc_vec
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            x = np.linspace(0, len(master) / samp_rate,
                            len(master))
            ax1.plot(x, master.data / float(master.data.max()), 'k',
                     label='Master')
            ax1.plot(x + (index / samp_rate), tr.data / float(tr.data.max()),
                     'r', label='Slave shifted')
            ax1.legend(loc="lower right", prop={'size': "small"})
            ax1.set_xlabel("time [s]")
            ax1.set_ylabel("norm. amplitude")
            ax2 = fig.add_subplot(212)
            print len(cc_vec)
            x = np.linspace(0, len(cc_vec) / samp_rate, len(cc_vec))
            ax2.plot(x, cc_vec, label='xcorr')
            # ax2.set_ylim(-1, 1)
            # ax2.set_xlim(0, len(master))
            plt.show()
        index, cc = xcorr(master, tr, shift_len)
        pick = PICK(station=tr.stats.station,
                    channel=tr.stats.channel,
                    impulsivity='E',
                    phase='S',
                    weight='1',
                    time=tr.stats.starttime + (index / tr.stats.sampling_rate))
        if debug > 2:
            print pick
        picks.append(pick)
    del st
    return picks


def stalta_pick(stream, stalen, ltalen, trig_on, trig_off, freqmin=False,
                freqmax=False, debug=0, show=False):
    r"""Simple sta-lta (short-term average/long-term average) picker, using \
    obspy's stalta routine to generate the characteristic function.

    Currently very basic quick wrapper, there are many other (better) options \
    in obspy, found \
    (here)[http://docs.obspy.org/packages/autogen/obspy.signal.trigger.html].

    :type stream: obspy.Stream
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
    :type debug: int
    :param debug: Debug output level from 0-5.
    :type show: bool
    :param show: Show picks on waveform.

    :returns: list of pick class.
    """
    from obspy.signal.trigger import classicSTALTA, triggerOnset, plotTrigger
    from Sfile_util import PICK
    import EQcorrscan_plotting as plotting
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
        cft = classicSTALTA(tr.data, int(stalen * df), int(ltalen * df))
        if debug > 3:
            plotTrigger(tr, cft, trig_on, trig_off)
        triggers = triggerOnset(cft, trig_on, trig_off)
        for trigger in triggers:
            on = tr.stats.starttime + (trigger[0] / df)
            # off = tr.stats.starttime + (trigger[1] / df)
            pick = PICK(station=tr.stats.station, channel=tr.stats.channel,
                        time=on, phase=phase)
            if debug > 2:
                print 'Pick made:'
                print pick
            picks.append(pick)
    # QC picks
    del pick
    pick_stations = list(set([pick.station for pick in picks]))
    for pick_station in pick_stations:
        station_picks = [pick for pick in picks if
                         pick.station == pick_station]
        # If P-pick is after S-picks, remove it.
        p_time = [pick.time for pick in station_picks if pick.phase == 'P']
        s_time = [pick.time for pick in station_picks if pick.phase == 'S']
        if p_time > s_time:
            p_pick = [pick for pick in station_picks if pick.phase == 'P']
            for pick in p_pick:
                print('P pick after S pick, removing P pick')
                picks.remove(pick)
    if show:
        plotting.pretty_template_plot(stream, picks=picks, title='Autopicks',
                                      size=(8, 9))
    return picks


if __name__ == "__main__":
    import doctest
    doctest.testmod()
