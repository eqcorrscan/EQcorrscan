#!/usr/bin/python
"""
Utility code for most of the plots used as part of the EQcorrscan package.

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
import numpy as np
import matplotlib.pylab as plt


def chunk_data(tr, samp_rate, state='mean'):
    r"""Function to downsample data for plotting by computing the maximum of \
    data within chunks, useful for plotting waveforms or cccsums, large \
    datasets that would otherwise exceed the complexity allowed, and overflow.

    :type tr: obspy.Trace
    :param tr: Trace to be chunked
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type state: str
    :param state: Either 'Min', 'Max', 'Mean' or 'Maxabs' to return one of \
        these for the chunks. Maxabs will return the largest (positive or \
        negative) for that chunk.

    :returns: :class: obspy.Trace
    """
    trout = tr.copy()  # Don't do it inplace on data
    x = np.arange(len(tr.data))
    y = tr.data

    chunksize = int(round(tr.stats.sampling_rate / samp_rate))
    # Wrap the array into a 2D array of chunks, truncating the last chunk if
    # chunksize isn't an even divisor of the total size.
    # (This part won't use _any_ additional memory)
    numchunks = y.size // chunksize
    ychunks = y[:chunksize*numchunks].reshape((-1, chunksize))
    xchunks = x[:chunksize*numchunks].reshape((-1, chunksize))

    # Calculate the max, min, and means of chunksize-element chunks...
    if state == 'Max':
        trout.data = ychunks.max(axis=1)
    elif state == 'Min':
        trout.data = ychunks.min(axis=1)
    elif state == 'Mean':
        trout.data = ychunks.mean(axis=1)
    elif state == 'Maxabs':
        max_env = ychunks.max(axis=1)
        min_env = ychunks.min(axis=1)
        indeces = np.argmax(np.vstack([np.abs(max_env), np.abs(min_env)]),
                            axis=0)
        stack = np.vstack([max_env, min_env]).T
        trout.data = np.array([stack[i][indeces[i]]
                              for i in xrange(len(stack))])
    xcenters = xchunks.mean(axis=1)
    trout.stats.starttime = tr.stats.starttime + xcenters[0] /\
        tr.stats.sampling_rate
    trout.stats.sampling_rate = samp_rate
    return trout


def triple_plot(cccsum, cccsum_hist, trace, threshold, save=False,
                savefile=''):
    r"""Main function to make a triple plot with a day-long seismogram, \
    day-long correlation sum trace and histogram of the correlation sum to \
    show normality.

    :type cccsum: numpy.ndarray
    :param cccsum: Array of the cross-channel cross-correlation sum
    :type cccsum_hist: numpy.ndarray
    :param cccsum_hist: cccsum for histogram plotting, can be the same as \
        cccsum but included if cccsum is just an envelope.
    :type trace: obspy.Trace
    :param trace: A sample trace from the same time as cccsum
    :type threshold: float
    :param threshold: Detection threshold within cccsum
    :type save: bool, optional
    :param save: If True will svae and not plot to screen, vice-versa if False
    :type savefile: str, optional
    :param savefile: Path to save figure to, only required if save=True
    """
    if len(cccsum) != len(trace.data):
        print 'cccsum is: '+str(len(cccsum))+' trace is: '+str(len(trace.data))
        msg = ' '.join(['cccsum and trace must have the',
                        'same number of data points'])
        raise ValueError(msg)
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / (df * 3600)
    # Generate the subplot for the seismic data
    ax1 = plt.subplot2grid((2, 5), (0, 0), colspan=4)
    ax1.plot(t, trace.data, 'k')
    ax1.axis('tight')
    ax1.set_ylim([-15 * np.mean(np.abs(trace.data)),
                  15 * np.mean(np.abs(trace.data))])
    # Generate the subplot for the correlation sum data
    ax2 = plt.subplot2grid((2, 5), (1, 0), colspan=4, sharex=ax1)
    # Plot the threshold values
    ax2.plot([min(t), max(t)], [threshold, threshold], color='r', lw=1,
             label="Threshold")
    ax2.plot([min(t), max(t)], [-threshold, -threshold], color='r', lw=1)
    ax2.plot(t, cccsum, 'k')
    ax2.axis('tight')
    ax2.set_ylim([-1.7 * threshold, 1.7 * threshold])
    ax2.set_xlabel("Time after %s [hr]" % trace.stats.starttime.isoformat())
    # ax2.legend()
    # Generate a small subplot for the histogram of the cccsum data
    ax3 = plt.subplot2grid((2, 5), (1, 4), sharey=ax2)
    ax3.hist(cccsum_hist, 200, normed=1, histtype='stepfilled',
             orientation='horizontal', color='black')
    ax3.set_ylim([-5, 5])
    fig = plt.gcf()
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(savefile)
    return


def peaks_plot(data, starttime, samp_rate, save=False, peaks=[(0, 0)],
               savefile=''):
    r"""Simple utility code to plot the correlation peaks to check that the \
    peak finding routine is running correctly, used in debugging for the \
    EQcorrscan module.

    :type data: numpy.array
    :param data: Numpy array of the data within which peaks have been found
    :type starttime: obspy.UTCDateTime
    :param starttime: Start time for the data
    :type samp_rate: float
    :param samp_rate: Sampling rate of data in Hz
    :type save: Boolean, optional
    :param save: Save figure or plot to screen (False)
    :type peaks: list of Tuple, optional
    :param peaks: List of peak locations and amplitudes (loc, amp)
    :type savefile: String, optional
    :param savefile: Path to save to, only used if save=True
    """
    npts = len(data)
    t = np.arange(npts, dtype=np.float32) / (samp_rate * 3600)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, data, 'k')
    ax1.scatter(peaks[0][1] / (samp_rate * 3600), abs(peaks[0][0]),
                color='r', label='Peaks')
    for peak in peaks:
        ax1.scatter(peak[1] / (samp_rate * 3600), abs(peak[0]), color='r')
    ax1.legend()
    ax1.set_xlabel("Time after %s [hr]" % starttime.isoformat())
    ax1.axis('tight')
    fig.suptitle('Peaks')
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(savefile)
    return


def cumulative_detections(dates, template_names, save=False, savefile=''):
    r"""Simple plotting function to take a list of datetime objects and plot \
    a cumulative detections list.  Can take dates as a list of lists and will \
    plot each list seperately, e.g. if you have dates from more than one \
    template it will overlay them in different colours.

    :type dates: list of lists of datetime.datetime
    :param dates: Must be a list of lists of datetime.datetime objects
    :type template_names: list of strings
    :param template_names: List of the template names in order of the dates
    :type save: bool
    :param save: Save figure or show to screen, optional
    :type savefile: str
    :param savefile: String to save to, optional
    """
    # Set up a default series of parameters for lines
    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black',
              'firebrick', 'purple', 'darkgoldenrod', 'gray']
    linestyles = ['-', '-.',  '--',  ':']
    # Check that dates is a list of lists
    if type(dates[0]) != list:
        dates = [dates]
    i = 0
    j = 0
    # This is an ugly way of looping through colours and linestyles, it would
    # be better with itertools functions...
    plothandles = []
    for k, template_dates in enumerate(dates):
        template_dates.sort()
        counts = np.arange(0, len(template_dates))
        print str(i)+' '+str(j)+' '+str(k)
        filename = plt.plot(template_dates, counts, linestyles[j],
                            color=colors[i], label=template_names[k],
                            linewidth=3.0)
        plothandles.append(filename)
        if i < len(colors) - 1:
            i += 1
        else:
            i = 0
            if j < len(linestyles) - 1:
                j += 1
            else:
                j = 0
    plt.xlabel('Date')
    plt.ylabel('Cumulative detections')
    plt.title('Cumulative detections for all templates')
    plt.legend(loc=2, prop={'size': 8}, ncol=2)  # handles=plothandles)
    if save:
        plt.savefig(savefile)
        plt.close()
    else:
        plt.show()
    return


def threeD_gridplot(nodes, save=False, savefile=''):
    r"""Function to plot in 3D a series of grid points.

    :type nodes: list of tuples
    :param nodes: List of tuples of the form (lat, long, depth)
    :type save: bool
    :param save: if True will save without plotting to screen, if False \
        (default) will plot to screen but not save
    :type savefile: str
    :param savefile: required if save=True, path to save figure to.
    """
    lats = []
    longs = []
    depths = []
    for node in nodes:
        lats.append(float(node[0]))
        longs.append(float(node[1]))
        depths.append(float(node[2]))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(lats, longs, depths)
    ax.set_ylabel("Latitude (deg)")
    ax.set_xlabel("Longitude (deg)")
    ax.set_zlabel("Depth(km)")
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(savefile)
    return


def multi_event_singlechan(streams, catalog, clip=10.0, pre_pick=2.0,
                           freqmin=False, freqmax=False, realign=False,
                           cut=(-3.0, 5.0), PWS=False, title=False):
    r"""Function to plot data from a single channel at a single station for \
    multiple events - data will be alligned by their pick-time given in the \
    picks.

    :type streams: list of :class:obspy.stream
    :param streams: List of the streams to use, can contain more traces than \
        you plan on plotting
    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog of events, one for each trace, with a single pick
    :type clip: float
    :param clip: Length in seconds to plot, defaults to 10.0
    :type pre_pick: float
    :param pre_pick: Length in seconds to extract and plot before the pick, \
        defaults to 2.0
    :type freqmin: float
    :param freqmin: Low cut for bandpass in Hz
    :type freqmax: float
    :param freqmax: High cut for bandpass in Hz
    :type realign: bool
    :param realign: To compute best alignement based on correlation or not.
    :type cut: tuple
    :param cut: tuple of start and end times for cut in seconds from the pick
    :type PWS: bool
    :param PWS: compute Phase Weighted Stack, if False, will compute linear \
        stack.
    :type title: str
    :param title: Plot title.

    :returns: Alligned and cut traces, and new picks
    """
    from eqcorrscan.utils import stacking
    import copy
    from eqcorrscan.core.match_filter import normxcorr2
    from obspy import Stream
    import warnings
    fig, axes = plt.subplots(len(catalog)+1, 1, sharex=True, figsize=(7, 12))
    axes = axes.ravel()
    traces = []
    al_traces = []
    # Keep input safe
    clist = copy.deepcopy(catalog)
    st_list = copy.deepcopy(streams)
    for i, event in enumerate(clist):
        if st_list[i].select(station=event.picks[0].waveform_id.station_code,
                             channel='*' +
                             event.picks[0].waveform_id.channel_code[-1]):
            tr = st_list[i].select(station=event.picks[0].waveforms_id.
                                   station_code,
                                   channel='*' +
                                   event.picks[0].waveform_id.
                                   channel_code[-1])[0]
        else:
            print 'No data for '+event.pick[0].waveform_id
            continue
        tr.detrend('linear')
        if freqmin:
            tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        if realign:
            tr_cut = tr.copy()
            tr_cut.trim(event.picks[0].time + cut[0],
                        event.picks[0].time + cut[1],
                        nearest_sample=False)
            if len(tr_cut.data) <= (0.5 * (cut[1] - cut[0]) *
                                    tr_cut.stats.sampling_rate):
                msg = ''.join(['Not enough in the trace for ',
                               tr.stats.station,
                               '.', tr.stats.channel, '\n',
                               'Suggest removing pick from sfile at time ',
                               str(event.picks[0].time)])
                warnings.warn(msg)
            else:
                al_traces.append(tr_cut)
        else:
            tr.trim(event.picks[0].time - pre_pick,
                    event.picks[0].time + clip - pre_pick,
                    nearest_sample=False)
        if len(tr.data) == 0:
            msg = ''.join(['No data in the trace for ', tr.stats.station,
                           '.', tr.stats.channel, '\n',
                           'Suggest removing pick from sfile at time ',
                           str(event.picks[0].time)])
            warnings.warn(msg)
            continue
        traces.append(tr)
    if realign:
        shift_len = int(0.25 * (cut[1] - cut[0]) *
                        al_traces[0].stats.sampling_rate)
        shifts = stacking.align_traces(al_traces, shift_len)
        for i in xrange(len(shifts)):
            print 'Shifting by '+str(shifts[i])+' seconds'
            event.picks[0].time -= shifts[i]
            traces[i].trim(event.picks[0].time - pre_pick,
                           event.picks[0].time + clip-pre_pick,
                           nearest_sample=False)
    # We now have a list of traces
    traces = [(trace, trace.stats.starttime.datetime) for trace in traces]
    traces.sort(key=lambda tup: tup[1])
    traces = [trace[0] for trace in traces]
    # Plot the traces
    for i, tr in enumerate(traces):
        y = tr.data
        x = np.arange(len(y))
        x = x / tr.stats.sampling_rate  # convert to seconds
        axes[i+1].plot(x, y, 'k', linewidth=1.1)
        axes[i+1].yaxis.set_ticks([])
    traces = [Stream(trace) for trace in traces]
    if PWS:
        linstack = stacking.PWS_stack(traces)
    else:
        linstack = stacking.linstack(traces)
    tr = linstack.select(station=event[0].picks[0].waveform_id.station_code,
                         channel='*' +
                         event[0].picks[0].waveform_id.channel_code[-1])[0]
    y = tr.data
    x = np.arange(len(y))
    x = x / tr.stats.sampling_rate
    axes[0].plot(x, y, 'r', linewidth=2.0)
    axes[0].set_ylabel('Stack', rotation=0)
    axes[0].yaxis.set_ticks([])
    for i, slave in enumerate(traces):
        cc = normxcorr2(tr.data, slave[0].data)
        axes[i+1].set_ylabel('cc='+str(round(np.max(cc), 2)), rotation=0)
        axes[i+1].text(0.9, 0.15, str(round(np.max(slave[0].data))),
                       bbox=dict(facecolor='white', alpha=0.95),
                       transform=axes[i+1].transAxes)
        axes[i+1].text(0.7, 0.85, slave[0].stats.starttime.datetime.
                       strftime('%Y/%m/%d %H:%M:%S'),
                       bbox=dict(facecolor='white', alpha=0.95),
                       transform=axes[i+1].transAxes)
    axes[-1].set_xlabel('Time (s)')
    if title:
        axes[0].set_title(title)
    plt.subplots_adjust(hspace=0)
    plt.show()
    return traces, clist


def detection_multiplot(stream, template, times, streamcolour='k',
                        templatecolour='r'):
    r"""Function to plot the stream of data that has been detected in, with\
    the template on top of it timed according to a list of given times, just\
    a pretty way to show a detection!

    :type stream: obspy.Stream
    :param stream: Stream of data to be plotted as the base (black)
    :type template: obspy.Stream
    :param template: Template to be plotted on top of the base stream (red)
    :type times: list of datetime.datetime
    :param times: list of times of detections in the order of the channels in
                template.
    :type streamcolour: str
    :param streamcolour: String of matplotlib colour types for the stream
    :type templatecolour: str
    :param templatecolour: Colour to plot the template in.
    """
    import datetime as dt
    from obspy import UTCDateTime
    fig, axes = plt.subplots(len(template), 1, sharex=True)
    axes = axes.ravel()
    mintime = min([tr.stats.starttime for tr in template])
    for i, template_tr in enumerate(template):
        image = stream.select(station=template_tr.stats.station,
                              channel='*'+template_tr.stats.channel[-1])
        if not image:
            msg = ' '.join(['No data for', template_tr.stats.station,
                            template_tr.stats.channel])
            print(msg)
            continue
        image = image.merge()[0]
        # Downsample if needed
        if image.stats.sampling_rate > 20:
            image.decimate(int(image.stats.sampling_rate/20))
        if template_tr.stats.sampling_rate > 20:
            template_tr.decimate(int(template_tr.stats.sampling_rate/20))
        # Get a list of datetime objects
        image_times = [image.stats.starttime.datetime +
                       dt.timedelta((j * image.stats.delta) / 86400)
                       for j in range(len(image.data))]
        axes[i].plot(image_times, image.data / max(image.data),
                     streamcolour, linewidth=1.2)
        for k, time in enumerate(times):
            lagged_time = UTCDateTime(time) + (template_tr.stats.starttime -
                                               mintime)
            lagged_time = lagged_time.datetime
            template_times = [lagged_time +
                              dt.timedelta((j * template_tr.stats.delta) /
                                           86400)
                              for j in range(len(template_tr.data))]
            axes[i].plot(template_times,
                         template_tr.data / max(template_tr.data),
                         templatecolour, linewidth=1.2)
        ylab = '.'.join([template_tr.stats.station,
                         template_tr.stats.channel])
        axes[i].set_ylabel(ylab, rotation=0,
                           horizontalalignment='right')
        axes[i].yaxis.set_ticks([])
    axes[len(axes) - 1].set_xlabel('Time')
    plt.subplots_adjust(hspace=0, left=0.175, right=0.95, bottom=0.07)
    plt.show()
    return


def interev_mag_sfiles(sfiles):
    r"""Function to plot interevent-time versus magnitude for series of events.
    **thin** Wrapper for interev_mag.

    :type sfiles: list
    :param sfiles: List of sfiles to read from
    """
    from eqcorrscan.utils import Sfile_util
    times = [Sfile_util.readheader(sfile)[0].origins[0].time
             for sfile in sfiles]
    mags = [Sfile_util.readheader(sfile)[0].magnitudes[0].mag
            for sfile in sfiles]
    interev_mag(times, mags)


def interev_mag(times, mags):
    r"""Function to plot interevent times against magnitude for given times
    and magnitudes.

    :type times: list of datetime
    :param times: list of the detection times, must be sorted the same as mags
    :type mags: list of float
    :param mags: list of magnitudes
    """
    l = [(times[i], mags[i]) for i in xrange(len(times))]
    l.sort(key=lambda tup: tup[0])
    times = [x[0] for x in l]
    mags = [x[1] for x in l]
    # Make two subplots next to each other of time before and time after
    fig, axes = plt.subplots(1, 2, sharey=True)
    axes = axes.ravel()
    pre_times = []
    post_times = []
    for i in range(len(times)):
        if i > 0:
            pre_times.append((times[i] - times[i - 1]) / 60)
        if i < len(times) - 1:
            post_times.append((times[i + 1] - times[i]) / 60)
    axes[0].scatter(pre_times, mags[1:])
    axes[0].set_title('Pre-event times')
    axes[0].set_ylabel('Magnitude')
    axes[0].set_xlabel('Time (Minutes)')
    plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=30)
    axes[1].scatter(pre_times, mags[:-1])
    axes[1].set_title('Post-event times')
    axes[1].set_xlabel('Time (Minutes)')
    plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=30)
    plt.show()


def threeD_seismplot(stations, nodes):
    r"""Function to plot seismicity and stations in a 3D, movable, zoomable \
    space using matplotlibs Axes3D package.

    :type stations: list of tuple
    :param stations: list of one tuple per station of (lat, long, elevation), \
        with up positive.
    :type nodes: list of tuple
    :param nodes: list of one tuple per event of (lat, long, depth) with down \
        positive.
    """
    stalats, stalongs, staelevs = zip(*stations)
    evlats, evlongs, evdepths = zip(*nodes)
    evdepths = [-1 * depth for depth in evdepths]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(evlats, evlongs, evdepths, marker="x", c="k")
    ax.scatter(stalats, stalongs, staelevs, marker="v", c="r")
    ax.set_ylabel("Latitude (deg)")
    ax.set_xlabel("Longitude (deg)")
    ax.set_zlabel("Depth(km)")
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    plt.show()
    return


def pretty_template_plot(template, size=(18.5, 10.5), save=False, title=False,
                         background=False, picks=False):
    r"""Function to make a pretty plot of a single template, designed to work \
    better than the default obspy plotting routine for short data lengths.

    :type template: :class: obspy.Stream
    :param template: Template stream to plot
    :type size: tuple
    :param size: tuple of plot size
    :type save: bool
    :param save: if False will plot to screen, if True will save
    :type title: bool
    :param title: String if set will be the plot title
    :type background: :class: obspy.stream
    :param background: Stream to plot the template within.
    :type picks: list of :class: eqcorrscan.utils.Sfile_util.PICK
    :param picks: List of eqcorrscan type picks.
    """
    fig, axes = plt.subplots(len(template), 1, sharex=True, figsize=size)
    if len(template) > 1:
        axes = axes.ravel()
    else:
        return
    if not background:
        mintime = template.sort(['starttime'])[0].stats.starttime
    else:
        mintime = background.sort(['starttime'])[0].stats.starttime
    template.sort(['network', 'station', 'starttime'])
    lengths = []
    for i, tr in enumerate(template):
        delay = tr.stats.starttime - mintime
        y = tr.data
        x = np.linspace(0, len(y) * tr.stats.delta, len(y))
        x += delay
        if background:
            btr = background.select(station=tr.stats.station,
                                    channel=tr.stats.channel)[0]
            bdelay = btr.stats.starttime - mintime
            by = btr.data
            bx = np.linspace(0, len(by) * btr.stats.delta, len(by))
            bx += bdelay
            axes[i].plot(bx, by, 'k', linewidth=1)
            axes[i].plot(x, y, 'r', linewidth=1.1)
            lengths.append(max(bx[-1], x[-1]))
        else:
            axes[i].plot(x, y, 'k', linewidth=1.1)
            lengths.append(x[-1])
        # print(' '.join([tr.stats.station, str(len(x)), str(len(y))]))
        axes[i].set_ylabel('.'.join([tr.stats.station, tr.stats.channel]),
                           rotation=0, horizontalalignment='right')
        axes[i].yaxis.set_ticks([])
        # Plot the picks if they are given
        if picks:
            tr_picks = [pick for pick in picks if
                        pick.station == tr.stats.station and
                        pick.channel[0] + pick.channel[-1] ==
                        tr.stats.channel[0] + tr.stats.channel[-1]]
            for pick in tr_picks:
                if pick.phase == 'P':
                    pcolor = 'red'
                elif pick.phase == 'S':
                    pcolor = 'blue'
                else:
                    pcolor = 'k'
                pdelay = pick.time - mintime
                # print(pdelay)
                axes[i].axvline(x=pdelay, color=pcolor, linewidth=2)
                # axes[i].plot([pdelay, pdelay], [])
    axes[i].set_xlim([0, max(lengths)])
    axes[len(template)-1].set_xlabel('Time (s) from start of template')
    plt.subplots_adjust(hspace=0, left=0.175, right=0.95, bottom=0.07)
    if title:
        axes[0].set_title(title)
    else:
        plt.subplots_adjust(top=0.98)
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(save)


def NR_plot(stream, NR_stream, detections, false_detections=False,
            size=(18.5, 10), save=False, title=False):
    r"""Function to plot the Network response alongside the streams used -\
    highlights detection times in the network response.

    :type stream: :class: obspy.Stream
    :param stream: Stream to plot
    :type NR_stream: :class: obspy.Stream
    :param NR_stream: Stream for the network response
    :type detections: list of datetime objects
    :param detections: List of the detections
    :type false_detections: list of datetime
    :param false_detections: Either False (default) or list of false detection\
     times
    :type size: tuple
    :param size: Size of figure, default is (18.5,10)
    :type save: bool
    :param save: Save figure or plot to screen, if not False, must be string\
        of save path.
    :type title: str
    :param title: String for the title of the plot, set to False
    """
    import datetime as dt
    import matplotlib.dates as mdates
    fig, axes = plt.subplots(len(stream)+1, 1, sharex=True, figsize=size)
    if len(stream) > 1:
        axes = axes.ravel()
    else:
        return
    mintime = stream.sort(['starttime'])[0].stats.starttime
    stream.sort(['network', 'station', 'starttime'])
    for i, tr in enumerate(stream):
        delay = tr.stats.starttime - mintime
        delay *= tr.stats.sampling_rate
        y = tr.data
        x = [tr.stats.starttime + dt.timedelta(seconds=s /
                                               tr.stats.sampling_rate)
             for s in xrange(len(y))]
        x = mdates.date2num(x)
        axes[i].plot(x, y, 'k', linewidth=1.1)
        axes[i].set_ylabel(tr.stats.station+'.'+tr.stats.channel, rotation=0)
        axes[i].yaxis.set_ticks([])
        axes[i].set_xlim(x[0], x[-1])
    # Plot the network response
    tr = NR_stream[0]
    delay = tr.stats.starttime - mintime
    delay *= tr.stats.sampling_rate
    y = tr.data
    x = [tr.stats.starttime + dt.timedelta(seconds=s / tr.stats.sampling_rate)
         for s in range(len(y))]
    x = mdates.date2num(x)
    axes[i].plot(x, y, 'k', linewidth=1.1)
    axes[i].set_ylabel(tr.stats.station+'.'+tr.stats.channel, rotation=0)
    axes[i].yaxis.set_ticks([])
    axes[-1].set_xlabel('Time')
    axes[-1].set_xlim(x[0], x[-1])
    # Plot the detections!
    ymin, ymax = axes[-1].get_ylim()
    if false_detections:
        for detection in false_detections:
            xd = mdates.date2num(detection)
            axes[-1].plot((xd, xd), (ymin, ymax), 'k--', linewidth=0.5,
                          alpha=0.5)
    for detection in detections:
        xd = mdates.date2num(detection)
        axes[-1].plot((xd, xd), (ymin, ymax), 'r--', linewidth=0.75)
    # Set formatters for x-labels
    mins = mdates.MinuteLocator()
    timedif = tr.stats.endtime.datetime - tr.stats.starttime.datetime
    if timedif.total_seconds() >= 10800 and timedif.total_seconds() <= 25200:
        hours = mdates.MinuteLocator(byminute=[0, 15, 30, 45])
    elif timedif.total_seconds() <= 1200:
        hours = mdates.MinuteLocator(byminute=range(0, 60, 2))
    elif timedif.total_seconds > 25200 and timedif.total_seconds() <= 172800:
        hours = mdates.HourLocator(byhour=range(0, 24, 3))
    elif timedif.total_seconds() > 172800:
        hours = mdates.DayLocator()
    else:
        hours = mdates.MinuteLocator(byminute=range(0, 60, 5))
    hrFMT = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
    axes[-1].xaxis.set_major_locator(hours)
    axes[-1].xaxis.set_major_formatter(hrFMT)
    axes[-1].xaxis.set_minor_locator(mins)
    plt.gcf().autofmt_xdate()
    axes[-1].fmt_xdata = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
    plt.subplots_adjust(hspace=0)
    if title:
        axes[0].set_title(title)
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(save)
    return


def SVD_plot(SVStreams, SValues, stachans, title=False):
    r"""Function to plot the singular vectors from the clustering routines, one\
    plot for each stachan

    :type SVStreams: list of :class:Obspy.Stream
    :param SVStreams: See clustering.SVD_2_Stream - will assume these are\
            ordered by power, e.g. first singular vector in the first stream
    :type SValues: list of float
    :param SValues: List of the singular values corresponding to the SVStreams
    :type stachans: list
    :param stachans: List of station.channel
    """
    for stachan in stachans:
        print(stachan)
        plot_traces = [SVStream.select(station=stachan.split('.')[0],
                                       channel=stachan.split('.')[1])[0]
                       for SVStream in SVStreams]
        fig, axes = plt.subplots(len(plot_traces), 1, sharex=True)
        axes = axes.ravel()
        for i, tr in enumerate(plot_traces):
            y = tr.data
            x = np.linspace(0, len(y) * tr.stats.delta, len(y))
            axes[i].plot(x, y, 'k', linewidth=1.1)
            ylab = 'SV '+str(i+1)+'='+str(round(SValues[i] / len(SValues), 2))
            axes[i].set_ylabel(ylab, rotation=0)
            axes[i].yaxis.set_ticks([])
            print(i)
        axes[-1].set_xlabel('Time (s)')
        plt.subplots_adjust(hspace=0)
        if title:
            axes[0].set_title(title)
        else:
            axes[0].set_title(stachan)
        plt.show()
    return


def plot_synth_real(real_template, synthetic, channels=False):
    r"""Plot multiple channels of data for real data and synthetic.

    :type real_template: obspy.Stream
    :param real_template: Stream of the real template
    :type synthetic: obspy.Stream
    :param synthetic: Stream of synthetic template
    :type channels: list of str
    :param channels: List of tuples of (station, channel) to plot, default is\
            False, which plots all.
    """
    from obspy.signal.cross_correlation import xcorr
    from obspy import Stream
    colours = ['k', 'r']
    labels = ['Real', 'Synthetic']
    if channels:
        real = []
        synth = []
        for stachan in channels:
            real.append(real_template.select(station=stachan[0],
                                             channel=stachan[1]))
            synth.append(synthetic.select(station=stachan[0],
                                          channel=stachan[1]))
        real_template = Stream(real)
        synthetic = Stream(synth)

    # Extract the station and channels
    stachans = list(set([(tr.stats.station, tr.stats.channel)
                         for tr in real_template]))
    fig, axes = plt.subplots(len(stachans), 1, sharex=True, figsize=(5, 10))
    axes = axes.ravel()
    for i, stachan in enumerate(stachans):
        real_tr = real_template.select(station=stachan[0],
                                       channel=stachan[1])[0]
        synth_tr = synthetic.select(station=stachan[0],
                                    channel=stachan[1])[0]
        shift, corr = xcorr(real_tr, synth_tr, 2)
        print 'Shifting by: '+str(shift)+' samples'
        if corr < 0:
            synth_tr.data = synth_tr.data * -1
            corr = corr * -1
        if shift < 0:
            synth_tr.data = synth_tr.data[abs(shift):]
            real_tr.data = real_tr.data[0:len(synth_tr.data)]
        elif shift > 0:
            real_tr.data = real_tr.data[abs(shift):]
            synth_tr.data = synth_tr.data[0:len(real_tr.data)]
        for j, tr in enumerate([real_tr, synth_tr]):
            y = tr.data
            y = y / float(max(abs(y)))
            x = np.linspace(0, len(y) * tr.stats.delta, len(y))
            axes[i].plot(x, y, colours[j], linewidth=2.0, label=labels[j])
            axes[i].get_yaxis().set_ticks([])
        ylab = stachan[0]+'.'+stachan[1]+' cc='+str(round(corr, 2))
        axes[i].set_ylabel(ylab, rotation=0)
    plt.subplots_adjust(hspace=0)
    # axes[0].legend()
    axes[-1].set_xlabel('Time (s)')
    plt.show()


def freq_mag(magnitudes, completeness, max_mag, binsize=0.2):
    r"""Function to make a frequency-magnitude histogram and cumulative \
    density plot.  This can compute a b-value, but not a completeness at \
    the moment.  B-value is computed by linear fitting to section of curve \
    between completeness and max_mag.

    :type magnitudes: list
    :param magnitudes: list of float of magnitudes
    :type completeness: float
    :param completeness: Level to compute the b-value above
    :type max_mag: float
    :param max_mag: Maximum magnitude to try and fit a b-value to
    :type binsize: float
    :param binsize: Width of histogram bins, defaults to 0.2
    """
    # Ensure magnitudes are sorted
    magnitudes.sort()
    fig, ax1 = plt.subplots()
    # Set up the bins, the bin-size could be a variables
    bins = np.arange(min(magnitudes), max(magnitudes), binsize)
    n, bins, patches = ax1.hist(magnitudes, bins, facecolor='Black',
                                alpha=0.5, label='Magnitudes')
    ax1.set_ylabel('Frequency')
    ax1.set_ylim([0, max(n) + 0.5 * max(n)])
    plt.xlabel('Magnitude')
    # Now make the cumulative density function
    cdf = np.arange(len(magnitudes)) / float(len(magnitudes))
    cdf = ((cdf * -1.0) + 1.0) * len(magnitudes)
    ax2 = ax1.twinx()
    ax2.scatter(magnitudes, np.log10(cdf), c='k', marker='+', s=20, lw=2,
                label='Magnitude cumulative density')
    # Now we want to calculate the b-value and plot the fit
    x = []
    y = []
    for i, magnitude in enumerate(magnitudes):
        if magnitude >= completeness <= max_mag:
            x.append(magnitude)
            y.append(cdf[i])
    fit = np.polyfit(x, np.log10(y), 1)
    fit_fn = np.poly1d(fit)
    ax2.plot(magnitudes, fit_fn(magnitudes), '--k',
             label='GR trend, b-value = ' + str(abs(fit[0]))[0:4] +
             '\n $M_C$ = ' + str(completeness))
    ax2.set_ylabel('$Log_{10}$ of cumulative density')
    plt.xlim([min(magnitudes) - 0.5, max(np.log10(cdf)) + 0.2])
    plt.ylim([min(magnitudes) - 0.5, max(np.log10(cdf)) + 1.0])
    plt.legend(loc=2)
    plt.show()


def spec_trace(traces, cmap=None, wlen=0.4, log=False, trc='k',
               tralpha=0.9, size=(10, 13), Fig=None, title=None, show=True):
    r"""Wrapper for _spec_trace, take a stream or list of traces and plots \
    the trace with the spectra beneath it - this just does the overseeing to \
    work out if it needs to add subplots or not.

    :type traces: either stream or list of traces
    :param traces: Traces to be plotted, can be a single obspy.Stream, or a \
        list of obspy.Trace
    :type cmap: str
    :param cmap: [Matplotlib colormap](http://matplotlib.org/examples/color/ \
        colormaps_reference.html)
    :type wlen: float
    :param wlen: Window length for fft in seconds
    :type log: bool
    :param log: Use a log frequency scale
    :type trc: str
    :param trc: Color for the trace.
    :type tralpha: float
    :param tralpha: Opacity level for the seismogram, from transparent (0.0) \
        to opaque (1.0).
    :type size: tuple
    :param size: Plot size, tuple of floats, inches
    :type Fig: matplotlib Fig
    :param axes: Figure to plot onto, defaults to self generating.
    :type show: bool
    :param show: To show plot or not, if false, will return Fig.
    """
    from obspy import Stream
    if isinstance(traces, Stream):
        traces.sort(['station', 'channel'])
    if not Fig:
        Fig = plt.figure(figsize=size)
    for i, tr in enumerate(traces):
        if i == 0:
            ax = Fig.add_subplot(len(traces), 1, i+1)
        else:
            ax = Fig.add_subplot(len(traces), 1, i+1, sharex=ax)
        ax1, ax2 = _spec_trace(tr, wlen=wlen, log=log, trc=trc,
                               tralpha=tralpha, axes=ax)
        ax2.set_yticks([])
        if i < len(traces) - 1:
            plt.setp(ax1.get_xticklabels(), visible=False)
        if type(traces) == list:
            ax2.text(0.005, 0.85, tr.stats.starttime.datetime.
                     strftime('%Y/%m/%d %H:%M:%S'),
                     bbox=dict(facecolor='white', alpha=0.8),
                     transform=ax2.transAxes)
        else:
            ax2.text(0.005, 0.85, '.'.join([tr.stats.station,
                                            tr.stats.channel]),
                     bbox=dict(facecolor='white', alpha=0.8),
                     transform=ax2.transAxes)
        ax2.text(0.005, 0.02, str(np.max(tr.data).round(1)),
                 bbox=dict(facecolor='white', alpha=0.95),
                 transform=ax2.transAxes)
    ax1.set_xlabel('Time (s)')
    Fig.subplots_adjust(hspace=0)
    Fig.text(0.04, 0.5, 'Frequency (Hz)', va='center', rotation='vertical')
    if title:
        plt.suptitle(title)
    if show:
        plt.show()
    else:
        return Fig


def _spec_trace(trace, cmap=None, wlen=0.4, log=False, trc='k',
                tralpha=0.9, size=(10, 2.5), axes=None, title=None):
    r"""Function to plot a trace over that traces spectrogram.
    Uses obspys spectrogram routine.

    :type trace: obspy.Trace
    :param trace: trace to plot
    :type cmap: str
    :param cmap: [Matplotlib colormap](http://matplotlib.org/examples/color/
        colormaps_reference.html)
    :type wlen: float
    :param wlen: Window length for fft in seconds
    :type log: bool
    :param log: Use a log frequency scale
    :type trc: str
    :param trc: Color for the trace.
    :type tralpha: float
    :param tralpha: Opacity level for the seismogram, from transparent (0.0) \
        to opaque (1.0).
    :type size: tuple
    :param size: Plot size, tuple of floats, inches
    :type axes: matplotlib axes
    :param axes: Axes to plot onto, defaults to self generating.
    :type title: str
    :param title: Title for the plot.
    """
    if not axes:
        Fig = plt.figure(figsize=size)
        ax1 = Fig.add_subplot(111)
    else:
        ax1 = axes
    trace.spectrogram(wlen=wlen, log=log, show=False, cmap=cmap, axes=ax1)
    Fig = plt.gcf()
    ax2 = ax1.twinx()
    y = trace.data
    x = np.linspace(0, len(y) / trace.stats.sampling_rate, len(y))
    ax2.plot(x, y, color=trc, linewidth=2.0, alpha=tralpha)
    ax2.set_xlim(min(x), max(x))
    ax2.set_ylim(min(y) * 2, max(y) * 2)
    if title:
        ax1.set_title(' '.join([trace.stats.station, trace.stats.channel,
                                trace.stats.starttime.datetime.
                                strftime('%Y/%m/%d %H:%M:%S')]))
    if not axes:
        Fig.set_size_inches(size)
        Fig.show()
    else:
        return ax1, ax2


if __name__ == "__main__":
    import doctest
    doctest.testmod()
