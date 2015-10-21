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
    """
    Function to downsample data for plotting by computing the maximum of
    data within chunks, usefil for plotting waveforms or cccsums, large datasets
    that ould otherwise exceed the complexity allowed, and overflow.

    :type tr: :class: obspy.Trace
    :param tr: Trace to be chunked
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type state: str
    :param state: Either 'Min', 'Max', 'Mean' or 'Maxabs' to return one of\
            these for the chunks. Maxabs will return the largest (positive or\
            negative) for that chunk.

    :returns: :class: obspy.Trace
    """
    trout=tr.copy() # Don't do it inplace on data
    x = np.arange(len(tr.data))
    y = tr.data

    chunksize=int(round(tr.stats.sampling_rate/samp_rate))
    # Wrap the array into a 2D array of chunks, truncating the last chunk if
    # chunksize isn't an even divisor of the total size.
    # (This part won't use _any_ additional memory)
    numchunks = y.size // chunksize
    ychunks = y[:chunksize*numchunks].reshape((-1, chunksize))
    xchunks = x[:chunksize*numchunks].reshape((-1, chunksize))

    # Calculate the max, min, and means of chunksize-element chunks...
    if state=='Max':
        trout.data = ychunks.max(axis=1)
    elif state=='Min':
        trout.data = ychunks.min(axis=1)
    elif state=='Mean':
        trout.data = ychunks.mean(axis=1)
    elif state=='Maxabs':
        max_env=ychunks.max(axis=1)
        min_env=ychunks.min(axis=1)
        indeces=np.argmax(np.vstack([np.abs(max_env), np.abs(min_env)]),axis=0)
        stack=np.vstack([max_env, min_env]).T
        trout.data=np.array([stack[i][indeces[i]] for i in xrange(len(stack))])
    xcenters = xchunks.mean(axis=1)
    trout.stats.starttime=tr.stats.starttime+xcenters[0]/tr.stats.sampling_rate
    trout.stats.sampling_rate=samp_rate
    return trout


def triple_plot(cccsum, cccsum_hist, trace, threshold, save=False, savefile=''):
    """
    Main function to make a triple plot with a day-long seismogram, day-long
    correlation sum trace and histogram of the correlation sum to show normality

    :type cccsum: numpy.ndarray
    :param cccsum: Array of the cross-channel cross-correlation sum
    :type cccsum_hist: numpy.ndarray
    :param ccsum_hist: cccsum for histogram plotting, can be the same as cccsum\
                    but included if cccsum is just an envelope.
    :type trace: obspy.Trace
    :param trace: A sample trace from the same time as cccsum
    :type threshold: float
    :param threshold: Detection threshold within cccsum
    :type save: Bool, optional
    :param save: If True will svae and not plot to screen, vice-versa if False
    :type savefile: String, optional
    :param savefile: Path to save figure to, only required if save=True
    """
    if len(cccsum) != len(trace.data):
        print 'cccsum is: '+str(len(cccsum))+' trace is: '+str(len(trace.data))
        raise ValueError('cccsum and trace must have the same number of data points')
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / (df*3600)
    # Generate the subplot for the seismic data
    ax1 = plt.subplot2grid((2,5), (0,0), colspan=4)
    ax1.plot(t, trace.data, 'k')
    ax1.axis('tight')
    ax1.set_ylim([-15*np.mean(np.abs(trace.data)),15*np.mean(np.abs(trace.data))])
    # Generate the subplot for the correlation sum data
    ax2 = plt.subplot2grid((2,5), (1,0), colspan=4, sharex=ax1)
    # Plot the threshold values
    ax2.plot([min(t), max(t)], [threshold, threshold], color='r', lw=1, label="Threshold")
    ax2.plot([min(t), max(t)], [-threshold,-threshold], color='r', lw=1)
    ax2.plot(t, cccsum, 'k')
    ax2.axis('tight')
    ax2.set_ylim([-1.7*threshold, 1.7*threshold])
    ax2.set_xlabel("Time after %s [hr]" % trace.stats.starttime.isoformat())
    # ax2.legend()
    # Generate a small subplot for the histogram of the cccsum data
    ax3 = plt.subplot2grid((2,5), (1,4), sharey=ax2)
    ax3.hist(cccsum_hist, 200, normed=1, histtype='stepfilled', \
             orientation='horizontal', color='black')
    ax3.set_ylim([-5, 5])
    fig=plt.gcf()
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(savefile)
    return

def peaks_plot(data, starttime, samp_rate, save=False, peaks=[(0,0)], \
               savefile=''):
    """
    Simple utility code to plot the correlation peaks to check that the peak
    finding routine is running correctly, used in debugging for the EQcorrscan
    module.

    :type data: numpy.array
    :param data: Numpy array of the data within which peaks have been found
    :type starttime: obspy.UTCDateTime
    :param starttime: Start time for the data
    :type samp_rate: float
    :param samp_rate: Sampling rate of data in Hz
    :type save: Boolean, optional
    :param save: Save figure or plot to screen (False)
    :type peaks: List of Tuple, optional
    :param peaks: List of peak locations and amplitudes (loc, amp)
    :type savefile: String, optional
    :param savefile: Path to save to, only used if save=True
    """
    npts=len(data)
    t = np.arange(npts, dtype=np.float32) / (samp_rate*3600)
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.plot(t, data, 'k')
    ax1.scatter(peaks[0][1]/(samp_rate*3600),abs(peaks[0][0]),color='r', label='Peaks')
    for peak in peaks:
        ax1.scatter(peak[1]/(samp_rate*3600),abs(peak[0]),color='r')
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
    """
    Simple plotting function to take a list of datetime objects and plot
    a cumulative detections list.  Can take dates as a list of lists and will
    plot each list seperately, e.g. if you have dates from more than one
    template it will overlay them in different colours.

    :type dates: list of lists of datetime.datetime
    :param dates: Must be a list of lists of datetime.datetime objects
    :type template_names: list of strings
    :param template_names: List of the template names in order of the dates
    :type save: Boolean, optional
    :param save: Save figure or show to screen
    :type savefile: String, optional
    :param savefile: String to save to.
    """
    # Set up a default series of parameters for lines
    colors=['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', \
            'firebrick', 'purple', 'darkgoldenrod', 'gray']
    linestyles=['-','-.', '--', ':']
    # Check that dates is a list of lists
    if type(dates[0]) != list:
        dates=[dates]
    i=0
    j=0
    k=0
    plothandles=[]
    for template_dates in dates:
        template_dates.sort()
        counts=np.arange(0,len(template_dates))
        print str(i)+' '+str(j)+' '+str(k)
        filename,=plt.plot(template_dates,counts, linestyles[j], \
                           color=colors[i], label=template_names[k],\
                           linewidth=3.0)
        k+=1
        plothandles.append(filename)
        if i < len(colors)-1:
            i+=1
        else:
            i=0
            if j < len(linestyles)-1:
                j+=1
            else:
                j=0
    plt.xlabel('Date')
    plt.ylabel('Cumulative detections')
    plt.title('Cumulative detections for all templates')
    plt.legend(loc=2,prop={'size':8},ncol=2)#handles=plothandles)
    if save:
        plt.savefig(savefile)
        plt.close()
    else:
        plt.show()
    return

def threeD_gridplot(nodes, save=False, savefile=''):
    """
    Function to plot in 3D a series of grid points.

    :type nodes: List of tuples
    :param nodes: List of tuples of the form (lat, long, depth)
    :type save: bool
    :param save: if True will save without plotting to screen, if False\
        (default) will plot to screen but not save
    :type savefile: str
    :param savefile: required if save=True, path to save figure to.
    """
    lats=[]
    longs=[]
    depths=[]
    for node in nodes:
        lats.append(float(node[0]))
        longs.append(float(node[1]))
        depths.append(float(node[2]))
    from mpl_toolkits.mplot3d import Axes3D
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

def multi_event_singlechan(streams, picks, clip=10.0, pre_pick=2.0,\
                           freqmin=False, freqmax=False, realign=False, \
                           cut=(-3.0,5.0), PWS=False, title=False):
    """
    Function to plot data from a single channel at a single station for multiple
    events - data will be alligned by their pick-time given in the picks

    :type streams: List of :class:obspy.stream
    :param streams: List of the streams to use, can contain more traces than\
        you plan on plotting
    :type picks: List of :class:PICK
    :param picks: List of picks, one for each stream
    :type clip: float
    :param clip: Length in seconds to plot, defaults to 10.0
    :type pre_pick: Float
    :param pre_pick: Length in seconds to extract and plot before the pick,\
        defaults to 2.0
    :type freqmin: float
    :param freqmin: Low cut for bandpass in Hz
    :type freqmax: float
    :param freqmax: High cut for bandpass in Hz
    :type realign: Bool
    :param realign: To compute best alignement based on correlation or not.
    :type cut: tuple:
    :param cut: tuple of start and end times for cut in seconds from the pick
    :type PWS: bool
    :param PWS: compute Phase Weighted Stack, if False, will compute linear stack
    :type title: str
    :param title: Plot title.

    :returns: Alligned and cut traces, and new picks
    """
    import stacking, copy
    from eqcorrscan.core.match_filter import normxcorr2
    from obspy import Stream
    fig, axes = plt.subplots(len(picks)+1, 1, sharex=True, figsize=(7, 12))
    axes = axes.ravel()
    traces=[]
    al_traces=[]
    # Keep input safe
    plist=copy.deepcopy(picks)
    st_list=copy.deepcopy(streams)
    for i in xrange(len(plist)):
        if st_list[i].select(station=plist[i].station, \
            channel='*'+plist[i].channel[-1]):
            tr=st_list[i].select(station=plist[i].station, \
                channel='*'+plist[i].channel[-1])[0]
        else:
            print 'No data for '+plist[i].station+'.'+plist[i].channel
            continue
        tr.detrend('linear')
        if freqmin:
            tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        if realign:
            tr_cut=tr.copy()
            tr_cut.trim(plist[i].time+cut[0], plist[i].time+cut[1],\
                        nearest_sample=False)
            if len(tr_cut.data)<=0.5*(cut[1]-cut[0])*tr_cut.stats.sampling_rate:
                print 'Not enough in the trace for '+plist[i].station+'.'+plist[i].channel
                print 'Suggest removing pick from sfile at time '+str(plist[i].time)
            else:
                al_traces.append(tr_cut)
        else:
            tr.trim(plist[i].time-pre_pick, plist[i].time+clip-pre_pick,\
                    nearest_sample=False)
        if len(tr.data)==0:
            print 'No data in the trace for '+plist[i].station+'.'+plist[i].channel
            print 'Suggest removing pick from sfile at time '+str(plist[i].time)
            continue
        traces.append(tr)
    if realign:
        shift_len=int(0.25*(cut[1]-cut[0])*al_traces[0].stats.sampling_rate)
        shifts=stacking.align_traces(al_traces, shift_len)
        for i in xrange(len(shifts)):
            print 'Shifting by '+str(shifts[i])+' seconds'
            plist[i].time-=shifts[i]
            traces[i].trim(plist[i].time-pre_pick, plist[i].time+clip-pre_pick,\
                           nearest_sample=False)
    # We now have a list of traces
    traces=[(trace, trace.stats.starttime.datetime) for trace in traces]
    traces.sort(key=lambda tup:tup[1])
    traces=[trace[0] for trace in traces]
    # Plot the traces
    for i in xrange(len(traces)):
        tr=traces[i]
        y = tr.data
        x = np.arange(len(y))
        x = x/tr.stats.sampling_rate # convert to seconds
        axes[i+1].plot(x, y, 'k', linewidth=1.1)
        # axes[i+1].set_ylabel(tr.stats.starttime.datetime.strftime('%Y/%m/%d %H:%M'),\
                             # rotation=0)
        axes[i+1].yaxis.set_ticks([])
        i+=1
    traces=[Stream(trace) for trace in traces]
    if PWS:
        linstack=stacking.PWS_stack(traces)
    else:
        linstack=stacking.linstack(traces)
    tr=linstack.select(station=picks[0].station, \
            channel='*'+picks[0].channel[-1])[0]
    y = tr.data
    x = np.arange(len(y))
    x = x/tr.stats.sampling_rate
    axes[0].plot(x, y, 'r', linewidth=2.0)
    axes[0].set_ylabel('Stack', rotation=0)
    axes[0].yaxis.set_ticks([])
    for i in xrange(len(traces)):
        cc=normxcorr2(tr.data, traces[i][0].data)
        axes[i+1].set_ylabel('cc='+str(round(np.max(cc),2)), rotation=0)
        axes[i+1].text(0.9, 0.15, str(round(np.max(traces[i][0].data))), \
                       bbox=dict(facecolor='white', alpha=0.95),\
                       transform=axes[i+1].transAxes)
        axes[i+1].text(0.7, 0.85, traces[i][0].stats.starttime.datetime.strftime('%Y/%m/%d %H:%M:%S'), \
                       bbox=dict(facecolor='white', alpha=0.95),\
                       transform=axes[i+1].transAxes)
    axes[-1].set_xlabel('Time (s)')
    if title:
        axes[0].set_title(title)
    plt.subplots_adjust(hspace=0)
    plt.show()
    return traces, plist

def detection_multiplot(stream, template, times, streamcolour='k',\
        templatecolour='r'):
    """
    Function to plot the stream of data that has been detected in, with the
    template on top of it timed according to a list of given times, just a
    pretty way to show a detection!

    :type stream: obspy.Stream
    :param stream: Stream of data to be plotted as the base (black)
    :type template: obspy.Stream
    :param template: Template to be plotted on top of the base stream (red)
    :type times: List of datetime.datetime
    :param times: list of times of detections in the order of the channels in
                template.
    :type streamcolour: str
    :param streamcolour: String of matplotlib colour types for the stream
    :type templatecolour: str
    :param templatecolour: Colour to plot the template in.
    """
    import datetime as dt
    fig, axes = plt.subplots(len(template), 1, sharex=True)
    axes = axes.ravel()
    print 'Template has '+str(len(template))+' channels'
    for i in xrange(len(template)):
        template_tr=template[i]
        print 'Working on: '+template_tr.stats.station+' '+\
                template_tr.stats.channel
        image=stream.select(station=template_tr.stats.station,\
                            channel='*'+template_tr.stats.channel[-1])
        if not image:
            print 'No data for '+template_tr.stats.station+' '+\
                    template_tr.stats.channel
            continue
        image=image.merge()[0]
        t_start=(times[i]-image.stats.starttime.datetime) # Gives a timedelta
        print t_start
        # Downsample if needed
        if image.stats.sampling_rate > 20:
            image.decimate(int(image.stats.sampling_rate/20))
        if template_tr.stats.sampling_rate > 20:
            template_tr.decimate(int(template_tr.stats.sampling_rate/20))
        image_times=[image.stats.starttime.datetime+dt.timedelta((j*image.stats.delta)/86400)\
                for j in xrange(len(image.data))] # Give list of datetime objects
        t_start=t_start+image.stats.starttime.datetime
        template_times=[t_start+dt.timedelta((j*template_tr.stats.delta)/86400)\
                for j in xrange(len(template_tr.data))]
        print image_times[0]
        print template_times[0]
        axes[i].plot(image_times, image.data,'k')
        axes[i].plot(template_times, template_tr.data,'r')
        axes[i].set_ylabel(template_tr.stats.station+'.'+template_tr.stats.channel)
    axes[len(axes)-1].set_xlabel('Time')
    plt.show()
    return

def interev_mag_sfiles(sfiles):
    """
    Function to plot interevent-time versus magnitude for series of events.

    :type sfiles: List
    :param sfiles: List of sfiles to read from
    """
    import Sfile_util
    times=[Sfile_util.readheader(sfile).time for sfile in sfiles]
    mags=[Sfile_util.readheader(sfile).Mag_1 for sfile in sfiles]
    interev_mag(times, mags)

def interev_mag(times, mags):
    """
    Function to plot interevent times against magnitude for given times
    and magnitudes.

    :type times: list of datetime
    :param times: list of the detection times, must be sorted the same as mags
    :type mags: list of float
    :param mags: list of magnitudes
    """
    l = [(times[i], mags[i]) for i in xrange(len(times))]
    l.sort(key=lambda tup:tup[0])
    times=[x[0] for x in l]
    mags=[x[1] for x in l]
    # Make two subplots next to each other of time before and time after
    fig, axes = plt.subplots(1,2, sharey=True)
    axes = axes.ravel()
    pre_times=[]
    post_times=[]
    for i in xrange(len(times)):
        if i > 0:
            pre_times.append((times[i]-times[i-1])/60)
        if i < len(times)-1:
            post_times.append((times[i+1]-times[i])/60)
    axes[0].scatter(pre_times, mags[1:])
    axes[0].set_title('Pre-event times')
    axes[0].set_ylabel('Magnitude')
    axes[0].set_xlabel('Time (Minutes)')
    # axes[0].set_xlim([0, max(pre_times)+(0.1*(max(pre_times)-min(pre_times)))])
    plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=30 )
    axes[1].scatter(pre_times, mags[:-1])
    axes[1].set_title('Post-event times')
    axes[1].set_xlabel('Time (Minutes)')
    # axes[1].set_xlim([0, max(post_times)+(0.1*(max(post_times)-min(post_times)))])
    plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=30 )
    plt.show()

def threeD_seismplot(stations, nodes):
    """
    Function to plot seismicity and stations in a 3D, movable, zoomable space
    using matplotlibs Axes3D package.

    :type stations: list of tuple
    :param stations: list of one tuple per station of (lat, long, elevation),
                    with up positive
    :type nodes: list of tuple
    :param nodes: list of one tuple per event of (lat, long, depth) with down
                positive
    """
    stalats=[]
    stalongs=[]
    staelevs=[]
    evlats=[]
    evlongs=[]
    evdepths=[]
    for station in stations:
        stalats+=[station[0]]
        stalongs+=[station[1]]
        staelevs+=[station[2]]
    for node in nodes:
        evlats+=[node[0]]
        evlongs+=[node[1]]
        evdepths+=[-1*node[2]]
    from mpl_toolkits.mplot3d import Axes3D
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

def Noise_plotting(station, channel, PAZ, datasource):
    """
    Function to make use of obspy's PPSD functionality to read in data from
    a single station and the poles-and-zeros for that station before plotting
    the PPSD for this station.  See McNamara(2004) for more details.

    :type station: String
    :param station: Station name as it is in the filenames in the database
    :type channel: String
    :param channel: Channel name as it is in the filenames in the database
    :type PAZ: Dict
    :param PAZ: Must contain, Poles, Zeros, Sensitivity, Gain
        :type Poles: List of Complex
        :type Zeros: List of Complex
        :type Sensitivity: Float
        :type Gain: Float
    :type datasource: String
    :param datasource: The directory in which data can be found, can contain
                        wildcards.

    :returns: PPSD object
    """
    from obspy.signal import PPSD
    from obspy import read as obsread
    import glob

    stafiles=glob.glob(datasource+'/*'+station+'*'+channel+'*')
    stafiles.sort()
    # Initialize PPSD
    st=obsread(stafiles[0])
    ppsd = PPSD(st[0].stats, PAZ)
    for stafile in stafiles[1:]:
        print 'Adding waveform from: '+stafile
        st=obsread(stafile)
        # Add after read to conserve memory
        ppsd.add(st)
    # Plot the PPSD
    ppsd.plot()
    return ppsd

def pretty_template_plot(template, size=(18.5, 10.5), save=False, title=False,\
                        background=False):
    """
    Function to make a pretty plot of a single template, designed to work better
    than the default obspy plotting routine for short data lengths.

    :type template: :class: obspy.Stream
    :param template: Template stream to plot
    :type size: tuple
    :param size: tuple of plot size
    :type save: Boolean
    :param save: if False will plot to screen, if True will save
    :type title: Boolean
    :param title: String if set will be the plot title
    :type backrgound: :class: obspy.stream
    :param background: Stream to plot the template within.
    """
    fig, axes = plt.subplots(len(template), 1, sharex=True, figsize=size)
    if len(template) > 1:
        axes = axes.ravel()
    else:
        return
    if not background:
        mintime=template.sort(['starttime'])[0].stats.starttime
    else:
        mintime=background.sort(['starttime'])[0].stats.starttime
    i=0
    template.sort(['network','station', 'starttime'])
    for tr in template:
        delay=tr.stats.starttime-mintime
        delay*=tr.stats.sampling_rate
        y=tr.data
        x=np.arange(len(y))
        x+=delay
        x=x/tr.stats.sampling_rate
        # x=np.arange(delay, (len(y)*tr.stats.sampling_rate)+delay,\
            # tr.stats.sampling_rate)
        if background:
            btr=background.select(station=tr.stats.station, \
                                channel=tr.stats.channel)[0]
            bdelay=btr.stats.starttime-mintime
            bdelay*=btr.stats.sampling_rate
            by=btr.data
            bx=np.arange(len(by))
            bx+=bdelay
            bx=bx/btr.stats.sampling_rate
            axes[i].plot(bx,by,'k',linewidth=1)
            axes[i].plot(x, y, 'r', linewidth=1.1)
        else:
            axes[i].plot(x, y, 'k', linewidth=1.1)
        print tr.stats.station+' '+str(len(x))+' '+str(len(y))
        axes[i].set_ylabel(tr.stats.station+'.'+tr.stats.channel, rotation=0)
        axes[i].yaxis.set_ticks([])
        i+=1
    axes[i-1].set_xlabel('Time (s) from start of template')
    plt.subplots_adjust(hspace=0)
    if title:
        axes[0].set_title(title)
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(save)

def NR_plot(stream, NR_stream, detections, false_detections=False,\
            size=(18.5,10), save=False, title=False):
    """
    Function to plot the Network response alongside the streams used - highlights
    detection times in the network response

    :type stream: :class: obspy.Stream
    :param stream: Stream to plot
    :type NR_stream: :class: obspy.Stream
    :param NR_stream: Stream for the network response
    :type detections: List of datetime objects
    :param detections: List of the detections
    :type false_detections: List of datetime
    :param false_detections: Either False (default) or list of false detection\
     times
    :type size: tuple
    :param size: Size of figure, default is (18.5,10)
    :type save: bool
    :param save: Save figure or plot to screen, if not False, must be string of\
        save path
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
    mintime=stream.sort(['starttime'])[0].stats.starttime
    i=0
    stream.sort(['network','station', 'starttime'])
    for tr in stream:
        delay=tr.stats.starttime-mintime
        delay*=tr.stats.sampling_rate
        y=tr.data
        x=[tr.stats.starttime + dt.timedelta(seconds=s/tr.stats.sampling_rate)\
           for s in xrange(len(y))]
        x=mdates.date2num(x)
        axes[i].plot(x, y, 'k', linewidth=1.1)
        axes[i].set_ylabel(tr.stats.station+'.'+tr.stats.channel, rotation=0)
        axes[i].yaxis.set_ticks([])
        axes[i].set_xlim(x[0], x[-1])
        i+=1
    # Plot the network response
    tr=NR_stream[0]
    delay=tr.stats.starttime-mintime
    delay*=tr.stats.sampling_rate
    y=tr.data
    x=[tr.stats.starttime + dt.timedelta(seconds=s/tr.stats.sampling_rate)\
        for s in xrange(len(y))]
    x=mdates.date2num(x)
    axes[i].plot(x, y, 'k', linewidth=1.1)
    axes[i].set_ylabel(tr.stats.station+'.'+tr.stats.channel, rotation=0)
    axes[i].yaxis.set_ticks([])
    axes[-1].set_xlabel('Time')
    axes[-1].set_xlim(x[0], x[-1])
    # Plot the detections!
    ymin, ymax = axes[-1].get_ylim()
    if false_detections:
        for detection in false_detections:
            xd=mdates.date2num(detection)
            axes[-1].plot((xd, xd), (ymin, ymax), 'k--', linewidth=0.5, alpha=0.5)
    for detection in detections:
        xd=mdates.date2num(detection)
        axes[-1].plot((xd, xd), (ymin, ymax), 'r--', linewidth=0.75)
    # Set formatters for x-labels
    mins=mdates.MinuteLocator()
    if (tr.stats.endtime.datetime-tr.stats.starttime.datetime).total_seconds() >= 10800\
       and (tr.stats.endtime.datetime-tr.stats.starttime.datetime).total_seconds() <= 25200:
        hours=mdates.MinuteLocator(byminute=[0,15,30,45])
    elif(tr.stats.endtime.datetime-tr.stats.starttime.datetime).total_seconds() <= 1200:
        hours=mdates.MinuteLocator(byminute=range(0,60,2))
    elif (tr.stats.endtime.datetime-tr.stats.starttime.datetime).total_seconds() > 25200\
        and (tr.stats.endtime.datetime-tr.stats.starttime.datetime).total_seconds() <= 172800:
        hours=mdates.HourLocator(byhour=range(0,24,3))
    elif (tr.stats.endtime.datetime-tr.stats.starttime.datetime).total_seconds() > 172800:
        hours=mdates.DayLocator()
    else:
        hours=mdates.MinuteLocator(byminute=range(0,60,5))
    hrFMT=mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
    axes[-1].xaxis.set_major_locator(hours)
    axes[-1].xaxis.set_major_formatter(hrFMT)
    axes[-1].xaxis.set_minor_locator(mins)
    plt.gcf().autofmt_xdate()
    axes[-1].fmt_xdata=mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
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
    """
    Function to plot the singular vectors from the clustering routines, one
    plot for each stachan

    :type SVStreams: List of :class:Obspy.Stream
    :param SVStreams: See clustering.SVD_2_Stream - will assume these are\
            ordered by power, e.g. first singular vector in the first stream
    :type SValues: List of float
    :param SValues: List of the singular values corresponding to the SVStreams
    :type stachans: List
    :param stachans: List of station.channel
    """
    for stachan in stachans:
        print stachan
        plot_traces=[SVStream.select(station=stachan.split('.')[0],\
                                     channel=stachan.split('.')[1])[0]\
                     for SVStream in SVStreams]
        fig, axes = plt.subplots(len(plot_traces), 1, sharex=True)
        axes = axes.ravel()
        i=0
        for tr in plot_traces:
            y = tr.data
            x = np.arange(len(y))
            x = x*tr.stats.delta
            axes[i].plot(x,y,'k', linewidth=1.1)
            axes[i].set_ylabel('SV '+str(i+1)+'='+\
                str(round(SValues[i]/len(SValues),2)), rotation=0)
            axes[i].yaxis.set_ticks([])
            print i
            i+=1
        axes[-1].set_xlabel('Time (s)')
        plt.subplots_adjust(hspace=0)
        if title:
            axes[0].set_title(title)
        else:
            axes[0].set_title(stachan)
        plt.show()
    return
