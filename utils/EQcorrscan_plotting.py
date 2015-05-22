#!/usr/bin/python
"""
Utility code for most of the plots used as part of the EQcorrscan package.
"""
import numpy as np
import matplotlib.pylab as plt
def triple_plot(cccsum, trace, threshold, save=False, savefile=''):
    """
    Main function to make a triple plot with a day-long seismogram, day-long
    correlation sum trace and histogram of the correlation sum to show normality

    :type cccsum: numpy.array
    :type trace: obspy.Trace
    :type threshold: float
    :type save: Bool, optional
    :type savefile: String, optional
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
    ax1.set_ylim([-100,100])
    # Generate the subplot for the correlation sum data
    ax2 = plt.subplot2grid((2,5), (1,0), colspan=4, sharex=ax1)
    # Plot the threshold values
    ax2.plot([min(t), max(t)], [threshold, threshold], color='r', lw=1, label="Threshold")
    ax2.plot([min(t), max(t)], [-threshold,-threshold], color='r', lw=1)
    ax2.plot(t, cccsum, 'k')
    ax2.axis('tight')
    ax2.set_ylim([-5, 5])
    ax2.set_xlabel("Time after %s [hr]" % trace.stats.starttime.isoformat())
    # ax2.legend()
    # Generate a small subplot for the histogram of the cccsum data
    ax3 = plt.subplot2grid((2,5), (1,4), sharey=ax2)
    ax3.hist(cccsum, 200, normed=1, histtype='stepfilled', \
             orientation='horizontal', color='black')
    ax3.set_ylim([-5, 5])
    fig=plt.gcf()
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if not save:
        plt.show()
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
    :type starttime: obspy.UTCDateTime
    :type samp_rate: float
    :type save: Boolean, optional
    :type peaks: List of Tuple, optional
    :type savefile: String, optional
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
    else:
        plt.savefig(savefile)
    return

def cumulative_detections(dates, template_names, save=False, savefile=''):
    """
    Simple plotting function to take a list of UTCDateTime objects and plot
    a cumulative detections list.  Can take dates as a list of lists and will
    plot each list seperately, e.g. if you have dates from more than one
    template it will overlay them in different colours.

    :type dates: list of lists of datetime.datetime
    :param dates: Must be a list of lists of datetime.datetime objects
    :type template_names: list of strings
    :param template_names: List of the template names in order of the dates
    :type save: Boolean, optional
    :type savefile: String, optional
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
        filename,=plt.plot(template_dates,counts, linestyles[j], \
                           color=colors[i], label=template_names[k],\
                           linewidth=3.0)
        k+=1
        plothandles.append(filename)
        if i < len(colors):
            i+=1
        else:
            i=0
            j+=1
    plt.xlabel('Date')
    plt.ylabel('Cumulative detections')
    plt.title('Cumulative detections for all templates')
    plt.legend()#handles=plothandles)
    if save:
        plt.savefig(savefile)
    else:
        plt.show()
    return

def threeD_gridplot(nodes):
    """
    Function to plot in 3D a series of grid points.

    :type nodes: List of tuples
    :param nodes: List of tuples of the form (lat, long, depth)
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
    plt.show()
    return

def detection_timeseries(stream, detector, detections):
    """
    Function to plot the data and detector with detections labelled in red,
    will downsample if too many data points.

    :type stream: obspy.Stream
    :type detector: np.array
    :type detections: np.array
    :param detections: array of positions of detections in samples
    """
    from obspy import Trace
    fig, axes = plt.subplots((len(stream)+1), 1, sharex=True)
    axes=axes.ravel()
    samp_rate=stream[0].stats.sampling_rate
    for i in xrange(len(stream)):
        tr=stream[i]
        if tr.stats.sampling_rate > 10:
            tr.decimate(int(tr.stats.sampling_rate/10))
        time=np.arange(0, tr.stats.npts)/float(tr.stats.delta)
        axes[i].plot(tr.data, time)
        axes[i].set_ylabel(tr.stats.station+'.'+tr.stats.channel)
    detector=Trace(detector)
    detector.stats.sampling_rate=samp_rate
    if detector.stats.sampling_rate > 10:
        detector.decimate(int(detector.stats.sampling_rate/10))
    time=np.arange(0, detector.stats.npts)/float(detector.stats.delta)
    detector=detector.data
    axes[len(axes)].plot(detector, time)
    axes[len(axes)].set_xlabel('Time')
    axes[len(axes)].set_ylabel('Detector')
    return
