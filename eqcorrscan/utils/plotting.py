"""
Utility code for most of the plots used as part of the EQcorrscan package.

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

import numpy as np
import warnings
import datetime as dt

import matplotlib.pylab as plt
import matplotlib.dates as mdates
from copy import deepcopy
from collections import Counter
from obspy import UTCDateTime, Stream, Catalog
from obspy.signal.cross_correlation import xcorr

from eqcorrscan.core.match_filter import DETECTION, normxcorr2
from eqcorrscan.utils import stacking, sfile_util


def _check_save_args(save, savefile):
    if save and not savefile:
        raise IOError('save is set to True, but no savefile is given')
    else:
        return


def chunk_data(tr, samp_rate, state='mean'):
    """
    Downsample data for plotting.

    Computes the maximum of data within chunks, useful for plotting waveforms
    or cccsums, large datasets that would otherwise exceed the complexity
    allowed, and overflow.

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to be chunked
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type state: str
    :param state:
        Either 'Min', 'Max', 'Mean' or 'Maxabs' to return one of these for the
        chunks. Maxabs will return the largest (positive or negative) for
        that chunk.

    :returns: :class:`obspy.core.trace.Trace`
    """
    trout = tr.copy()  # Don't do it inplace on data
    x = np.arange(len(tr.data))
    y = tr.data

    chunksize = int(round(tr.stats.sampling_rate / samp_rate))
    # Wrap the array into a 2D array of chunks, truncating the last chunk if
    # chunksize isn't an even divisor of the total size.
    # (This part won't use _any_ additional memory)
    numchunks = int(y.size // chunksize)
    ychunks = y[:chunksize * numchunks].reshape((-1, chunksize))
    xchunks = x[:chunksize * numchunks].reshape((-1, chunksize))

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
                              for i in range(len(stack))])
    xcenters = xchunks.mean(axis=1)
    trout.stats.starttime = tr.stats.starttime + xcenters[0] /\
        tr.stats.sampling_rate
    trout.stats.sampling_rate = samp_rate
    return trout


def xcorr_plot(template, image, shift=None, cc=None, cc_vec=None, save=False,
               savefile=None):
    """
    Plot a template overlying an image aligned by correlation.

    :type template: numpy.ndarray
    :param template: Short template image
    :type image: numpy.ndarray
    :param image: Long master image
    :type shift: int
    :param shift: Shift to apply to template relative to image, in samples
    :type cc: float
    :param cc: Cross-correlation at shift
    :type cc_vec: numpy.ndarray
    :param cc_vec: Cross-correlation vector.
    :type save: bool
    :param save: Whether to save the plot or not.
    :type savefile: str
    :param savefile: File name to save to

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.plotting import xcorr_plot
    >>> from eqcorrscan.utils.stacking import align_traces
    >>> st = read().detrend('simple').filter('bandpass', freqmin=2, freqmax=15)
    >>> shifts, ccs = align_traces([st[0], st[1]], 40)
    >>> shift = shifts[1] * st[1].stats.sampling_rate
    >>> cc = ccs[1]
    >>> xcorr_plot(template=st[1].data, image=st[0].data, shift=shift,
    ...            cc=cc) # doctest: +SKIP

    .. image:: ../../plots/xcorr_plot.png
    """
    _check_save_args(save, savefile)
    if cc is None or shift is None:
        if not isinstance(cc_vec, np.ndarray):
            print('Given cc: %s and shift: %s' % (cc, shift))
            raise IOError('Must provide either cc_vec, or cc and shift')
        shift = np.abs(cc_vec).argmax()
        cc = cc_vec[shift]
    x = np.arange(len(image))
    plt.plot(x, image / abs(image).max(), 'k', lw=1.3, label='Image')
    x = np.arange(len(template)) + shift
    plt.plot(x, template / abs(template).max(), 'r', lw=1.1, label='Template')
    plt.title('Shift=%s, Correlation=%s' % (shift, cc))
    fig = plt.gcf()
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def triple_plot(cccsum, cccsum_hist, trace, threshold, save=False,
                savefile=None):
    """
    Plot a seismogram, correlogram and histogram.

    :type cccsum: numpy.ndarray
    :param cccsum: Array of the cross-channel cross-correlation sum
    :type cccsum_hist: numpy.ndarray
    :param cccsum_hist: cccsum for histogram plotting, can be the same as \
        cccsum but included if cccsum is just an envelope.
    :type trace: obspy.core.trace.Trace
    :param trace: A sample trace from the same time as cccsum
    :type threshold: float
    :param threshold: Detection threshold within cccsum
    :type save: bool
    :param save: If True will save and not plot to screen, vice-versa if False
    :type savefile: str
    :param savefile: Path to save figure to, only required if save=True

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.core.match_filter import normxcorr2
    >>> from eqcorrscan.utils.plotting import triple_plot
    >>> st = read()
    >>> template = st[0].copy().trim(st[0].stats.starttime + 8,
    ...                            st[0].stats.starttime + 12)
    >>> tr = st[0]
    >>> ccc = normxcorr2(template=template.data, image=tr.data)
    >>> tr.data = tr.data[0:len(ccc[0])]
    >>> triple_plot(cccsum=ccc, cccsum_hist=ccc, trace=tr,
    ...             threshold=0.8) # doctest: +SKIP


    .. image:: ../../plots/triple_plot.png
    """
    _check_save_args(save, savefile)
    if len(cccsum) != len(trace.data):
        print('cccsum is: ' +
              str(len(cccsum)) + ' trace is: ' + str(len(trace.data)))
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
        plt.close()
    return fig


def peaks_plot(data, starttime, samp_rate, save=False, peaks=[(0, 0)],
               savefile=None):
    """
    Plot peaks to check that the peak finding routine is running correctly.

    Used in debugging for the EQcorrscan module.

    :type data: numpy.array
    :param data: Numpy array of the data within which peaks have been found
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Start time for the data
    :type samp_rate: float
    :param samp_rate: Sampling rate of data in Hz
    :type save: bool
    :param save: Save figure or plot to screen.
    :type peaks: list
    :param peaks: List of tuples of peak locations and amplitudes (loc, amp)
    :type savefile: str
    :param savefile: Path to save to, only used if save=True

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> import numpy as np
    >>> from eqcorrscan.utils import findpeaks
    >>> from eqcorrscan.utils.plotting import peaks_plot
    >>> from obspy import UTCDateTime
    >>> data = np.random.randn(200)
    >>> data[30]=100
    >>> data[60]=40
    >>> threshold = 10
    >>> peaks = findpeaks.find_peaks2_short(data, threshold, 3)
    >>> peaks_plot(data=data, starttime=UTCDateTime("2008001"),
    ...            samp_rate=10, peaks=peaks)  # doctest: +SKIP

    .. plot::

        import matplotlib.pyplot as plt
        import numpy as np
        from eqcorrscan.utils import findpeaks
        from eqcorrscan.utils.plotting import peaks_plot
        from obspy import UTCDateTime
        data = np.random.randn(200)
        data[30]=100
        data[60]=40
        threshold = 10
        peaks = findpeaks.find_peaks2_short(data, threshold, 3)
        peaks_plot(data=data, starttime=UTCDateTime("2008001"),
                   samp_rate=10, peaks=peaks)
    """
    _check_save_args(save, savefile)
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
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def cumulative_detections(dates=None, template_names=None, detections=None,
                          plot_grouped=False, show=True, plot_legend=True,
                          save=False, savefile=None):
    """
    Plot cumulative detections in time.

    Simple plotting function to take a list of datetime objects and plot
    a cumulative detections list.  Can take dates as a list of lists and will
    plot each list separately, e.g. if you have dates from more than one
    template it will overlay them in different colours.

    :type dates: list
    :param dates: Must be a list of lists of datetime.datetime objects
    :type template_names: list
    :param template_names: List of the template names in order of the dates
    :type detections: list
    :param detections: List of :class:`eqcorrscan.core.match_filter.DETECTION`
    :type plot_grouped: bool
    :param plot_grouped: Plot detections for each template individually, or \
        group them all together - set to False (plot template detections \
        individually) by default.
    :type show: bool
    :param show: Whether or not to show the plot, defaults to True.
    :type plot_legend: bool
    :param plot_legend: Specify whether to plot legend of template names. \
        Defaults to True.
    :type save: bool
    :param save: Save figure or show to screen, optional
    :type savefile: str
    :param savefile: String to save to, required is save=True

    :returns: :class:`matplotlib.figure.Figure`

    .. note::
        Can either take lists of
        :class:`eqcorrscan.core.match_filter.DETECTION` objects directly, or
        two lists of dates and template names - either/or, not both.

    .. rubric:: Example

    >>> import datetime as dt
    >>> import numpy as np
    >>> from eqcorrscan.utils.plotting import cumulative_detections
    >>> dates = []
    >>> for i in range(3):
    ...     dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
    ...                   for n in np.random.randn(100)])
    >>> cumulative_detections(dates, ['a', 'b', 'c'],
    ...                       show=True) # doctest: +SKIP

    .. plot::

        import datetime as dt
        import numpy as np
        from eqcorrscan.utils.plotting import cumulative_detections
        dates = []
        for i in range(3):
            dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
                          for n in np.random.randn(100)])
        cumulative_detections(dates, ['a', 'b', 'c'], show=True)
    """
    _check_save_args(save, savefile)
    # Set up a default series of parameters for lines
    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black',
              'firebrick', 'purple', 'darkgoldenrod', 'gray']
    linestyles = ['-', '-.', '--', ':']
    # Check that dates is a list of lists
    if not detections:
        if type(dates[0]) != list:
            dates = [dates]
    else:
        dates = []
        template_names = []
        for detection in detections:
            if not type(detection) == DETECTION:
                msg = 'detection not of type: ' +\
                    'eqcorrscan.core.match_filter.DETECTION'
                raise IOError(msg)
            dates.append(detection.detect_time.datetime)
            template_names.append(detection.template_name)
        _dates = []
        _template_names = []
        for template_name in list(set(template_names)):
            _template_names.append(template_name)
            _dates.append([date for i, date in enumerate(dates)
                           if template_names[i] == template_name])
        dates = _dates
        template_names = _template_names
    if plot_grouped:
        _dates = []
        for template_dates in dates:
            _dates += template_dates
        dates = [_dates]
        template_names = ['All templates']
    i = 0
    j = 0
    # This is an ugly way of looping through colours and linestyles, it would
    # be better with itertools functions...
    fig, ax1 = plt.subplots()
    min_date = min([min(_d) for _d in dates])
    for k, template_dates in enumerate(dates):
        template_dates.sort()
        plot_dates = deepcopy(template_dates)
        plot_dates.insert(0, min_date)
        counts = np.arange(-1, len(template_dates))
        ax1.step(plot_dates, counts, linestyles[j],
                 color=colors[i], label=template_names[k],
                 linewidth=3.0)
        if i < len(colors) - 1:
            i += 1
        else:
            i = 0
            if j < len(linestyles) - 1:
                j += 1
            else:
                j = 0
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Cumulative detections')
    plt.title('Cumulative detections for all templates')
    # Set formatters for x-labels
    mins = mdates.MinuteLocator()
    max_date = dates[0][0]
    min_date = max_date
    for date_list in dates:
        if max(date_list) > max_date:
            max_date = max(date_list)
        if min(date_list) < min_date:
            min_date = min(date_list)
    timedif = max_date - min_date
    if 10800 <= timedif.total_seconds() <= 25200:
        hours = mdates.MinuteLocator(byminute=[0, 30])
        mins = mdates.MinuteLocator(byminute=range(0, 60, 10))
    elif 7200 <= timedif.total_seconds() < 10800:
        hours = mdates.MinuteLocator(byminute=[0, 15, 30, 45])
        mins = mdates.MinuteLocator(byminute=range(0, 60, 5))
    elif timedif.total_seconds() <= 1200:
        hours = mdates.MinuteLocator(byminute=range(0, 60, 2))
        mins = mdates.MinuteLocator(byminute=range(0, 60, 0.5))
    elif 25200 < timedif.total_seconds() <= 86400:
        hours = mdates.HourLocator(byhour=range(0, 24, 3))
        mins = mdates.HourLocator(byhour=range(0, 24, 1))
    elif 86400 < timedif.total_seconds() <= 172800:
        hours = mdates.HourLocator(byhour=range(0, 24, 6))
        mins = mdates.HourLocator(byhour=range(0, 24, 1))
    elif timedif.total_seconds() > 172800:
        hours = mdates.AutoDateLocator()
        mins = mdates.HourLocator(byhour=range(0, 24, 3))
    else:
        hours = mdates.MinuteLocator(byminute=range(0, 60, 5))
    # Minor locator overruns maxticks for ~year-long datasets
    if timedif.total_seconds() < 172800:
        ax1.xaxis.set_minor_locator(mins)
        hrFMT = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
    else:
        hrFMT = mdates.DateFormatter('%Y/%m/%d')
    ax1.xaxis.set_major_locator(hours)
    ax1.xaxis.set_major_formatter(hrFMT)
    plt.gcf().autofmt_xdate()
    locs, labels = plt.xticks()
    ax1.set_ylim([0, max([len(_d) for _d in dates])])
    plt.setp(labels, rotation=15)
    if plot_legend:
        ax1.legend(loc=2, prop={'size': 8}, ncol=2)
    if save:
        fig.savefig(savefile)
        plt.close()
    else:
        if show:
            plt.show()
    return fig


def threeD_gridplot(nodes, save=False, savefile=None):
    r"""Plot in a series of grid points in 3D.

    :type nodes: list
    :param nodes: List of tuples of the form (lat, long, depth)
    :type save: bool
    :param save: if True will save without plotting to screen, if False \
        (default) will plot to screen but not save
    :type savefile: str
    :param savefile: required if save=True, path to save figure to.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from eqcorrscan.utils.plotting import threeD_gridplot
    >>> nodes = [(-43.5, 170.4, 4), (-43.3, 170.8, 12), (-43.4, 170.3, 8)]
    >>> threeD_gridplot(nodes=nodes)  # doctest: +SKIP

    .. plot::

        from eqcorrscan.utils.plotting import threeD_gridplot
        nodes = [(-43.5, 170.4, 4), (-43.3, 170.8, 12), (-43.4, 170.3, 8)]
        threeD_gridplot(nodes=nodes)
    """
    _check_save_args(save, savefile)
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
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def multi_event_singlechan(streams, catalog, station, channel,
                           clip=10.0, pre_pick=2.0,
                           freqmin=False, freqmax=False, realign=False,
                           cut=(-3.0, 5.0), PWS=False, title=False,
                           save=False, savefile=None):
    """
    Plot data from a single channel for multiple events.

    Data will be aligned by their pick-time given in the appropriate picks.
    Requires an individual stream for each event you want to plot,
    events are stored in the
    :class:`obspy.core.event.Catalog` object, and there must be picks present
    for the streams you wish to plot.  Events will be aligned if
    `realign=True`, in this case the traces will be aligned using the window
    defined by `cut`.

    :type streams: list
    :param streams:
        List of the :class:`obspy.core.stream.Stream` objects to use, can
        contain more traces than you plan on plotting (e.g. from more channels)
         - must be in the same order as events in catalog.
    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog of events, one for each stream.
    :type station: str
    :param station: Station to plot.
    :type channel: str
    :param channel: Channel to plot.
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
    :param realign:
        To compute best alignment based on correlation with the stack or not.
    :type cut: tuple
    :param cut: tuple of start and end times for cut in seconds from the \
        pick, used for alignment.  Will only use this window to align the \
        traces.
    :type PWS: bool
    :param PWS: compute Phase Weighted Stack, if False, will compute linear \
        stack for alignment.
    :type title: str
    :param title: Plot title.
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True

    :returns: Aligned and cut :class:`obspy.core.trace.Trace`
    :rtype: list
    :returns:
        New picks in based on alignment (if alignment is performed, if not
        will return the same as input)
    :rtype: :class:`obspy.core.event.Catalog`
    :returns: Figure object for further editing
    :rtype: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read, Catalog
    >>> from eqcorrscan.utils.sfile_util import read_event, readwavename
    >>> from eqcorrscan.utils.plotting import multi_event_singlechan
    >>> import glob
    >>> sfiles = glob.glob('eqcorrscan/tests/test_data/REA/TEST_/*')
    >>> catalog = Catalog()
    >>> streams = []
    >>> for sfile in sfiles:
    ...     catalog.append(read_event(sfile))
    ...     wavfile = readwavename(sfile)[0]
    ...     stream_path = 'eqcorrscan/tests/test_data/WAV/TEST_/' + wavfile
    ...     stream = read(stream_path)
    ...     # Annoying coping with seisan 2 letter channels
    ...     for tr in stream:
    ...         tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
    ...     streams.append(stream)
    >>> multi_event_singlechan(streams=streams, catalog=catalog,
    ...                        station='GCSZ', channel='EZ') # doctest: +SKIP

    .. image:: ../../plots/multi_event_singlechan.png
    """
    _check_save_args(save, savefile)
    # Work out how many picks we should have...
    short_cat = Catalog()
    short_streams = []
    for i, event in enumerate(catalog):
        event_stachans = [(pick.waveform_id.station_code,
                           pick.waveform_id.channel_code)
                          for pick in event.picks]
        if (station, channel) in event_stachans:
            short_cat.append(event)
            short_streams.append(streams[i])
    if len(short_cat) == 0:
        raise IOError('No picks for ' + station + ' ' + channel)
    traces = []
    al_traces = []
    if isinstance(short_streams, Stream):
        short_streams = [short_streams]
    st_list = deepcopy(short_streams)
    for i, event in enumerate(short_cat):
        # Extract the appropriate pick
        _pick = [pick for pick in event.picks if
                 pick.waveform_id.station_code == station and
                 pick.waveform_id.channel_code == channel]
        if len(_pick) == 0:
            print('No pick for channel')
            continue
        else:
            _pick = _pick[0]
        if st_list[i].select(station=station, channel=channel):
            tr = st_list[i].select(station=station, channel=channel)[0]
        else:
            print('No data for ' + _pick.waveform_id.station_code)
            continue
        tr.detrend('linear')
        if freqmin:
            tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        if realign:
            tr_cut = tr.copy()
            tr_cut.trim(_pick.time + cut[0],
                        _pick.time + cut[1],
                        nearest_sample=False)
            if len(tr_cut.data) <= (0.5 * (cut[1] - cut[0]) *
                                    tr_cut.stats.sampling_rate):
                msg = ''.join(['Not enough in the trace for ',
                               tr.stats.station,
                               '.', tr.stats.channel, '\n',
                               'Suggest removing pick from event at time ',
                               str(_pick.time)])
                warnings.warn(msg)
            else:
                al_traces.append(tr_cut)
        else:
            tr.trim(_pick.time - pre_pick,
                    _pick.time + clip - pre_pick,
                    nearest_sample=False)
        if len(tr.data) == 0:
            msg = ''.join(['No data in the trace for ', tr.stats.station,
                           '.', tr.stats.channel, '\n',
                           'Suggest removing pick from event at time ',
                           str(event.picks[0].time)])
            warnings.warn(msg)
            continue
        traces.append(tr)
    if realign:
        shift_len = int(0.25 * (cut[1] - cut[0]) *
                        al_traces[0].stats.sampling_rate)
        shifts = stacking.align_traces(al_traces, shift_len)
        for i in xrange(len(shifts)):
            print('Shifting by ' + str(shifts[i]) + ' seconds')
            _pick.time -= shifts[i]
            traces[i].trim(_pick.time - pre_pick,
                           _pick.time + clip - pre_pick,
                           nearest_sample=False)
    # We now have a list of traces
    if PWS:
        stack = 'PWS'
    else:
        stack = 'linstack'
    fig = multi_trace_plot(traces=traces, corr=True, stack=stack)
    if title:
        fig.suptitle(title)
    plt.subplots_adjust(hspace=0)
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return traces, short_cat, fig


def multi_trace_plot(traces, corr=True, stack='linstack', size=(7, 12),
                     show=True, title=None):
    """
    Plot multiple traces (usually from the same station) on the same plot.

    Differs somewhat to obspys stream.plot in that only relative time within \
    traces is worried about, it will not merge traces together.

    :type traces: list
    :param traces: List of obspy.core.Trace
    :type corr: bool
    :param corr: To calculate the correlation or not, if True, will add this \
        to the axes
    :type stack: str
    :param stack: To plot the stack as the first trace or not, select type of \
        stack: 'linstack' or 'PWS', or None.
    :type size: tuple
    :param size: Size of figure.
    :type show: bool
    :param show: Whether to plot the figure to screen or not.
    :type title: str
    :param title: Title to plot
    """
    if stack in ['linstack', 'PWS']:
        fig, axes = plt.subplots(len(traces) + 1, 1, sharex=True,
                                 figsize=size)
    else:
        fig, axes = plt.subplots(len(traces), 1, sharex=True,
                                 figsize=size)
    if len(traces) > 1:
        axes = axes.ravel()
    traces = [(trace, trace.stats.starttime.datetime) for trace in traces]
    traces.sort(key=lambda tup: tup[1])
    traces = [trace[0] for trace in traces]
    # Plot the traces
    for i, tr in enumerate(traces):
        y = tr.data
        x = np.arange(len(y))
        x = x / tr.stats.sampling_rate  # convert to seconds
        if not stack:
            ind = i
        else:
            ind = i + 1
        axes[ind].plot(x, y, 'k', linewidth=1.1)
        axes[ind].yaxis.set_ticks([])
    traces = [Stream(trace) for trace in traces]
    if stack == 'PWS':
        linstack = stacking.PWS_stack(traces)
    elif stack == 'linstack':
        linstack = stacking.linstack(traces)
    if stack in ['linstack', 'PWS']:
        tr = linstack[0]
        y = tr.data
        x = np.arange(len(y))
        x = x / tr.stats.sampling_rate
        axes[0].plot(x, y, 'r', linewidth=2.0)
        axes[0].set_ylabel('Stack', rotation=0)
        axes[0].yaxis.set_ticks([])
    for i, slave in enumerate(traces):
        if corr:
            cc = normxcorr2(tr.data, slave[0].data)
        if not stack:
            ind = i
        else:
            ind = i + 1
        if corr:
            axes[ind].set_ylabel('cc=' + str(round(np.max(cc), 2)), rotation=0)
        axes[ind].text(0.9, 0.15, str(round(np.max(slave[0].data))),
                       bbox=dict(facecolor='white', alpha=0.95),
                       transform=axes[ind].transAxes)
        axes[ind].text(0.7, 0.85, slave[0].stats.starttime.datetime.
                       strftime('%Y/%m/%d %H:%M:%S'),
                       bbox=dict(facecolor='white', alpha=0.95),
                       transform=axes[ind].transAxes)
    axes[-1].set_xlabel('Time (s)')
    if title:
        fig.suptitle(title)
    if show:
        plt.show()
    return fig


def detection_multiplot(stream, template, times, streamcolour='k',
                        templatecolour='r', save=False, savefile=None,
                        size=(10.5, 7.5), title=None):
    """
    Plot a stream of data with a template on top of it at detection times.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of data to be plotted as the background.
    :type template: obspy.core.stream.Stream
    :param template: Template to be plotted on top of the base stream.
    :type times: list
    :param times: list of times of detections in the order of the channels in
                template.
    :type streamcolour: str
    :param streamcolour: String of matplotlib colour types for the stream
    :type templatecolour: str
    :param templatecolour: Colour to plot the template in.
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True
    :type size: tuple
    :param size: Figure size.
    :type title: str
    :param title: Title for plot.

    :returns: :class:`matplotlib.figure.Figure`


    .. image:: ../../plots/detection_multiplot.png

    """
    _check_save_args(save, savefile)
    # Sort before plotting
    template = template.sort()
    # Only take traces that match in both
    template_stachans = [(tr.stats.station, tr.stats.channel)
                         for tr in template]
    stream = Stream([tr for tr in stream
                     if (tr.stats.station,
                         tr.stats.channel) in template_stachans])
    ntraces = min(len(template), len(stream))
    fig, axes = plt.subplots(ntraces, 1, sharex=True, figsize=size)
    if len(template) > 1:
        axes = axes.ravel()
    mintime = min([tr.stats.starttime for tr in template])
    for i, template_tr in enumerate(template):
        if len(template) > 1:
            axis = axes[i]
        else:
            axis = axes
        image = stream.select(station=template_tr.stats.station,
                              channel='*' + template_tr.stats.channel[-1])
        if not image:
            msg = ' '.join(['No data for', template_tr.stats.station,
                            template_tr.stats.channel])
            print(msg)
            continue
        image = image.merge()[0]
        # Downsample if needed
        if image.stats.sampling_rate > 20 and image.stats.npts > 10000:
            image.decimate(int(image.stats.sampling_rate // 20))
            template_tr.decimate(int(template_tr.stats.sampling_rate // 20))
        # Get a list of datetime objects
        image_times = [image.stats.starttime.datetime +
                       dt.timedelta((j * image.stats.delta) / 86400)
                       for j in range(len(image.data))]
        axis.plot(image_times, image.data / max(image.data),
                  streamcolour, linewidth=1.2)
        for k, time in enumerate(times):
            lagged_time = UTCDateTime(time) + (template_tr.stats.starttime -
                                               mintime)
            lagged_time = lagged_time.datetime
            template_times = [lagged_time +
                              dt.timedelta((j * template_tr.stats.delta) /
                                           86400)
                              for j in range(len(template_tr.data))]
            # Normalize the template according to the data detected in
            try:
                normalizer = max(image.data[int((template_times[0] -
                                                image_times[0]).
                                                total_seconds() /
                                                image.stats.delta):
                                            int((template_times[-1] -
                                                 image_times[0]).
                                                total_seconds() /
                                                image.stats.delta)] /
                                 max(image.data))
            except ValueError:
                # Occurs when there is no data in the image at this time...
                normalizer = max(image.data)
            normalizer /= max(template_tr.data)
            axis.plot(template_times,
                      template_tr.data * normalizer,
                      templatecolour, linewidth=1.2)
        ylab = '.'.join([template_tr.stats.station,
                         template_tr.stats.channel])
        axis.set_ylabel(ylab, rotation=0,
                        horizontalalignment='right')
        # axis.yaxis.set_ticks([])
    if len(template) > 1:
        axes[len(axes) - 1].set_xlabel('Time')
    else:
        axis.set_xlabel('Time')
    plt.subplots_adjust(hspace=0, left=0.175, right=0.95, bottom=0.07)
    plt.xticks(rotation=10)
    if title:
        plt.suptitle(title)
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def interev_mag_sfiles(sfiles, save=False, savefile=None, size=(10.5, 7.5)):
    """
    Plot inter-event time versus magnitude for series of events.

    Wrapper for :func:`eqcorrscan.utils.plotting.interev_mag`.

    :type sfiles: list
    :param sfiles: List of sfiles to read from
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True
    :type size: tuple
    :param size: Size of figure in inches.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> import glob
    >>> from eqcorrscan.utils.plotting import interev_mag_sfiles
    >>> sfiles = glob.glob('eqcorrscan/tests/test_data/REA/TEST_/*')
    >>> interev_mag_sfiles(sfiles=sfiles) # doctest: +SKIP

    .. plot::

        import glob, os
        from eqcorrscan.utils.plotting import interev_mag_sfiles
        sfiles = glob.glob(os.path.
                           realpath('../../../tests/test_data/REA/TEST_') +
                           os.sep + '*')
        print(sfiles)
        interev_mag_sfiles(sfiles=sfiles)
    """
    _check_save_args(save, savefile)
    times = []
    mags = []
    for sfile in sfiles:
        head = sfile_util.readheader(sfile)
        if head.preferred_origin():
            origin = head.preferred_origin()
        elif len(head.origins) > 0:
            origin = head.origins[0]
        else:
            origin = False
        if head.preferred_magnitude():
            magnitude = head.preferred_magnitude()
        elif len(head.magnitudes) > 0:
            magnitude = head.magnitudes[0]
        else:
            magnitude = False
        if origin and magnitude:
            times.append(origin.time)
            mags.append(magnitude.mag)
    fig = interev_mag(times=times, mags=mags, save=save, savefile=savefile,
                      size=size)
    return fig


def interev_mag(times, mags, save=False, savefile=None, size=(10.5, 7.5)):
    r"""Plot inter-event times against magnitude.

    :type times: list
    :param times: list of the detection times, must be sorted the same as mags
    :type mags: list
    :param mags: list of magnitudes
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True
    :type size: tuple
    :param size: Size of figure in inches.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.plotting import interev_mag
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime('2012-03-26T00:00:00')
    >>> t2 = t1 + (3 * 86400)
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
    >>> magnitudes = [event.preferred_magnitude().mag for event in catalog]
    >>> times = [event.preferred_origin().time for event in catalog]
    >>> interev_mag(times, magnitudes) # doctest: +SKIP

    .. plot::

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from eqcorrscan.utils.plotting import interev_mag
        client = Client('IRIS')
        t1 = UTCDateTime('2012-03-26T00:00:00')
        t2 = t1 + (3 * 86400)
        catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
        magnitudes = [event.preferred_magnitude().mag for event in catalog]
        times = [event.preferred_origin().time for event in catalog]
        interev_mag(times, magnitudes)
    """
    _check_save_args(save, savefile)
    l = [(times[i], mags[i]) for i in xrange(len(times))]
    l.sort(key=lambda tup: tup[0])
    times = [x[0] for x in l]
    mags = [x[1] for x in l]
    # Make two subplots next to each other of time before and time after
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=size)
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
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def obspy_3d_plot(inventory, catalog, save=False, savefile=None,
                  size=(10.5, 7.5)):
    """
    Plot obspy Inventory and obspy Catalog classes in three dimensions.

    :type inventory: obspy.core.inventory.inventory.Inventory
    :param inventory: Obspy inventory class containing station metadata
    :type catalog: obspy.core.event.catalog.Catalog
    :param catalog: Obspy catalog class containing event metadata
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True
    :type size: tuple
    :param size: Size of figure in inches.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example:

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.plotting import obspy_3d_plot
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> t2 = t1 + 86400
    >>> catalog = client.get_events(starttime=t1, endtime=t2, latitude=-43,
    ...                             longitude=170, maxradius=5)
    >>> inventory = client.get_stations(starttime=t1, endtime=t2, latitude=-43,
    ...                                 longitude=170, maxradius=10)
    >>> obspy_3d_plot(inventory=inventory, catalog=catalog) # doctest: +SKIP

    .. plot::

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from eqcorrscan.utils.plotting import obspy_3d_plot
        client = Client('IRIS')
        t1 = UTCDateTime(2012, 3, 26)
        t2 = t1 + 86400
        catalog = client.get_events(starttime=t1, endtime=t2, latitude=-43,
                                    longitude=170, maxradius=5)
        inventory = client.get_stations(starttime=t1, endtime=t2, latitude=-43,
                                        longitude=170, maxradius=10)
        obspy_3d_plot(inventory=inventory, catalog=catalog)
    """
    _check_save_args(save, savefile)
    nodes = []
    for ev in catalog:
        nodes.append((ev.preferred_origin().latitude,
                      ev.preferred_origin().longitude,
                      ev.preferred_origin().depth / 1000))
    # Will plot borehole instruments at elevation - depth if provided
    all_stas = []
    for net in inventory:
        for sta in net:
            if len(sta.channels) > 0:
                all_stas.append((sta.latitude, sta.longitude,
                                 sta.elevation / 1000 -
                                 sta.channels[0].depth / 1000))
            else:
                warnings.warn('No channel information attached, '
                              'setting elevation without depth')
                all_stas.append((sta.latitude, sta.longitude,
                                 sta.elevation / 1000))
    fig = threeD_seismplot(stations=all_stas, nodes=nodes, save=save,
                           savefile=savefile, size=size)
    return fig


def threeD_seismplot(stations, nodes, save=False, savefile=None,
                     size=(10.5, 7.5)):
    """
    Plot seismicity and stations in a 3D, movable, zoomable space.

    Uses matplotlibs Axes3D package.

    :type stations: list
    :param stations: list of one tuple per station of (lat, long, elevation), \
        with up positive.
    :type nodes: list
    :param nodes: list of one tuple per event of (lat, long, depth) with down \
        positive.
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True
    :type size: tuple
    :param size: Size of figure in inches.

    :returns: :class:`matplotlib.figure.Figure`

    .. Note::
        See :func:`eqcorrscan.utils.plotting.obspy_3d_plot` for example output.
    """
    _check_save_args(save, savefile)
    stalats, stalongs, staelevs = zip(*stations)
    evlats, evlongs, evdepths = zip(*nodes)
    # Cope with +/-180 latitudes...
    _evlongs = []
    for evlong in evlongs:
        if evlong < 0:
            evlong = float(evlong)
            evlong += 360
        _evlongs.append(evlong)
    evlongs = _evlongs
    _stalongs = []
    for stalong in stalongs:
        if stalong < 0:
            stalong = float(stalong)
            stalong += 360
        _stalongs.append(stalong)
    stalongs = _stalongs
    evdepths = [-1 * depth for depth in evdepths]
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(evlats, evlongs, evdepths, marker="x", c="k",
               label='Hypocenters')
    ax.scatter(stalats, stalongs, staelevs, marker="v", c="r",
               label='Stations')
    ax.set_ylabel("Longitude (deg)")
    ax.set_xlabel("Latitude (deg)")
    ax.set_zlabel("Elevation (km)")
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    plt.legend()
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def pretty_template_plot(template, size=(10.5, 7.5), save=False,
                         savefile=None, title=False, background=False,
                         picks=False):
    """
    Plot of a single template, possibly within background data.

    :type template: obspy.core.stream.Stream
    :param template: Template stream to plot
    :type size: tuple
    :param size: tuple of plot size
    :type save: bool
    :param save: if False will plot to screen, if True will save
    :type savefile: str
    :param savefile: String to save plot as, required if save=True.
    :type title: bool
    :param title: String if set will be the plot title
    :type background: obspy.core.stream.stream
    :param background: Stream to plot the template within.
    :type picks: list
    :param picks: List of :class:`obspy.core.event.origin.Pick` picks.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read
    >>> import os
    >>> from eqcorrscan.core import template_gen
    >>> from eqcorrscan.utils.plotting import pretty_template_plot
    >>> from eqcorrscan.utils import sfile_util
    >>>
    >>> test_file = os.path.join('eqcorrscan', 'tests', 'test_data', 'REA',
    ...                          'TEST_', '01-0411-15L.S201309')
    >>> test_wavefile = os.path.join('eqcorrscan', 'tests', 'test_data', 'WAV',
    ...                              'TEST_',
    ...                              '2013-09-01-0410-35.DFDPC_024_00')
    >>> event = sfile_util.readpicks(test_file)
    >>> st = read(test_wavefile)
    >>> st = st.filter('bandpass', freqmin=2.0, freqmax=15.0)
    >>> for tr in st:
    ...     tr = tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
    >>> template = template_gen._template_gen(event.picks, st, 2)
    >>> pretty_template_plot(template, background=st, # doctest +SKIP
    ...                      picks=event.picks) # doctest: +SKIP

    .. plot::

        from obspy import read
        from eqcorrscan.core import template_gen
        from eqcorrscan.utils.plotting import pretty_template_plot
        from eqcorrscan.utils import sfile_util
        import os
        test_file = os.path.realpath('../../..') + \
            '/tests/test_data/REA/TEST_/01-0411-15L.S201309'
        test_wavefile = os.path.realpath('../../..') +\
            '/tests/test_data/WAV/TEST_/' +\
            '2013-09-01-0410-35.DFDPC_024_00'
        event = sfile_util.readpicks(test_file)
        st = read(test_wavefile)
        st.filter('bandpass', freqmin=2.0, freqmax=15.0)
        for tr in st:
            tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
        template = template_gen._template_gen(event.picks, st, 2)
        pretty_template_plot(template, background=st,
                             picks=event.picks)
    """
    _check_save_args(save, savefile)
    fig, axes = plt.subplots(len(template), 1, sharex=True, figsize=size)
    if len(template) > 1:
        axes = axes.ravel()
    if not background:
        mintime = template.sort(['starttime'])[0].stats.starttime
    else:
        mintime = background.sort(['starttime'])[0].stats.starttime
    template.sort(['network', 'station', 'starttime'])
    lengths = []
    lines = []
    labels = []
    for i, tr in enumerate(template):
        # Cope with a single channel template case.
        if len(template) > 1:
            axis = axes[i]
        else:
            axis = axes
        delay = tr.stats.starttime - mintime
        y = tr.data
        x = np.linspace(0, (len(y) - 1) * tr.stats.delta, len(y))
        x += delay
        if background:
            btr = background.select(station=tr.stats.station,
                                    channel=tr.stats.channel)[0]
            bdelay = btr.stats.starttime - mintime
            by = btr.data
            bx = np.linspace(0, (len(by) - 1) * btr.stats.delta, len(by))
            bx += bdelay
            axis.plot(bx, by, 'k', linewidth=1)
            template_line, = axis.plot(x, y, 'r', linewidth=1.1,
                                       label='Template')
            if i == 0:
                lines.append(template_line)
                labels.append('Template')
            lengths.append(max(bx[-1], x[-1]))
        else:
            template_line, = axis.plot(x, y, 'k', linewidth=1.1,
                                       label='Template')
            if i == 0:
                lines.append(template_line)
                labels.append('Template')
            lengths.append(x[-1])
        # print(' '.join([tr.stats.station, str(len(x)), str(len(y))]))
        axis.set_ylabel('.'.join([tr.stats.station, tr.stats.channel]),
                        rotation=0, horizontalalignment='right')
        axis.yaxis.set_ticks([])
        # Plot the picks if they are given
        if picks:
            tr_picks = [pick for pick in picks if
                        pick.waveform_id.station_code == tr.stats.station and
                        pick.waveform_id.channel_code[0] +
                        pick.waveform_id.channel_code[-1] ==
                        tr.stats.channel[0] + tr.stats.channel[-1]]
            for pick in tr_picks:
                if not pick.phase_hint:
                    pcolor = 'k'
                    label = 'Unknown pick'
                elif 'P' in pick.phase_hint.upper():
                    pcolor = 'red'
                    label = 'P-pick'
                elif 'S' in pick.phase_hint.upper():
                    pcolor = 'blue'
                    label = 'S-pick'
                else:
                    pcolor = 'k'
                    label = 'Unknown pick'
                pdelay = pick.time - mintime
                # print(pdelay)
                line = axis.axvline(x=pdelay, color=pcolor, linewidth=2,
                                    linestyle='--', label=label)
                if label not in labels:
                    lines.append(line)
                    labels.append(label)
                # axes[i].plot([pdelay, pdelay], [])
    axis.set_xlim([0, max(lengths)])
    if len(template) > 1:
        axis = axes[len(template) - 1]
    else:
        axis = axes
    axis.set_xlabel('Time (s) from start of template')
    plt.figlegend(lines, labels, 'upper right')
    if title:
        if len(template) > 1:
            axes[0].set_title(title)
        else:
            axes.set_title(title)
    else:
        plt.subplots_adjust(top=0.98)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def plot_repicked(template, picks, det_stream, size=(10.5, 7.5), save=False,
                  savefile=None, title=False):
    """
    Plot a template over a detected stream, with picks corrected by lag-calc.

    :type template: obspy.core.stream.Stream
    :param template: Template used to make the detection, will be aligned \
        according to picks.
    :param picks:
        list of corrected picks, each pick must be an
        :class:`obspy.core.event.origin.Pick` object.
    :type picks: list
    :param det_stream: Stream to plot in the background, should be the \
        detection, data should encompass the time the picks are made.
    :type det_stream: obspy.core.stream.Stream
    :param size: Plot size.
    :type size: tuple
    :param save: To save figure or not, if false, will show to screen.
    :type save: bool
    :param savefile: File name to save file, required if save==True.
    :type savefile: str
    :param title: Title for plot, defaults to None.
    :type title: str

    :return: Figure handle which can be edited.
    :rtype: :class:`matplotlib.figure.Figure`

    .. image:: ../../plots/plot_repicked.png
    """
    _check_save_args(save, savefile)
    fig, axes = plt.subplots(len(template), 1, sharex=True, figsize=size)
    if len(template) > 1:
        axes = axes.ravel()
    mintime = det_stream.sort(['starttime'])[0].stats.starttime
    template.sort(['network', 'station', 'starttime'])
    lengths = []
    lines = []
    labels = []
    n_templates_plotted = 0
    for i, tr in enumerate(template.sort(['starttime'])):
        # Cope with a single channel template case.
        if len(template) > 1:
            axis = axes[i]
        else:
            axis = axes
        tr_picks = [pick for pick in picks if
                    pick.waveform_id.station_code == tr.stats.station and
                    pick.waveform_id.channel_code[0] +
                    pick.waveform_id.channel_code[-1] ==
                    tr.stats.channel[0] + tr.stats.channel[-1]]
        if len(tr_picks) > 1:
            msg = 'Multiple picks on channel %s' % tr.stats.station + ', ' + \
                  tr.stats.channel
            raise NotImplementedError(msg)
        if len(tr_picks) == 0:
            msg = 'No pick for chanel %s' % tr.stats.station + ', ' + \
                  tr.stats.channel
            print(msg)
        else:
            pick = tr_picks[0]
            delay = pick.time - mintime
            y = tr.data
            # Normalise
            y = y / max(y)
            x = np.linspace(0, (len(y) - 1) * tr.stats.delta, len(y))
            x += delay
        btr = det_stream.select(station=tr.stats.station,
                                channel=tr.stats.channel)[0]
        bdelay = btr.stats.starttime - mintime
        by = btr.data
        if len(tr_picks) > 0:
            by = by / max(by[int((delay - bdelay) * btr.stats.sampling_rate):
                             int((delay - bdelay) * btr.stats.sampling_rate) +
                             len(x)])
        else:
            by = by / max(by)
        bx = np.linspace(0, (len(by) - 1) * btr.stats.delta, len(by))
        bx += bdelay
        axis.plot(bx, by, 'k', linewidth=1.5)
        if len(tr_picks) > 0:
            template_line, = axis.plot(x, y, 'r', linewidth=1.6,
                                       label='Template')
            if not pick.phase_hint:
                pcolor = 'k'
                label = 'Unknown pick'
            elif 'P' in pick.phase_hint.upper():
                pcolor = 'red'
                label = 'P-pick'
            elif 'S' in pick.phase_hint.upper():
                pcolor = 'blue'
                label = 'S-pick'
            else:
                pcolor = 'k'
                label = 'Unknown pick'
            pdelay = pick.time - mintime
            line = axis.axvline(x=pdelay, color=pcolor, linewidth=2,
                                linestyle='--', label=label)
            if label not in labels:
                lines.append(line)
                labels.append(label)
            if n_templates_plotted == 0:
                lines.append(template_line)
                labels.append('Template')
            n_templates_plotted += 1
            lengths.append(max(bx[-1], x[-1]))
        else:
            lengths.append(bx[1])
        axis.set_ylabel('.'.join([tr.stats.station, tr.stats.channel]),
                        rotation=0, horizontalalignment='right')
        axis.yaxis.set_ticks([])
    axis.set_xlim([0, max(lengths)])
    if len(template) > 1:
        axis = axes[len(template) - 1]
    else:
        axis = axes
    axis.set_xlabel('Time (s) from %s' %
                    mintime.datetime.strftime('%Y/%m/%d %H:%M:%S.%f'))
    plt.figlegend(lines, labels, 'upper right')
    if title:
        if len(template) > 1:
            axes[0].set_title(title)
        else:
            axes.set_title(title)
    else:
        plt.subplots_adjust(top=0.98)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    if not save:
        plt.show()
        plt.close()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def NR_plot(stream, NR_stream, detections, false_detections=False,
            size=(18.5, 10), save=False, savefile=None, title=False):
    """
    Plot Network response alongside the stream used.

    Highlights detection times in the network response.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream to plot
    :type NR_stream: obspy.core.stream.Stream
    :param NR_stream: Stream for the network response
    :type detections: list
    :param detections: List of the detection time as :class:`datetime.datetime`
    :type false_detections: list
    :param false_detections:
        Either False (default) or list of false detection times
        (:class:`datetime.datetime`).
    :type size: tuple
    :param size: Size of figure, default is (18.5, 10)
    :type save: bool
    :param save:
        Save figure or plot to screen, if not False, must be string of save
        path.
    :type title: str
    :param title: String for the title of the plot, set to False

    :returns: :class:`matplotlib.figure.Figure`

    .. Note::
        Called by :mod:`eqcorrscan.core.bright_lights`, not a general use
        plot (hence no example)
    """
    _check_save_args(save, savefile)
    fig, axes = plt.subplots(len(stream) + 1, 1, sharex=True, figsize=size)
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
        axes[i].set_ylabel('.'.join([tr.stats.station, tr.stats.channel]),
                           rotation=0)
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
    axes[i].set_ylabel('.'.join([tr.stats.station, tr.stats.channel]),
                       rotation=0)
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
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def SVD_plot(SVStreams, SValues, stachans, title=False, save=False,
             savefile=None):
    """
    Plot singular vectors from the :mod:`eqcorrscan.utils.clustering` routines.

    One plot for each channel.

    :type SVStreams: list
    :param SVStreams:
        See :func:`eqcorrscan.utils.clustering.SVD_2_stream` - these should be
        ordered by power, e.g. first singular vector in the first stream.
    :type SValues: list
    :param SValues:
        List of floats of the singular values corresponding to the SVStreams
    :type stachans: list
    :param stachans: List of station.channel
    :type save: bool
    :param save:
        False will plot to screen, true will save plot and not show to screen.
    :type savefile: str
    :param savefile:
        Filename to save to, required for save=True, will label additionally
        according to station and channel.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read
    >>> import glob
    >>> from eqcorrscan.utils.plotting import SVD_plot
    >>> from eqcorrscan.utils.clustering import svd, SVD_2_stream
    >>> wavefiles = glob.glob('eqcorrscan/tests/test_data/WAV/TEST_/*')
    >>> streams = [read(w) for w in wavefiles[1:10]]
    >>> stream_list = []
    >>> for st in streams:
    ...     tr = st.select(station='GCSZ', channel='EHZ')
    ...     tr = tr.detrend('simple').resample(100).filter('bandpass',
    ...                                                    freqmin=2,
    ...                                                    freqmax=8)
    ...     stream_list.append(tr)
    >>> svec, sval, uvec, stachans = svd(stream_list=stream_list)
    >>> SVstreams = SVD_2_stream(SVectors=svec, stachans=stachans, k=3,
    ...                          sampling_rate=100)
    >>> SVD_plot(SVStreams=SVstreams, SValues=sval,
    ...          stachans=stachans) # doctest: +SKIP

    .. plot::

        from obspy import read
        import glob, os
        from eqcorrscan.utils.plotting import SVD_plot
        from eqcorrscan.utils.clustering import svd, SVD_2_stream
        wavefiles = glob.glob(os.path.realpath('../../..') +
                             '/tests/test_data/WAV/TEST_/*')
        streams = [read(w) for w in wavefiles[1:10]]
        stream_list = []
        for st in streams:
            tr = st.select(station='GCSZ', channel='EHZ')
            st.detrend('simple').resample(100).filter('bandpass', freqmin=5,
                                                      freqmax=40)
            stream_list.append(tr)
        svec, sval, uvec, stachans = svd(stream_list=stream_list)
        SVstreams = SVD_2_stream(SVectors=svec, stachans=stachans, k=3,
                                 sampling_rate=100)
        SVD_plot(SVStreams=SVstreams, SValues=sval,
                 stachans=stachans)
    """
    _check_save_args(save, savefile)
    for sval, stachan in zip(SValues, stachans):
        print(stachan)
        plot_traces = [SVStream.select(station=stachan.split('.')[0],
                                       channel=stachan.split('.')[1])[0]
                       for SVStream in SVStreams]
        fig, axes = plt.subplots(len(plot_traces), 1, sharex=True)
        if len(plot_traces) > 1:
            axes = axes.ravel()
        for i, tr in enumerate(plot_traces):
            y = tr.data
            x = np.linspace(0, len(y) * tr.stats.delta, len(y))
            axes[i].plot(x, y, 'k', linewidth=1.1)
            ylab = 'SV %s = %s' % (i + 1, round(sval[i] / len(sval), 2))
            axes[i].set_ylabel(ylab, rotation=0)
            axes[i].yaxis.set_ticks([])
            print(i)
        axes[-1].set_xlabel('Time (s)')
        plt.subplots_adjust(hspace=0)
        if title:
            axes[0].set_title(title)
        else:
            axes[0].set_title(stachan)
        if not save:
            plt.show()
        else:
            plt.savefig(savefile.split('.') + '_stachan.' +
                        savefile.split('.')[-1])
            plt.close()
    return fig


def plot_synth_real(real_template, synthetic, channels=False, size=(5, 10),
                    save=False, savefile=None):
    """
    Plot multiple channels of data for real data and synthetic.

    :type real_template: obspy.core.stream.Stream
    :param real_template: Stream of the real template
    :type synthetic: obspy.core.stream.Stream
    :param synthetic: Stream of synthetic template
    :type channels: list
    :param channels: List of tuples of (station, channel) to plot, default is \
            False, which plots all.
    :type size: tuple
    :param size: Plot size.
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True

    :returns: :class:`matplotlib.figure.Figure`

    >>> from obspy import read, Stream, Trace
    >>> from eqcorrscan.utils.synth_seis import seis_sim
    >>> from eqcorrscan.utils.plotting import plot_synth_real
    >>> real = read()
    >>> synth = Stream(Trace(seis_sim(sp=100, flength=200)))
    >>> synth[0].stats.station = 'RJOB'
    >>> synth[0].stats.channel = 'EHZ'
    >>> synth[0].stats.sampling_rate = 100
    >>> synth = synth.filter('bandpass', freqmin=2, freqmax=8)
    >>> real = real.select(station='RJOB',
    ...                    channel='EHZ').detrend('simple').filter('bandpass',
    ...                                                            freqmin=2,
    ...                                                            freqmax=8)
    >>> real = real.trim(starttime=real[0].stats.starttime + 43,
    ...                  endtime=real[0].stats.starttime +
    ...                  45).detrend('simple')
    >>> plot_synth_real(real_template=real, synthetic=synth,
    ...                 size=(7, 4)) # doctest: +SKIP

    .. plot::

        from eqcorrscan.utils.plotting import plot_synth_real
        from obspy import read, Stream, Trace
        from eqcorrscan.utils.synth_seis import seis_sim
        import os
        real = read()
        synth = Stream(Trace(seis_sim(SP=100, flength=200)))
        synth[0].stats.station = 'RJOB'
        synth[0].stats.channel = 'EHZ'
        synth[0].stats.sampling_rate = 100
        synth.filter('bandpass', freqmin=2, freqmax=8)
        real = real.select(station='RJOB', channel='EHZ').detrend('simple').\
            filter('bandpass', freqmin=2, freqmax=8)
        real.trim(starttime=real[0].stats.starttime + 4.9,
                  endtime=real[0].stats.starttime + 6.9).detrend('simple')
        plot_synth_real(real_template=real, synthetic=synth, size=(7, 4))
    """
    _check_save_args(save, savefile)
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
    fig, axes = plt.subplots(len(stachans), 1, sharex=True, figsize=size)
    if len(stachans) > 1:
        axes = axes.ravel()
    for i, stachan in enumerate(stachans):
        if len(stachans) > 1:
            axis = axes[i]
        else:
            axis = axes
        real_tr = real_template.select(station=stachan[0],
                                       channel=stachan[1])[0]
        synth_tr = synthetic.select(station=stachan[0],
                                    channel=stachan[1])[0]
        shift, corr = xcorr(real_tr, synth_tr, 2)
        print('Shifting by: ' + str(shift) + ' samples')
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
            axis.plot(x, y, colours[j], linewidth=2.0, label=labels[j])
            axis.get_yaxis().set_ticks([])
        ylab = stachan[0] + '.' + stachan[1] + ' cc=' + str(round(corr, 2))
        axis.set_ylabel(ylab, rotation=0)
    plt.subplots_adjust(hspace=0)
    # axes[0].legend()
    if len(stachans) > 1:
        axes[-1].set_xlabel('Time (s)')
    else:
        axis.set_xlabel('Time (s)')
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def freq_mag(magnitudes, completeness, max_mag, binsize=0.2, save=False,
             savefile=None):
    """
    Plot a frequency-magnitude histogram and cumulative density plot.

    Currently this will compute a b-value, for a goiven completeness.
    B-value is computed by linear fitting to section of curve between
    completeness and max_mag.

    :type magnitudes: list
    :param magnitudes: list of float of magnitudes
    :type completeness: float
    :param completeness: Level to compute the b-value above
    :type max_mag: float
    :param max_mag: Maximum magnitude to try and fit a b-value to
    :type binsize: float
    :param binsize: Width of histogram bins, defaults to 0.2
    :type save: bool
    :param save: False will plot to screen, true will save plot and not show \
        to screen.
    :type savefile: str
    :param savefile: Filename to save to, required for save=True

    :returns: :class:`matplotlib.figure.Figure`

    .. Note::
        See :func:`eqcorrscan.utils.mag_calc.calc_b_value` for a least-squares
        method of estimating completeness and b-value. For estimating maximum
        curvature see :func:`eqcorrscan.utils.mag_calc.calc_max_curv`.

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.plotting import freq_mag
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime('2012-03-26T00:00:00')
    >>> t2 = t1 + (3 * 86400)
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
    >>> magnitudes = [event.preferred_magnitude().mag for event in catalog]
    >>> freq_mag(magnitudes, completeness=4, max_mag=7) # doctest: +SKIP

    .. plot::

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from eqcorrscan.utils.plotting import freq_mag
        client = Client('IRIS')
        t1 = UTCDateTime('2012-03-26T00:00:00')
        t2 = t1 + (3 * 86400)
        catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
        magnitudes = [event.preferred_magnitude().mag for event in catalog]
        freq_mag(magnitudes, completeness=4, max_mag=7)
    """
    _check_save_args(save, savefile)
    # Ensure magnitudes are sorted
    magnitudes.sort()
    # Check that there are no nans or infs
    if np.isnan(magnitudes).any():
        warnings.warn('Found nan values, removing them')
        magnitudes = [mag for mag in magnitudes if not np.isnan(mag)]
    if np.isinf(magnitudes).any():
        warnings.warn('Found inf values, removing them')
        magnitudes = [mag for mag in magnitudes if not np.isinf(mag)]
    fig, ax1 = plt.subplots()
    # Set up the bins, the bin-size could be a variables
    bins = np.arange(int(min(magnitudes) - 1), int(max(magnitudes) + 1),
                     binsize)
    n, bins, patches = ax1.hist(magnitudes, bins, facecolor='Black',
                                alpha=0.5, label='Magnitudes')
    ax1.set_ylabel('Frequency')
    ax1.set_ylim([0, max(n) + 0.5 * max(n)])
    plt.xlabel('Magnitude')
    # Now make the cumulative density function
    counts = Counter(magnitudes)
    cdf = np.zeros(len(counts))
    mag_steps = np.zeros(len(counts))
    for i, magnitude in enumerate(sorted(counts.keys(), reverse=True)):
        mag_steps[i] = magnitude
        if i > 0:
            cdf[i] = cdf[i - 1] + counts[magnitude]
        else:
            cdf[i] = counts[magnitude]
    ax2 = ax1.twinx()
    # ax2.scatter(magnitudes, np.log10(cdf), c='k', marker='+', s=20, lw=2,
    ax2.scatter(mag_steps, np.log10(cdf), c='k', marker='+', s=20, lw=2,
                label='Magnitude cumulative density')
    # Now we want to calculate the b-value and plot the fit
    x = []
    y = []
    for i, magnitude in enumerate(mag_steps):
        if magnitude >= completeness <= max_mag:
            x.append(magnitude)
            y.append(cdf[i])
    fit = np.polyfit(x, np.log10(y), 1)
    fit_fn = np.poly1d(fit)
    ax2.plot(magnitudes, fit_fn(magnitudes), '--k',
             label='GR trend, b-value = ' + str(abs(fit[0]))[0:4] +
             '\n $M_C$ = ' + str(completeness))
    ax2.set_ylabel('$Log_{10}$ of cumulative density')
    plt.xlim([min(magnitudes) - 0.1, max(magnitudes) + 0.2])
    plt.ylim([min(np.log10(cdf)) - 0.5, max(np.log10(cdf)) + 1.0])
    plt.legend(loc=2)
    if not save:
        plt.show()
    else:
        plt.savefig(savefile)
        plt.close()
    return fig


def spec_trace(traces, cmap=None, wlen=0.4, log=False, trc='k', tralpha=0.9,
               size=(10, 13), fig=None, title=None, show=True):
    """
    Plots seismic data with spectrogram behind.

    Takes a stream or list of traces and plots the trace with the spectra
    beneath it.

    :type traces: list
    :param traces: Traces to be plotted, can be a single
        :class:`obspy.core.stream.Stream`, or a list of
        :class:`obspy.core.trace.Trace`.
    :type cmap: str
    :param cmap:
        `Matplotlib colormap
        <http://matplotlib.org/examples/color/colormaps_reference.html>`_.
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
    :type fig: matplotlib.figure.Figure
    :param fig: Figure to plot onto, defaults to self generating.
    :type title: str
    :param title: Title of plot.
    :type show: bool
    :param show: To show plot or not, if false, will return Fig.

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.plotting import spec_trace
    >>> st = read()
    >>> spec_trace(st, trc='white') # doctest: +SKIP


    .. plot::

        from obspy import read
        from eqcorrscan.utils.plotting import spec_trace
        st = read()
        spec_trace(st, trc='white')

    """
    if isinstance(traces, Stream):
        traces.sort(['station', 'channel'])
    if not fig:
        fig = plt.figure()
    for i, tr in enumerate(traces):
        if i == 0:
            ax = fig.add_subplot(len(traces), 1, i + 1)
        else:
            ax = fig.add_subplot(len(traces), 1, i + 1, sharex=ax)
        ax1, ax2 = _spec_trace(tr, cmap=cmap, wlen=wlen, log=log, trc=trc,
                               tralpha=tralpha, axes=ax)
        ax.set_yticks([])
        if i < len(traces) - 1:
            plt.setp(ax1.get_xticklabels(), visible=False)
        if type(traces) == list:
            ax.text(0.005, 0.85, tr.stats.starttime.datetime.
                    strftime('%Y/%m/%d %H:%M:%S'),
                    bbox=dict(facecolor='white', alpha=0.8),
                    transform=ax2.transAxes)
        else:
            ax.text(0.005, 0.85, '.'.join([tr.stats.station,
                                           tr.stats.channel]),
                    bbox=dict(facecolor='white', alpha=0.8),
                    transform=ax2.transAxes)
        ax.text(0.005, 0.02, str(np.max(tr.data).round(1)),
                bbox=dict(facecolor='white', alpha=0.95),
                transform=ax2.transAxes)
    ax.set_xlabel('Time (s)')
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(w=size[0], h=size[1], forward=True)
    fig.text(0.04, 0.5, 'Frequency (Hz)', va='center', rotation='vertical')
    if title:
        plt.suptitle(title)
    if show:
        plt.show()
    else:
        return fig


def _spec_trace(trace, cmap=None, wlen=0.4, log=False, trc='k',
                tralpha=0.9, size=(10, 2.5), axes=None, title=None):
    """
    Function to plot a trace over that traces spectrogram.

    Uses obspys spectrogram routine.

    :type trace: obspy.core.trace.Trace
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
        # Fig.close()
    else:
        return ax1, ax2


if __name__ == "__main__":
    import doctest
    doctest.testmod()
