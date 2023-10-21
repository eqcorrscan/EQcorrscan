"""
Utility code for most of the plots used as part of the EQcorrscan package.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import logging
import datetime as dt
import copy
import os

import matplotlib.dates as mdates
from copy import deepcopy
from collections import Counter
from itertools import cycle
from scipy.linalg import diagsvd
from scipy import fftpack
from obspy import UTCDateTime, Stream, Catalog, Trace
from obspy.signal.cross_correlation import correlate, xcorr_max

from eqcorrscan.utils.stacking import align_traces, PWS_stack, linstack


Logger = logging.getLogger(__name__)


# A wrapper to add the same docs everywhere

def additional_docstring(**kwargs):
    def _wrapper(target):
        target.__doc__ = target.__doc__.format(**kwargs)
        return target
    return _wrapper


plotting_kwargs = """
    :type title: str
    :param title: Title of figure
    :type show: bool
    :param show: Whether to show the figure or not (defaults to True)
    :type save: bool
    :param save: Whether to save the figure or not (defaults to False)
    :type savefile: str
    :param savefile:
        Filename to save figure to, if `save==True` (defaults to
        "EQcorrscan_figure.png")
    :type return_figure: bool
    :param return_figure:
        Whether to return the figure or not (defaults to True), if False
        then the figure will be cleared and closed.
    :type size: tuple of float
    :param size: Figure size as (width, height) in inches. Defaults to
        (10.5, 7.5)"""


@additional_docstring(plotting_kwargs=plotting_kwargs)
def _finalise_figure(fig, **kwargs):  # pragma: no cover
    """
    Internal function to wrap up a figure.
    {plotting_kwargs}
    """
    import matplotlib.pyplot as plt

    title = kwargs.get("title")
    show = kwargs.get("show", True)
    save = kwargs.get("save", False)
    savefile = kwargs.get("savefile", "EQcorrscan_figure.png")
    return_figure = kwargs.get("return_figure", False)
    size = kwargs.get("size", (10.5, 7.5))
    fig.set_size_inches(size)
    if title:
        fig.suptitle(title)
    if save:
        fig.savefig(savefile, bbox_inches="tight")
        Logger.info("Saved figure to {0}".format(savefile))
    if show:
        plt.show(block=True)
    if return_figure:
        return fig
    fig.clf()
    plt.close(fig)
    return None


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


@additional_docstring(plotting_kwargs=plotting_kwargs)
def xcorr_plot(template, image, shift=None, cc=None, cc_vec=None, **kwargs):
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
    {plotting_kwargs}

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
    import matplotlib.pyplot as plt
    if cc is None or shift is None:
        if not isinstance(cc_vec, np.ndarray):
            Logger.error('Given cc: %s and shift: %s' % (cc, shift))
            raise IOError('Must provide either cc_vec, or cc and shift')
        shift = np.abs(cc_vec).argmax()
        cc = cc_vec[shift]
    x = np.arange(len(image))
    plt.plot(x, image / abs(image).max(), 'k', lw=1.3, label='Image')
    x = np.arange(len(template)) + shift
    plt.plot(x, template / abs(template).max(), 'r', lw=1.1, label='Template')
    plt.title('Shift=%s, Correlation=%s' % (shift, cc))
    fig = plt.gcf()
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def triple_plot(cccsum, cccsum_hist, trace, threshold, **kwargs):
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
    {plotting_kwargs}

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
    >>> triple_plot(cccsum=ccc[0], cccsum_hist=ccc[0], trace=tr,
    ...             threshold=0.8) # doctest: +SKIP


    .. image:: ../../plots/triple_plot.png
    """
    import matplotlib.pyplot as plt
    if len(cccsum) != len(trace.data):
        Logger.error(
            'cccsum is: ' + str(len(cccsum)) + ' trace is: ' +
            str(len(trace.data)))
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
    ax3.hist(cccsum_hist, 200, density=True, histtype='stepfilled',
             orientation='horizontal', color='black')
    ax3.set_ylim([-5, 5])
    fig = plt.gcf()
    fig.suptitle(trace.id)
    fig.canvas.draw()
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def peaks_plot(data, starttime, samp_rate, peaks=None, **kwargs):
    """
    Plot peaks to check that the peak finding routine is running correctly.

    Used in debugging for the EQcorrscan module.

    :type data: numpy.array
    :param data: Numpy array of the data within which peaks have been found
    :type starttime: obspy.core.utcdatetime.UTCDateTime
    :param starttime: Start time for the data
    :type samp_rate: float
    :param samp_rate: Sampling rate of data in Hz
    :type peaks: list
    :param peaks: List of tuples of peak locations and amplitudes (loc, amp)
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> import numpy as np
    >>> from eqcorrscan.utils import findpeaks
    >>> from eqcorrscan.utils.plotting import peaks_plot
    >>> from obspy import UTCDateTime
    >>> data = np.random.randn(200)
    >>> data[30] = 100
    >>> data[60] = 40
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
    import matplotlib.pyplot as plt
    peaks = peaks or [(0, 0)]
    npts = len(data)
    t = np.arange(npts, dtype=np.float32) / (samp_rate * 3600)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, data, 'k')
    ax1.scatter(peaks[0][1] / (samp_rate * 3600), peaks[0][0],
                color='r', label='Peaks')
    for peak in peaks:
        ax1.scatter(peak[1] / (samp_rate * 3600), peak[0], color='r')
    ax1.legend()
    ax1.set_xlabel("Time after %s [hr]" % starttime.isoformat())
    ax1.axis('tight')
    fig.suptitle('Peaks')
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def cumulative_detections(dates=None, template_names=None, detections=None,
                          plot_grouped=False, group_name=None, rate=False,
                          binsize=None, plot_legend=True, ax=None, **kwargs):
    """
    Plot cumulative detections or detection rate in time.

    Simple plotting function to take a list of either datetime objects or
    :class:`eqcorrscan.core.match_filter.Detection` objects and plot
    a cumulative detections list.  Can take dates as a list of lists and will
    plot each list separately, e.g. if you have dates from more than one
    template it will overlay them in different colours.

    :type dates: list
    :param dates: Must be a list of lists of datetime.datetime objects
    :type template_names: list
    :param template_names: List of the template names in order of the dates
    :type detections: list
    :param detections: List of :class:`eqcorrscan.core.match_filter.Detection`
    :type plot_grouped: bool
    :param plot_grouped:
        Plot detections for each template individually, or group them all
        together - set to False (plot template detections individually) by
        default.
    :type group_name: str
    :param group_name:
        Name to put in legend for the group, only used if `plot_grouped=True`
    :type rate: bool
    :param rate:
        Whether or not to plot the rate of detection per day. Only works for
        plot_grouped=True
    :type binsize: int
    :param binsize: Bin size for rate plotting in seconds.
    :type plot_legend: bool
    :param plot_legend:
        Specify whether to plot legend of template names. Defaults to True.
    :type ax: `matplotlib.pyplot.Axis`
    :param ax: Axis to plot into, if you want to re-use a figure.
    {plotting_kwargs}


    :returns: :class:`matplotlib.figure.Figure`

    .. note::
        Can either take lists of
        :class:`eqcorrscan.core.match_filter.Detection` objects directly, or
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

    .. rubric:: Example 2: Rate plotting

    >>> import datetime as dt
    >>> import numpy as np
    >>> from eqcorrscan.utils.plotting import cumulative_detections
    >>> dates = []
    >>> for i in range(3):
    ...     dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
    ...                   for n in np.random.randn(100)])
    >>> cumulative_detections(dates, ['a', 'b', 'c'], plot_grouped=True,
    ...                       rate=True, show=True) # doctest: +SKIP

    .. plot::

        import datetime as dt
        import numpy as np
        from eqcorrscan.utils.plotting import cumulative_detections
        dates = []
        for i in range(3):
            dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
                          for n in np.random.randn(100)])
        cumulative_detections(dates, ['a', 'b', 'c'], plot_grouped=True,
                              rate=True, show=True)

    """
    import matplotlib.pyplot as plt
    from eqcorrscan.core.match_filter import Detection
    # Set up a default series of parameters for lines
    colors = cycle(['red', 'green', 'blue', 'cyan', 'magenta', 'yellow',
                    'black', 'firebrick', 'purple', 'darkgoldenrod', 'gray'])
    linestyles = cycle(['-', '-.', '--', ':'])
    # Check that dates is a list of lists
    if detections is None:
        if not isinstance(dates, list):
            raise IndexError("No detections or dates given")
        if not isinstance(dates[0], list):
            dates = [dates]
    else:
        dates = []
        template_names = []
        for detection in detections:
            if not type(detection) == Detection:
                raise IOError(
                    'detection not of type: eqcorrscan.core.match_filter'
                    '.Detection')
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
        if group_name is not None:
            template_names = group_name
        else:
            template_names = ['All templates']
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
    else:
        fig = ax.get_figure()
    # Make sure not to pad at edges
    ax.margins(0, 0)
    min_date = min([min(_d) for _d in dates])
    max_date = max([max(_d) for _d in dates])
    for k, template_dates in enumerate(dates):
        template_dates.sort()
        plot_dates = deepcopy(template_dates)
        plot_dates.insert(0, min_date)
        plot_dates.insert(-1, template_dates[-1])
        color = next(colors)
        if color == 'red':
            linestyle = next(linestyles)
        counts = np.arange(-1, len(template_dates) + 1)
        if rate:
            if not plot_grouped:
                msg = 'Plotting rate only implemented for plot_grouped=True'
                raise NotImplementedError(msg)
            if binsize is None:
                if 365 <= (max_date - min_date).days:
                    binsize = 7 * 86400
                    ax.set_ylabel('Detections per week')
                elif 31 < (max_date - min_date).days < 365:
                    binsize = 86400
                    ax.set_ylabel('Detections per day')
                elif 1 < (max_date - min_date).days <= 31:
                    binsize = 6 * 3600
                    ax.set_ylabel('Detections per 6 hours')
                elif 3600 < (max_date - min_date).total_seconds() <= 86400:
                    binsize = 900
                    ax.set_ylabel('Detections per 15 minutes')
                else:
                    binsize = 60
                    ax.set_ylabel('Detections per minute')
            else:
                ax.set_ylabel('Detections per {0} seconds'.format(binsize))
            bins = np.arange(
                min_date, max_date + dt.timedelta(seconds=binsize),
                dt.timedelta(seconds=binsize))
            ax.hist(mdates.date2num(plot_dates), bins=mdates.date2num(bins),
                    label='Rate of detections', color='darkgrey',
                    alpha=0.5)
        else:
            ax.plot(plot_dates, counts, linestyle,
                    color=color, label=template_names[k],
                    linewidth=2.0, drawstyle='steps')
            ax.set_ylabel('Cumulative detections')
    ax.set_xlabel('Date')
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
        mins = mdates.MinuteLocator(byminute=np.arange(0, 60, 10))
    elif 7200 <= timedif.total_seconds() < 10800:
        hours = mdates.MinuteLocator(byminute=[0, 15, 30, 45])
        mins = mdates.MinuteLocator(byminute=np.arange(0, 60, 5))
    elif timedif.total_seconds() <= 1200:
        hours = mdates.MinuteLocator(byminute=np.arange(0, 60, 2))
        mins = mdates.MinuteLocator(byminute=np.arange(0, 60, 0.5))
    elif 25200 < timedif.total_seconds() <= 86400:
        hours = mdates.HourLocator(byhour=np.arange(0, 24, 3))
        mins = mdates.HourLocator(byhour=np.arange(0, 24, 1))
    elif 86400 < timedif.total_seconds() <= 172800:
        hours = mdates.HourLocator(byhour=np.arange(0, 24, 6))
        mins = mdates.HourLocator(byhour=np.arange(0, 24, 1))
    elif timedif.total_seconds() > 172800:
        hours = mdates.AutoDateLocator()
        mins = mdates.HourLocator(byhour=np.arange(0, 24, 3))
    else:
        hours = mdates.MinuteLocator(byminute=np.arange(0, 60, 5))
    # Minor locator overruns maxticks for ~year-long datasets
    if timedif.total_seconds() < 172800:
        ax.xaxis.set_minor_locator(mins)
        hrFMT = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
    else:
        hrFMT = mdates.DateFormatter('%Y/%m/%d')
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hrFMT)
    fig.autofmt_xdate()
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=15)
    if not rate:
        ax.set_ylim([0, max([len(_d) for _d in dates])])
    if plot_legend:
        if ax.legend() is not None:
            leg = ax.legend(loc=2, prop={'size': 8}, ncol=2)
            leg.get_frame().set_alpha(0.5)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def threeD_gridplot(nodes, **kwargs):
    """Plot in a series of grid points in 3D.

    :type nodes: list
    :param nodes: List of tuples of the form (lat, long, depth)
    {plotting_kwargs}

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
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    import matplotlib.pyplot as plt
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
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def multi_event_singlechan(streams, catalog, station, channel,
                           clip=10.0, pre_pick=2.0,
                           freqmin=False, freqmax=False, realign=False,
                           cut=(-3.0, 5.0), PWS=False, **kwargs):
    """
    Plot data from a single channel for multiple events.

    Data will be aligned by their pick-time given in the appropriate picks.
    Requires an individual stream for each event you want to plot,
    events are stored in the :class:`obspy.core.event.Catalog` object, and
    there must be picks present for the streams you wish to plot.  Events will
    be aligned if `realign=True`, in this case the traces will be aligned
    using the window defined by `cut`.

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
    {plotting_kwargs}

    :returns: Aligned and cut :class:`obspy.core.trace.Trace`
    :rtype: list
    :returns:
        New picks in based on alignment (if alignment is performed, if not
        will return the same as input)
    :rtype: :class:`obspy.core.event.Catalog`
    :returns: Figure object for further editing
    :rtype: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read, Catalog, read_events
    >>> from obspy.io.nordic.core import readwavename
    >>> from eqcorrscan.utils.plotting import multi_event_singlechan
    >>> import glob
    >>> sfiles = glob.glob('eqcorrscan/tests/test_data/REA/TEST_/*.S??????')
    >>> catalog = Catalog()
    >>> streams = []
    >>> for sfile in sfiles:
    ...     catalog += read_events(sfile)
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
    import matplotlib.pyplot as plt
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
    al_picks = []
    if isinstance(short_streams, Stream):
        short_streams = [short_streams]
    st_list = deepcopy(short_streams)
    Logger.debug(short_cat)
    for i, event in enumerate(short_cat):
        # Extract the appropriate pick
        _pick = [pick for pick in event.picks if
                 pick.waveform_id.station_code == station and
                 pick.waveform_id.channel_code == channel]
        if len(_pick) == 0:
            Logger.info('No pick for channel')
            continue
        else:
            _pick = _pick[0]
        if st_list[i].select(station=station, channel=channel):
            tr = st_list[i].select(station=station, channel=channel)[0]
        else:
            Logger.info('No data for ' + _pick.waveform_id.station_code)
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
                Logger.warning(msg)
            else:
                al_traces.append(tr_cut)
                al_picks.append(_pick)
        else:
            tr.trim(_pick.time - pre_pick,
                    _pick.time + clip - pre_pick,
                    nearest_sample=False)
        if len(tr.data) == 0:
            msg = ''.join(['No data in the trace for ', tr.stats.station,
                           '.', tr.stats.channel, '\n',
                           'Suggest removing pick from event at time ',
                           str(event.picks[0].time)])
            Logger.warning(msg)
            continue
        traces.append(tr)
    if realign:
        shift_len = int(0.25 * (cut[1] - cut[0]) *
                        al_traces[0].stats.sampling_rate)
        shifts = align_traces(al_traces, shift_len)[0]
        for i in range(len(shifts)):
            Logger.info('Shifting by ' + str(shifts[i]) + ' seconds')
            _pick.time -= shifts[i]
            traces[i].trim(al_picks[i].time - pre_pick,
                           al_picks[i].time + clip - pre_pick,
                           nearest_sample=True)
    # We now have a list of traces
    if PWS:
        stack = 'PWS'
    else:
        stack = 'linstack'
    for tr in traces:
        Logger.debug(tr)
    fig = multi_trace_plot(
        traces=traces, corr=True, stack=stack, show=False, return_figure=True)
    plt.subplots_adjust(hspace=0)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return traces, short_cat, fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def multi_trace_plot(traces, corr=True, stack='linstack', **kwargs):
    """
    Plot multiple traces (usually from the same station) on the same plot.

    Differs somewhat to obspy's stream.plot in that only relative time within
    traces is worried about, it will not merge traces together.

    :type traces: list
    :param traces: List of obspy.core.Trace
    :type corr: bool
    :param corr:
        To calculate the correlation or not, if True, will add this to the
        axes
    :type stack: str
    :param stack:
        To plot the stack as the first trace or not, select type of
         stack: 'linstack' or 'PWS', or None.
    {plotting_kwargs}
    """
    import matplotlib.pyplot as plt
    from eqcorrscan.core.match_filter import normxcorr2
    n_axes = len(traces)
    if stack in ['linstack', 'PWS']:
        n_axes += 1
    fig, axes = plt.subplots(n_axes, 1, sharex=True)
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
    if stack in ['linstack', 'PWS']:
        if stack == "PWS":
            tr = PWS_stack(traces)[0]
        else:
            tr = linstack(traces)[0]
        y = tr.data
        x = np.arange(len(y))
        x = x / tr.stats.sampling_rate
        axes[0].plot(x, y, 'r', linewidth=2.0)
        axes[0].set_ylabel('Stack', rotation=0)
        axes[0].yaxis.set_ticks([])
    else:
        tr = traces[0]
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
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def detection_multiplot(stream, template, times, streamcolour='k',
                        templatecolour='r', **kwargs):
    """
    Plot a stream of data with a template on top of it at detection times.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of data to be plotted as the background.
    :type template: obspy.core.stream.Stream
    :param template: Template to be plotted on top of the base stream.
    :type times: list
    :param times: list of detection times, one for each event
    :type streamcolour: str
    :param streamcolour: String of matplotlib colour types for the stream
    :type templatecolour: str
    :param templatecolour: Colour to plot the template in.
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read, read_events
    >>> import os
    >>> from eqcorrscan.core import template_gen
    >>> from eqcorrscan.utils.plotting import detection_multiplot
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>>
    >>> test_file = os.path.join(TEST_PATH, 'REA',
    ...                          'TEST_', '01-0411-15L.S201309')
    >>> test_wavefile = os.path.join(
    ...     TEST_PATH, 'WAV', 'TEST_', '2013-09-01-0410-35.DFDPC_024_00')
    >>> event = read_events(test_file)[0]
    >>> st = read(test_wavefile)
    >>> st = st.filter('bandpass', freqmin=2.0, freqmax=15.0)
    >>> for tr in st:
    ...     tr = tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
    ...     # Hack around seisan 2-letter channel naming
    ...     tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
    >>> template = template_gen._template_gen(event.picks, st, 2)
    >>> times = [min([pk.time -0.05 for pk in event.picks])]
    >>> detection_multiplot(stream=st, template=template,
    ...                     times=times) # doctest: +SKIP

    .. plot::

        from obspy import read, read_events
        import os
        from eqcorrscan.core import template_gen
        from eqcorrscan.utils.plotting import detection_multiplot
        test_file = os.path.realpath('../../..') + \
            '/tests/test_data/REA/TEST_/01-0411-15L.S201309'
        test_wavefile = os.path.realpath('../../..') +\
            '/tests/test_data/WAV/TEST_/' +\
            '2013-09-01-0410-35.DFDPC_024_00'
        event = read_events(test_file)[0]
        st = read(test_wavefile)
        st.filter('bandpass', freqmin=2.0, freqmax=15.0)
        for tr in st:
            tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        template = template_gen._template_gen(event.picks, st, 2)
        times = [min([pk.time -0.05 for pk in event.picks])]
        detection_multiplot(stream=st, template=template,
                            times=times)

    """
    import matplotlib.pyplot as plt
    # Only take traces that match in both accounting for streams shorter than
    # templates
    template_stachans = [(tr.stats.station, tr.stats.channel)
                         for tr in template]
    stream_stachans = [(tr.stats.station, tr.stats.channel)
                       for tr in stream]
    temp = Stream([tr for tr in template
                   if (tr.stats.station,
                       tr.stats.channel) in stream_stachans])
    st = Stream([tr for tr in stream
                 if (tr.stats.station,
                     tr.stats.channel) in template_stachans])
    ntraces = len(temp)
    fig, axes = plt.subplots(ntraces, 1, sharex=True)
    if len(temp) > 1:
        axes = axes.ravel()
    mintime = min([tr.stats.starttime for tr in temp])
    temp.sort(keys=['starttime'])
    for i, template_tr in enumerate(temp):
        if len(temp) > 1:
            axis = axes[i]
        else:
            axis = axes
        image = st.select(station=template_tr.stats.station,
                          channel='*' + template_tr.stats.channel[-1])
        if not image:
            msg = ' '.join(['No data for', template_tr.stats.station,
                            template_tr.stats.channel])
            Logger.info(msg)
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
        for time in times:
            lagged_time = UTCDateTime(time) + (template_tr.stats.starttime -
                                               mintime)
            lagged_time = lagged_time.datetime
            template_times = [lagged_time +
                              dt.timedelta((j * template_tr.stats.delta) /
                                           86400)
                              for j in range(len(template_tr.data))]
            # Normalize the template according to the data detected in
            try:
                normalizer = max(
                    image.data[int(
                        (template_times[0] - image_times[0]
                         ).total_seconds() / image.stats.delta):
                               int(
                        (template_times[-1] - image_times[0]
                         ).total_seconds() / image.stats.delta)
                    ] / max(image.data))
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
    if len(template) > 1:
        axes[len(axes) - 1].set_xlabel('Time')
    else:
        axis.set_xlabel('Time')
    plt.subplots_adjust(hspace=0, left=0.175, right=0.95, bottom=0.07)
    plt.xticks(rotation=10)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def interev_mag(times, mags, **kwargs):
    """
    Plot inter-event times against magnitude.

    :type times: list
    :param times: list of the detection times, must be sorted the same as mags
    :type mags: list
    :param mags: list of magnitudes
    {plotting_kwargs}

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
    import matplotlib.pyplot as plt
    info = [(times[i], mags[i]) for i in range(len(times))]
    info.sort(key=lambda tup: tup[0])
    times = [x[0] for x in info]
    mags = [x[1] for x in info]
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
    axes[0].autoscale(enable=True, tight=True)
    axes[1].autoscale(enable=True, tight=True)
    plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=30)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def obspy_3d_plot(inventory, catalog, **kwargs):
    """
    Plot obspy Inventory and obspy Catalog classes in three dimensions.

    :type inventory: obspy.core.inventory.inventory.Inventory
    :param inventory: Obspy inventory class containing station metadata
    :type catalog: obspy.core.event.catalog.Catalog
    :param catalog: Obspy catalog class containing event metadata
    {plotting_kwargs}

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
                Logger.warning('No channel information attached, '
                               'setting elevation without depth')
                all_stas.append((sta.latitude, sta.longitude,
                                 sta.elevation / 1000))
    fig = threeD_seismplot(
        stations=all_stas, nodes=nodes, **kwargs)
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def threeD_seismplot(stations, nodes, **kwargs):
    """
    Plot seismicity and stations in a 3D, movable, zoomable space.

    Uses matplotlibs Axes3D package.

    :type stations: list
    :param stations:
        list of one tuple per station of (lat, long, elevation), with up
        positive.
    :type nodes: list
    :param nodes:
        list of one tuple per event of (lat, long, depth) with down positive.
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. Note::
        See :func:`eqcorrscan.utils.plotting.obspy_3d_plot` for example output.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
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
    fig = plt.figure()
    ax = Axes3D(fig)
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
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def noise_plot(signal, noise, normalise=False, **kwargs):
    """
    Plot signal and noise fourier transforms and the difference.

    :type signal: `obspy.core.stream.Stream`
    :param signal: Stream of "signal" window
    :type noise: `obspy.core.stream.Stream`
    :param noise: Stream of the "noise" window.
    :type normalise: bool
    :param normalise: Whether to normalise the data before plotting or not.
    {plotting_kwargs}

    :return: `matplotlib.pyplot.Figure`
    """
    import matplotlib.pyplot as plt

    # Work out how many traces we can plot
    n_traces = 0
    for tr in signal:
        try:
            noise.select(id=tr.id)[0]
        except IndexError:  # pragma: no cover
            continue
        n_traces += 1

    fig, axes = plt.subplots(n_traces, 2, sharex=True, sharey="col")
    if len(signal) > 1:
        axes = axes.ravel()
    i = 0
    lines = []
    labels = []
    for tr in signal:
        try:
            noise_tr = noise.select(id=tr.id)[0]
        except IndexError:  # pragma: no cover
            continue
        ax1 = axes[i]
        ax2 = axes[i + 1]
        fft_len = fftpack.next_fast_len(
            max(noise_tr.stats.npts, tr.stats.npts))
        if not normalise:
            signal_fft = fftpack.rfft(tr.data, fft_len)
            noise_fft = fftpack.rfft(noise_tr.data, fft_len)
        else:
            signal_fft = fftpack.rfft(tr.data / max(tr.data), fft_len)
            noise_fft = fftpack.rfft(
                noise_tr.data / max(noise_tr.data), fft_len)
        frequencies = np.linspace(0, 1 / (2 * tr.stats.delta), fft_len // 2)
        noise_line, = ax1.semilogy(
            frequencies, 2.0 / fft_len * np.abs(noise_fft[0: fft_len // 2]),
            'k', label="noise")
        signal_line, = ax1.semilogy(
            frequencies, 2.0 / fft_len * np.abs(signal_fft[0: fft_len // 2]),
            'r', label="signal")
        if "signal" not in labels:
            labels.append("signal")
            lines.append(signal_line)
        if "noise" not in labels:
            labels.append("noise")
            lines.append(noise_line)
        ax1.set_ylabel(tr.id, rotation=0, horizontalalignment='right')
        ax2.plot(
            frequencies,
            (2.0 / fft_len * np.abs(signal_fft[0: fft_len // 2])) -
            (2.0 / fft_len * np.abs(noise_fft[0: fft_len // 2])), 'k')
        ax2.yaxis.tick_right()
        ax2.set_ylim(bottom=1e-6)
        i += 2
    axes[-1].set_xlabel("Frequency (Hz)")
    axes[-2].set_xlabel("Frequency (Hz)")
    axes[0].set_title("Spectra")
    axes[1].set_title("Signal - noise")
    fig.legend(lines, labels, 'upper left')
    fig.subplots_adjust(hspace=0, top=0.91)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def pretty_template_plot(template, background=False, event=False,
                         sort_by="distance", **kwargs):
    """
    Plot of a single template, possibly within background data.

    :type template: obspy.core.stream.Stream
    :param template: Template stream to plot
    :type background: obspy.core.stream.stream
    :param background: Stream to plot the template within.
    :type event: obspy.core.event.event.Event
    :param event:
        Event object containing picks, and optionally information on the origin
        and arrivals. When supplied, function tries to extract hypocentral
        distance from origin/arrivals, to sort the template traces by
        hypocentral distance.
    :type sort_by: string
    :param sort_by:
        "distance" (default) or "pick_time" (not relevant if no event supplied)
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read, read_events
    >>> import os
    >>> from eqcorrscan.core import template_gen
    >>> from eqcorrscan.utils.plotting import pretty_template_plot
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>>
    >>> test_file = os.path.join(TEST_PATH, 'REA', 'TEST_',
    ...                          '01-0411-15L.S201309')
    >>> test_wavefile = os.path.join(
    ...     TEST_PATH, 'WAV', 'TEST_', '2013-09-01-0410-35.DFDPC_024_00')
    >>> event = read_events(test_file)[0]
    >>> st = read(test_wavefile)
    >>> st = st.filter('bandpass', freqmin=2.0, freqmax=15.0)
    >>> for tr in st:
    ...     tr = tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
    ...     # Hack around seisan 2-letter channel naming
    ...     tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
    >>> template = template_gen._template_gen(event.picks, st, 2)
    >>> pretty_template_plot(template, background=st, # doctest +SKIP
    ...                      event=event) # doctest: +SKIP

    .. plot::

        from obspy import read, read_events
        from eqcorrscan.core import template_gen
        from eqcorrscan.utils.plotting import pretty_template_plot
        import os
        # Get the path to the test data
        import eqcorrscan
        import os
        TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
        test_file = os.path.join(
            TEST_PATH, 'REA', 'TEST_', '01-0411-15L.S201309')
        test_wavefile = os.path.join(
            TEST_PATH, 'WAV', 'TEST_', '2013-09-01-0410-35.DFDPC_024_00')
        event = read_events(test_file)[0]
        st = read(test_wavefile)
        st.filter('bandpass', freqmin=2.0, freqmax=15.0)
        for tr in st:
            tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        template = template_gen._template_gen(event.picks, st, 2)
        pretty_template_plot(template, background=st, event=event)
    """
    from eqcorrscan.utils.catalog_utils import get_ordered_trace_indices
    import matplotlib.pyplot as plt
    picks = kwargs.get("picks", None)
    if picks:
        Logger.warning(
            "The picks argument is depreciated, please use a full event. "
            "This argument will be removed in future versions.")
    elif event:
        picks = event.picks
    fig, axes = plt.subplots(len(template), 1, sharex=True)
    if len(template) > 1:
        axes = axes.ravel()
    if not background:
        mintime = template.sort(['starttime'])[0].stats.starttime
    else:
        mintime = background.sort(['starttime'])[0].stats.starttime
    template.sort(['network', 'station', 'starttime'])
    trace_indices_ordered = get_ordered_trace_indices(template, event=event,
                                                      sort_by=sort_by)
    if trace_indices_ordered is not None:
        template = Stream([template[j] for j in trace_indices_ordered])
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
    axes[0].legend(lines, labels, loc='upper right', framealpha=1)
    axes[0].set_zorder(2)
    title = kwargs.get("title") or None
    if title:
        if len(template) > 1:
            axes[0].set_title(title)
        else:
            axes.set_title(title)
        kwargs.pop("title")  # Do not give title to _finalise_figure
    else:
        plt.subplots_adjust(top=0.98)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def plot_repicked(template, picks, det_stream, **kwargs):
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
    {plotting_kwargs}

    :return: Figure handle which can be edited.
    :rtype: :class:`matplotlib.figure.Figure`

    .. image:: ../../plots/plot_repicked.png
    """
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(template), 1, sharex=True)
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
            Logger.info(msg)
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
            template_line, = axis.plot(
                x, y, 'r', linewidth=1.6, label='Template')
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
        axis_legend = axes[0]
    else:
        axis = axes
        axis_legend = axes
    axis.set_xlabel('Time (s) from %s' %
                    mintime.datetime.strftime('%Y/%m/%d %H:%M:%S.%f'))
    axis_legend.legend(lines, labels, loc='upper right', framealpha=1)
    axis_legend.set_zorder(2)
    title = kwargs.get("title") or None
    if title:
        if len(template) > 1:
            axes[0].set_title(title)
        else:
            axes.set_title(title)
        kwargs.pop("title")
    else:
        plt.subplots_adjust(top=0.98)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def svd_plot(svstreams, svalues, stachans, **kwargs):
    """
    Plot singular vectors from the :mod:`eqcorrscan.utils.clustering` routines.

    One plot for each channel.

    :type svstreams: list
    :param svstreams:
        See :func:`eqcorrscan.utils.clustering.svd_to_stream` - these should be
        ordered by power, e.g. first singular vector in the first stream.
    :type svalues: list
    :param svalues:
        List of floats of the singular values corresponding to the SVStreams
    :type stachans: list
    :param stachans: List of station.channel
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. rubric:: Example

    >>> from obspy import read
    >>> import glob
    >>> from eqcorrscan.utils.plotting import svd_plot
    >>> from eqcorrscan.utils.clustering import svd, svd_to_stream
    >>> wavefiles = glob.glob('eqcorrscan/tests/test_data/WAV/TEST_/2013-*')
    >>> streams = [read(w) for w in wavefiles[1:10]]
    >>> stream_list = []
    >>> for st in streams:
    ...     tr = st.select(station='GCSZ', channel='EHZ')
    ...     tr = tr.detrend('simple').resample(100).filter(
    ...        'bandpass', freqmin=2, freqmax=8)
    ...     stream_list.append(tr)
    >>> uvec, sval, svec, stachans = svd(stream_list=stream_list)
    >>> svstreams = svd_to_stream(uvectors=uvec, stachans=stachans, k=3,
    ...                           sampling_rate=100)
    >>> svd_plot(svstreams=svstreams, svalues=sval,
    ...          stachans=stachans) # doctest: +SKIP

    .. plot::

        from obspy import read
        import glob, os
        from eqcorrscan.utils.plotting import svd_plot
        from eqcorrscan.utils.clustering import svd, svd_to_stream
        wavefiles = glob.glob(os.path.realpath('../../..') +
                             '/tests/test_data/WAV/TEST_/2013-*')
        streams = [read(w) for w in wavefiles[1:10]]
        stream_list = []
        for st in streams:
            tr = st.select(station='GCSZ', channel='EHZ')
            st.detrend('simple').resample(100).filter('bandpass', freqmin=5,
                                                      freqmax=40)
            stream_list.append(tr)
        svec, sval, uvec, stachans = svd(stream_list=stream_list)
        svstreams = svd_to_stream(uvectors=uvec, stachans=stachans, k=3,
                                  sampling_rate=100)
        svd_plot(svstreams=svstreams, svalues=sval,
                 stachans=stachans)
    """
    import matplotlib.pyplot as plt
    figures = []
    for sval, stachan in zip(svalues, stachans):
        Logger.info(stachan)
        plot_traces = [SVStream.select(station=stachan[0],
                                       channel=stachan[1])[0]
                       for SVStream in svstreams]
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
            Logger.debug(i)
        axes[-1].set_xlabel('Time (s)')
        plt.subplots_adjust(hspace=0)
        fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
        figures.append(fig)
    return figures


@additional_docstring(plotting_kwargs=plotting_kwargs)
def plot_synth_real(real_template, synthetic, channels=False, **kwargs):
    """
    Plot multiple channels of data for real data and synthetic.

    :type real_template: obspy.core.stream.Stream
    :param real_template: Stream of the real template
    :type synthetic: obspy.core.stream.Stream
    :param synthetic: Stream of synthetic template
    :type channels: list
    :param channels: List of tuples of (station, channel) to plot, default is \
            False, which plots all.
    {plotting_kwargs}

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
    >>> real = real.select(
    ...    station='RJOB', channel='EHZ').detrend('simple').filter(
    ...       'bandpass', freqmin=2, freqmax=8)
    >>> real = real.trim(
    ...    starttime=real[0].stats.starttime + 4.9,
    ...    endtime=real[0].stats.starttime + 6.9).detrend('simple')
    >>> plot_synth_real(real_template=real, synthetic=synth,
    ...                 size=(7, 4)) # doctest: +SKIP

    .. plot::

        from eqcorrscan.utils.plotting import plot_synth_real
        from obspy import read, Stream, Trace
        from eqcorrscan.utils.synth_seis import seis_sim
        import os
        real = read()
        synth = Stream(Trace(seis_sim(sp=100, flength=200)))
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
    import matplotlib.pyplot as plt
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
    fig, axes = plt.subplots(len(stachans), 1, sharex=True)
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
        corr_fun = correlate(real_tr, synth_tr, 2)
        shift, corr = xcorr_max(corr_fun)
        shift = int(shift)
        Logger.info('Shifting by: ' + str(shift) + ' samples')
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
        axes.set_xlabel('Time (s)')
    if "size" not in kwargs.keys():
        kwargs.update({"size": (5, 10)})  # Backwards compat
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def freq_mag(magnitudes, completeness, max_mag, binsize=0.2, **kwargs):
    """
    Plot a frequency-magnitude histogram and cumulative density plot.

    Currently this will compute a b-value, for a given completeness.
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
    {plotting_kwargs}

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
    import matplotlib.pyplot as plt
    # Ensure magnitudes are sorted
    magnitudes.sort()
    # Check that there are no nans or infs
    if np.isnan(magnitudes).any():
        Logger.warning('Found nan values, removing them')
        magnitudes = [mag for mag in magnitudes if not np.isnan(mag)]
    if np.isinf(magnitudes).any():
        Logger.warning('Found inf values, removing them')
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
        if completeness <= magnitude <= max_mag:
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
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def spec_trace(traces, cmap=None, wlen=0.4, log=False, trc='k', tralpha=0.9,
               fig=None, **kwargs):
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
    :param fig: Figure to plot onto, defaults to self generating.
    {plotting_kwargs}

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
    import matplotlib.pyplot as plt
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
        if i < len(traces) - 1:
            plt.setp(ax1.get_xticklabels(), visible=False)
        if isinstance(traces, list):
            ax.text(0.005, 0.85, "{0}::{1}".format(tr.id, tr.stats.starttime),
                    bbox=dict(facecolor='white', alpha=0.8),
                    transform=ax2.transAxes)
        elif isinstance(traces, Stream):
            ax.text(0.005, 0.85, tr.id,
                    bbox=dict(facecolor='white', alpha=0.8),
                    transform=ax2.transAxes)
        ax.text(0.005, 0.02, str(np.max(tr.data).round(1)),
                bbox=dict(facecolor='white', alpha=0.95),
                transform=ax2.transAxes)
    ax.set_xlabel('Time (s)')
    fig.subplots_adjust(hspace=0)
    fig.text(0.04, 0.5, 'Frequency (Hz)', va='center', rotation='vertical')
    if "size" not in kwargs.keys():
        kwargs.update({"size": (10, 13)})  # backwards compat
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
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
    import matplotlib.pyplot as plt
    if not axes:
        fig = plt.figure(figsize=size)
        ax1 = fig.add_subplot(111)
    else:
        ax1 = axes
    trace.spectrogram(wlen=wlen, log=log, show=False, cmap=cmap, axes=ax1)
    fig = plt.gcf()
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
        fig.set_size_inches(size)
        fig.show()
    else:
        return ax1, ax2


@additional_docstring(plotting_kwargs=plotting_kwargs)
def subspace_detector_plot(detector, stachans, **kwargs):
    """
    Plotting for the subspace detector class.

    Plot the output basis vectors for the detector at the given dimension.

    Corresponds to the first n horizontal vectors of the V matrix.

    :type detector: :class:`eqcorrscan.core.subspace.Detector`
    :type stachans: list
    :param stachans: list of tuples of station, channel pairs to plot.
    :type stachans: list
    :param stachans: List of tuples of (station, channel) to use.  Can set\
        to 'all' to use all the station-channel pairs available. If \
        detector is multiplexed, will just plot that.
    {plotting_kwargs}

    :returns: Figure
    :rtype: matplotlib.pyplot.Figure

    .. rubric:: Example

    >>> from eqcorrscan.core import subspace
    >>> import os
    >>> detector = subspace.Detector()
    >>> detector.read(os.path.join(
    ...    os.path.abspath(os.path.dirname(__file__)),
    ...    '..', 'tests', 'test_data', 'subspace',
    ...    'stat_test_detector.h5'))
    Detector: Tester
    >>> subspace_detector_plot(detector=detector, stachans='all', size=(10, 7),
    ...                        show=True) # doctest: +SKIP

    .. plot::

        from eqcorrscan.core import subspace
        from eqcorrscan.utils.plotting import subspace_detector_plot
        import os
        print('running subspace plot')
        detector = subspace.Detector()
        detector.read(os.path.join('..', '..', '..', 'tests', 'test_data',
                                   'subspace', 'stat_test_detector.h5'))
        subspace_detector_plot(detector=detector, stachans='all', size=(10, 7),
                               show=True)
    """
    import matplotlib.pyplot as plt
    if stachans == 'all' and not detector.multiplex:
        stachans = detector.stachans
    elif detector.multiplex:
        stachans = [('multi', ' ')]
    if np.isinf(detector.dimension):
        msg = ' '.join(['Infinite subspace dimension. Only plotting as many',
                        'dimensions as events in design set'])
        Logger.warning(msg)
        nrows = detector.v[0].shape[1]
    else:
        nrows = detector.dimension
    fig, axes = plt.subplots(nrows=nrows, ncols=len(stachans),
                             sharex=True, sharey=True)
    x = np.arange(len(detector.u[0]), dtype=np.float32)
    if detector.multiplex:
        x /= len(detector.stachans) * detector.sampling_rate
    else:
        x /= detector.sampling_rate
    for column, stachan in enumerate(stachans):
        channel = detector.u[column]
        for row, vector in enumerate(channel.T[0:nrows]):
            if len(stachans) == 1:
                if nrows == 1:
                    axis = axes
                else:
                    axis = axes[row]
            else:
                axis = axes[row, column]
            if row == 0:
                axis.set_title('.'.join(stachan))
            axis.plot(x, vector, 'k', linewidth=1.1)
            if column == 0:
                axis.set_ylabel('Basis %s' % (row + 1), rotation=0)
            if row == nrows - 1:
                axis.set_xlabel('Time (s)')
            axis.set_yticks([])
    plt.subplots_adjust(hspace=0.05)
    plt.subplots_adjust(wspace=0.05)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def subspace_fc_plot(detector, stachans, **kwargs):
    """
    Plot the fractional energy capture of the detector for all events in
    the design set

    :type detector: :class:`eqcorrscan.core.subspace.Detector`
    :type stachans: list
    :param stachans: list of tuples of station, channel pairs to plot.
    :type stachans: list
    :param stachans: List of tuples of (station, channel) to use.  Can set\
        to 'all' to use all the station-channel pairs available. If \
        detector is multiplexed, will just plot that.
    {plotting_kwargs}

    :returns: Figure
    :rtype: matplotlib.pyplot.Figure

    .. rubric:: Example

    >>> from eqcorrscan.core import subspace
    >>> import os
    >>> detector = subspace.Detector()
    >>> detector.read(os.path.join(os.path.abspath(os.path.dirname(__file__)),
    ...                            '..', 'tests', 'test_data', 'subspace',
    ...                            'stat_test_detector.h5'))
    Detector: Tester
    >>> subspace_fc_plot(detector=detector, stachans='all', size=(10, 7),
    ...                        show=True) # doctest: +SKIP

    .. plot::

        from eqcorrscan.core import subspace
        from eqcorrscan.utils.plotting import subspace_detector_plot
        import os
        print('running subspace plot')
        detector = subspace.Detector()
        detector.read(os.path.join('..', '..', '..', 'tests', 'test_data',
                                   'subspace', 'stat_test_detector.h5'))
        subspace_fc_plot(detector=detector, stachans='all', size=(10, 7),
                               show=True)

    """
    import matplotlib.pyplot as plt
    if stachans == 'all' and not detector.multiplex:
        stachans = detector.stachans
    elif detector.multiplex:
        stachans = [('multi', ' ')]
    # Work out how many rows and columns are most 'square'
    pfs = []
    for x in range(1, len(stachans)):
        if len(stachans) % x == 0:
            pfs.append(x)
    if stachans == [('multi', ' ')]:
        ncols = 1
    else:
        ncols = min(pfs,
                    key=lambda x: abs((np.floor(np.sqrt(len(stachans))) - x)))
    nrows = len(stachans) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True,
                             sharey=True, squeeze=False)
    for column, axis in enumerate(axes.reshape(-1)):
        axis.set_title('.'.join(stachans[column]))
        sig = diagsvd(detector.sigma[column], detector.u[column].shape[0],
                      detector.v[column].shape[0])
        A = np.dot(sig, detector.v[column])  # v is v.H from scipy.svd
        if detector.dimension > max(
                detector.v[column].shape) or detector.dimension == np.inf:
            dim = max(detector.v[column].shape) + 1
        else:
            dim = detector.dimension + 1
        av_fc_dict = {i: [] for i in range(dim)}
        for ai in A.T:
            fcs = []
            for j in range(dim):
                av_fc_dict[j].append(float(np.dot(ai[:j].T, ai[:j])))
                fcs.append(float(np.dot(ai[:j].T, ai[:j])))
            axis.plot(fcs, color='grey')
        avg = [np.average(_dim[1]) for _dim in av_fc_dict.items()]
        axis.plot(avg, color='red', linewidth=3.)
        if column % ncols == 0 or column == 0:
            axis.set_ylabel('Frac. E Capture (Fc)')
        if column + 1 > len(stachans) - ncols:
            axis.set_xlabel('Subspace Dimension')
    plt.subplots_adjust(hspace=0.2)
    plt.subplots_adjust(wspace=0.2)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig


@additional_docstring(plotting_kwargs=plotting_kwargs)
def twoD_seismplot(catalog=None, locations=None, bgcolor='#909090',
                   method='depth', **kwargs):
    """
    Plot seismicity in a 2D map with two cross section along latitude and
    longitude.

    :type catalog: obspy.core.event.catalog.Catalog
    :param catalog: Obspy catalog class containing event metadata
    :type locations: list
    :param locations:
        list of one tuple per event of (lat, long, depth, time) with
        down positive.
    :type bgcolor: string
    :param bgcolor: Background's color of map and sections.
        all name or RGB code that acceptable in matplotlib.
    :type method: string
    :param method:
        making color palette of locations according to 'depth', 'time' or
        'sequence'.
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. note::
        If each location doesn't have time or depth, set them to zero.
    .. note::
        kwargs accepts all option that available in
        `matplotlib.axes.Axes.scatter`.
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    assert (catalog and locations) or catalog or locations,\
        "Requires catalog and/or locations"
    # set default parameters of plt.scatter()
    default_parameters = {'cmap': 'viridis_r', 'marker': ',', 's': 1, 'lw': 1}
    for key in default_parameters.keys():
        if key not in kwargs.keys():
            kwargs[key] = default_parameters[key]
    # get parameters of _finalise_figure
    _kwargs = {}
    for key in ['title', 'show', 'save', 'savefile', 'return_figure', 'size']:
        if key in kwargs.keys():
            _kwargs[key] = kwargs[key]
            del kwargs[key]
    # making coordinates
    locations = locations or []
    msg = "An event of the catalog got ignored, because it didn't have origin"
    if catalog:
        for event in catalog:
            try:
                origin = event.preferred_origin() or event.origins[0]
            except IndexError:  # No origin found
                Logger.warning(msg)
                continue
            _lat = origin.latitude
            _lon = origin.longitude
            _dep = origin.depth / 1000
            _time = origin.time
            locations.append((_lat, _lon, _dep, _time))
    # sort location according to method
    if method in ['time', 'sequence']:
        locations.sort(key=lambda ind: ind[3])
    elif method == 'depth':
        locations.sort(reverse=False, key=lambda ind: ind[2])
    lat, lon, dep, time = zip(*locations)
    if method == 'depth':
        c0, c1, c2 = dep, lon, lat
        label0, label1, label2 = 'Depth (km)', 'Longitude', 'Latitude'
    elif method == 'time':
        dt = [t - time[0] for t in time]
        c0 = c1 = c2 = dt
        label = f'Origin-time offset from {time[0]} (s)'
    elif method == 'sequence':
        c0 = c1 = c2 = range(len(dep))
        label = 'Event number'
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[3, 1],
                           wspace=0.01, hspace=0.01)
    # map view
    ax0 = plt.subplot(gs[0])
    ax0.set_facecolor(bgcolor)
    ax0.set_ylabel('Latitude')
    ax0.set_xticks([])
    map0 = ax0.scatter(lon, lat, c=c0, **kwargs)
    # cross section parallel to latitude (lat ,depth)
    ax1 = fig.add_subplot(gs[1])
    ax1.set_facecolor(bgcolor)
    ax1.set_yticks([])
    ax1.set_xlabel('Depth (km)')
    map1 = ax1.scatter(dep, lat, c=c1, **kwargs)
    # cross section parallel to longitude (lon ,depth)
    ax2 = plt.subplot(gs[2])
    ax2.set_facecolor(bgcolor)
    ax2.invert_yaxis()
    ax2.set_ylabel('Depth (km)')
    ax2.set_xlabel('Longitude')
    map2 = ax2.scatter(lon, dep, c=c2, **kwargs)
    # location of color bar
    if method == 'depth':
        #
        divider0 = make_axes_locatable(ax0)
        cax0 = divider0.append_axes("top", size="4%", pad="2%")
        cbar0 = fig.colorbar(map0, ax=ax0, cax=cax0, orientation="horizontal")
        cbar0.set_label(label0, rotation=0, labelpad=-45, y=1.05)
        cax0.xaxis.set_ticks_position("top")
        #
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("top", size="4%", pad="2%")
        cbar1 = fig.colorbar(map1, ax=ax1, cax=cax1, orientation="horizontal")
        cbar1.set_label(label1, rotation=0, labelpad=-45, y=1.03)
        cax1.xaxis.set_ticks_position("top")
        #
        divider2 = make_axes_locatable(ax2)
        cax2 = divider2.append_axes("bottom", size="7%", pad="35%")
        cbar2 = fig.colorbar(map2, ax=ax2, cax=cax2, orientation="horizontal",
                             pad=0.7)
        cbar2.set_label(label2, rotation=0, labelpad=-8, x=1.02)
        ax2.xaxis.set_label_coords(1.02, -0.1)
    elif method == 'time' or method == 'sequence':
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("right", size="4%", pad="2%")
        cbar1 = fig.colorbar(map1, ax=ax1, cax=cax1, orientation="vertical")
        cbar1.set_label(label)
    fig = _finalise_figure(fig=fig, **_kwargs)  # pragma: no cover
    return fig


def _match_filter_plot(stream, cccsum, template_names, rawthresh, plotdir,
                       plot_format, i):  # pragma: no cover
    """
    Plotting function for match_filter.

    :param stream: Stream to plot
    :param cccsum: Cross-correlation sum to plot
    :param template_names: Template names used
    :param rawthresh: Threshold level
    :param plotdir: Location to save plots
    :param plot_format: Output plot type (e.g. png, svg, eps, pdf...)
    :param i: Template index name to plot.
    """
    import matplotlib.pyplot as plt
    tr = stream[0]
    pad_len = len(tr.data) - len(cccsum)
    cccsum = np.pad(cccsum, (0, pad_len))
    if plotdir is not None:
        plt.ioff()
    stream_plot = copy.deepcopy(tr)
    # Downsample for plotting
    stream_plot = _plotting_decimation(stream_plot, 10e5, 4)
    samp_rate = stream_plot.stats.sampling_rate
    cccsum_plot = Trace(cccsum)
    cccsum_plot.stats.sampling_rate = tr.stats.sampling_rate
    # Resample here to maintain shape better
    cccsum_hist = cccsum_plot.copy()
    cccsum_hist = _plotting_decimation(cccsum_hist, 10e5, 4).data
    cccsum_plot = chunk_data(cccsum_plot, samp_rate, 'Maxabs').data
    # Enforce same length
    stream_plot.data = stream_plot.data[0:len(cccsum_plot)]
    cccsum_plot = cccsum_plot[0:len(stream_plot.data)]
    cccsum_hist = cccsum_hist[0:len(stream_plot.data)]
    plot_name = "{0}/cccsum_plot_{1}_{2}.{3}".format(
        plotdir, template_names[i], stream[0].stats.starttime.strftime(
                  "%Y-%m-%dT%H%M%S"), plot_format)
    plot_kwargs = dict(show=True)
    if plotdir is not None:
        if not os.path.isdir(plotdir):
            os.makedirs(plotdir)
        plot_kwargs.update(dict(show=False, save=True, savefile=plot_name))
    triple_plot(cccsum=cccsum_plot, cccsum_hist=cccsum_hist,
                trace=stream_plot, threshold=rawthresh, **plot_kwargs)


def _plotting_decimation(trace, max_len=10e5, decimation_step=4):
    """
    Decimate data until required length reached.

    :type trace: obspy.core.stream.Trace
    :param trace: Trace to decimate
    type max_len: int
    :param max_len: Maximum length in samples
    :type decimation_step: int
    :param decimation_step: Decimation factor to use for each step.

    :return: obspy.core.stream.Trace

    .. rubric: Example

    >>> from obspy import Trace
    >>> import numpy as np
    >>> trace = Trace(np.random.randn(1000))
    >>> trace = _plotting_decimation(trace, max_len=100, decimation_step=2)
    >>> print(trace.stats.npts)
    63
    """
    trace_len = trace.stats.npts
    while trace_len > max_len:
        trace.decimate(decimation_step)
        trace_len = trace.stats.npts
    return trace


if __name__ == "__main__":
    import doctest
    doctest.testmod()
