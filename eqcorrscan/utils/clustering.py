"""
Functions to cluster seismograms by a range of constraints.

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
import matplotlib.pyplot as plt

from multiprocessing import Pool, cpu_count
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from obspy.signal.cross_correlation import xcorr
from obspy import Stream, Catalog

from eqcorrscan.core.match_filter import normxcorr2
from eqcorrscan.utils.mag_calc import dist_calc
from eqcorrscan.utils import stacking


def cross_chan_coherence(st1, st2, allow_shift=False, shift_len=0.2, i=0):
    """
    Calculate cross-channel coherency.

    Determine the cross-channel coherency between two streams of multichannel
    seismic data.

    :type st1: obspy.core.stream.Stream
    :param st1: Stream one
    :type st2: obspy.core.stream.Stream
    :param st2: Stream two
    :type allow_shift: bool
    :param allow_shift:
        Whether to allow the optimum alignment to be found for coherence,
        defaults to `False` for strict coherence
    :type shift_len: int
    :param shift_len: Samples to shift, only used if `allow_shift=True`
    :type i: int
    :param i: index used for parallel async processing, returned unaltered

    :returns:
        cross channel coherence, float - normalized by number of channels,
        and i, where i is int, as input.
    :rtype: tuple
    """
    cccoh = 0.0
    kchan = 0
    if allow_shift:
        for tr in st1:
            tr2 = st2.select(station=tr.stats.station,
                             channel=tr.stats.channel)
            if tr2:
                index, corval = xcorr(tr, tr2[0], shift_len)
                cccoh += corval
                kchan += 1
    else:
        for tr in st1:
            tr1 = tr.data
            # Assume you only have one waveform for each channel
            tr2 = st2.select(station=tr.stats.station,
                             channel=tr.stats.channel)
            if tr2:
                cccoh += normxcorr2(tr1, tr2[0].data)[0][0]
                kchan += 1
    if kchan:
        cccoh /= kchan
        return cccoh, i
    else:
        warnings.warn('No matching channels')
        return 0, i


def distance_matrix(stream_list, allow_shift=False, shift_len=0, cores=1):
    """
    Compute distance matrix for waveforms based on cross-correlations.

    Function to compute the distance matrix for all templates - will give
    distance as 1-abs(cccoh), e.g. a well correlated pair of templates will
    have small distances, and an equally well correlated reverse image will
    have the same distance as a positively correlated image - this is an issue.

    :type stream_list: list
    :param stream_list:
        List of the :class:`obspy.core.stream.Stream`s to compute the distance
        matrix for
    :type allow_shift: bool
    :param allow_shift: To allow templates to shift or not?
    :type shift_len: int
    :param shift_len: How many samples for templates to shift in time
    :type cores: int
    :param cores: Number of cores to parallel process using, defaults to 1.

    :returns: distance matrix
    :rtype: :class:`numpy.ndarray`

    .. warning::
        Because distance is given as :math:`1-abs(coherence)`, negatively
        correlated and positively correlated objects are given the same
        distance.
    """
    # Initialize square matrix
    dist_mat = np.array([np.array([0.0] * len(stream_list))] *
                        len(stream_list))
    for i, master in enumerate(stream_list):
        # Start a parallel processing pool
        pool = Pool(processes=cores)
        # Parallel processing
        results = [pool.apply_async(cross_chan_coherence, args=(master,
                                                                stream_list[j],
                                                                allow_shift,
                                                                shift_len,
                                                                j))
                   for j in range(len(stream_list))]
        pool.close()
        # Extract the results when they are done
        dist_list = [p.get() for p in results]
        # Close and join all the processes back to the master process
        pool.join()
        # Sort the results by the input j
        dist_list.sort(key=lambda tup: tup[1])
        # Sort the list into the dist_mat structure
        for j in range(i, len(stream_list)):
            if i == j:
                dist_mat[i, j] = 0.0
            else:
                dist_mat[i, j] = 1 - dist_list[j][0]
    # Reshape the distance matrix
    for i in range(1, len(stream_list)):
        for j in range(i):
            dist_mat[i, j] = dist_mat.T[i, j]
    return dist_mat


def cluster(template_list, show=True, corr_thresh=0.3, allow_shift=False,
            shift_len=0, save_corrmat=False,
            cores='all', debug=1):
    """
    Cluster template waveforms based on average correlations.

    Function to take a set of templates and cluster them, will return groups
    as lists of streams.  Clustering is done by computing the cross-channel
    correlation sum of each stream in stream_list with every other stream in
    the list.  :mod:`scipy.cluster.hierarchy` functions are then used to
    compute the complete distance matrix, where distance is 1 minus the
    normalised cross-correlation sum such that larger distances are less
    similar events.  Groups are then created by clustering the distance matrix
    at distances less than 1 - corr_thresh.

    Will compute the distance matrix in parallel, using all available cores

    :type template_list: list
    :param template_list:
        List of tuples of the template (:class:`obspy.core.stream.Stream`)
        and the template id to compute clustering for
    :type show: bool
    :param show: plot linkage on screen if True, defaults to True
    :type corr_thresh: float
    :param corr_thresh: Cross-channel correlation threshold for grouping
    :type allow_shift: bool
    :param allow_shift:
        Whether to allow the templates to shift when correlating
    :type shift_len: int
    :param shift_len: How many samples to allow the templates to shift in time
    :type save_corrmat: bool
    :param save_corrmat: If True will save the distance matrix to \
        dist_mat.npy in the local directory.
    :type cores: int
    :param cores: number of cores to use when computing the distance matrix, \
        defaults to 'all' which will work out how many cpus are available \
        and hog them.
    :type debug: int
    :param debug: Level of debugging from 1-5, higher is more output, \
        currently only level 1 implemented.

    :returns:
        List of groups. Each group is a list of
        :class:`obspy.core.stream.Stream`s making up that group.
    """
    if cores == 'all':
        num_cores = cpu_count()
    else:
        num_cores = cores
    # Extract only the Streams from stream_list
    stream_list = [x[0] for x in template_list]
    # Compute the distance matrix
    if debug >= 1:
        print('Computing the distance matrix using %i cores' % num_cores)
    dist_mat = distance_matrix(stream_list, allow_shift, shift_len,
                               cores=num_cores)
    if save_corrmat:
        np.save('dist_mat.npy', dist_mat)
        if debug >= 1:
            print('Saved the distance matrix as dist_mat.npy')
    dist_vec = squareform(dist_mat)
    # plt.matshow(dist_mat, aspect='auto', origin='lower',
    #             cmap=pylab.cm.YlGnBu)
    if debug >= 1:
        print('Computing linkage')
    Z = linkage(dist_vec)
    if show:
        if debug >= 1:
            print('Plotting the dendrogram')
        dendrogram(Z, color_threshold=1 - corr_thresh,
                   distance_sort='ascending')
        plt.show()
    # Get the indices of the groups
    if debug >= 1:
        print('Clustering')
    indices = fcluster(Z, t=1 - corr_thresh, criterion='distance')
    # Indices start at 1...
    group_ids = list(set(indices))  # Unique list of group ids
    if debug >= 1:
        msg = ' '.join(['Found', str(len(group_ids)), 'groups'])
        print(msg)
    # Convert to tuple of (group id, stream id)
    indices = [(indices[i], i) for i in range(len(indices))]
    # Sort by group id
    indices.sort(key=lambda tup: tup[0])
    groups = []
    if debug >= 1:
        print('Extracting and grouping')
    for group_id in group_ids:
        group = []
        for ind in indices:
            if ind[0] == group_id:
                group.append(template_list[ind[1]])
            elif ind[0] > group_id:
                # Because we have sorted by group id, when the index is greater
                # than the group_id we can break the inner loop.
                # Patch applied by CJC 05/11/2015
                groups.append(group)
                break
    # Catch the final group
    groups.append(group)
    return groups


def group_delays(stream_list):
    """
    Group template waveforms according to their arrival times (delays).

    :type stream_list: list
    :param stream_list:
        List of :class:`obspy.core.stream.Stream` waveforms you want to group.

    :returns:
        list of List of :class:`obspy.core.stream.Stream` where each initial
        list is a group with the same delays.
    :rtype: list
    """
    groups = []
    group_delays = []
    group_chans = []
    # Sort templates by number of channels
    stream_list = [(st, len(st)) for st in stream_list]
    stream_list.sort(key=lambda tup: tup[1])
    stream_list = [st[0] for st in stream_list]
    for i, st in enumerate(stream_list):
        msg = ' '.join(['Working on waveform', str(i), 'of',
                        str(len(stream_list))])
        print(msg)
        # Calculate the delays
        starttimes = []
        chans = []
        for tr in st:
            starttimes.append(tr.stats.starttime)
            chans.append((tr.stats.station, tr.stats.channel))
        # This delays calculation will be an issue if we
        # have changes in channels
        delays = [starttimes[m] - min(starttimes)
                  for m in range(len(starttimes))]
        delays = [round(d, 2) for d in delays]
        if len(groups) == 0:
            groups.append([stream_list[i]])
            group_delays.append(delays)
            group_chans.append(chans)
        else:
            j = 0
            match = False
            while not match:
                kmatch = 0
                # Find the set of shared stations and channels
                shared_chans = []
                shared_delays_slave = []
                shared_delays_master = []
                for k, chan in enumerate(chans):
                    if chan in group_chans[j]:
                        shared_chans.append(chan)
                        shared_delays_slave.append(delays[k])
                        shared_delays_master.\
                            append(group_delays[j][group_chans[j].index(chan)])
                # Normalize master and slave delay times
                shared_delays_slave = [delay - min(shared_delays_slave)
                                       for delay in shared_delays_slave]
                shared_delays_master = [delay - min(shared_delays_master)
                                        for delay in shared_delays_master]
                for k in range(len(shared_chans)):
                    # Check if the channel and delay match another group
                    if shared_delays_slave[k] == shared_delays_master[k]:
                        kmatch += 1  # increase the match index
                if kmatch == len(shared_chans):
                    # If all the channels match, add it to the group
                    groups[j].append(stream_list[i])
                    match = True
                elif j < len(groups) - 1:
                    j += 1
                else:
                    # Create a new group and break the loop
                    groups.append([stream_list[i]])
                    group_delays.append(delays)
                    group_chans.append(chans)
                    match = True  # Use this to break the loop
    return groups


def SVD(stream_list, full=False):
    """
    Depreciated. Use svd.
    """
    warnings.warn('Depreciated, use svd instead.')
    return svd(stream_list=stream_list, full=full)


def svd(stream_list, full=False):
    """
    Compute the SVD of a number of templates.

    Returns the singular vectors and singular values of the templates.

    :type stream_list: List of :class: obspy.Stream
    :param stream_list: List of the templates to be analysed
    :type full: bool
    :param full: Whether to compute the full input vector matrix or not.

    :return: SValues(list) for each channel, SVectors(list of ndarray),  \
        UVectors(list of ndarray) for each channel, \
        stachans, List of String (station.channel)

    .. note:: We recommend that you align the data before computing the \
        SVD, e.g., the P-arrival on all templates for the same channel \
        should appear at the same time in the trace.  See the \
        stacking.align_traces function for a way to do this.

    .. note:: Uses the numpy.linalg.svd function, their U, s and V are mapped \
        to UVectors, SValues and SVectors respectively.  Their V (and ours) \
        corresponds to V.H.
    """
    # Convert templates into ndarrays for each channel
    # First find all unique channels:
    stachans = []
    for st in stream_list:
        for tr in st:
            stachans.append('.'.join([tr.stats.station, tr.stats.channel]))
    stachans = list(set(stachans))
    stachans.sort()
    # Initialize a list for the output matrices, one matrix per-channel
    svalues = []
    svectors = []
    uvectors = []
    for stachan in stachans:
        lengths = []
        for st in stream_list:
            tr = st.select(station=stachan.split('.')[0],
                           channel=stachan.split('.')[1])
            if len(tr) > 0:
                tr = tr[0]
            else:
                warnings.warn('Stream does not contain ' + stachan)
                continue
            lengths.append(len(tr.data))
        min_length = min(lengths)
        for stream in stream_list:
            chan = stream.select(station=stachan.split('.')[0],
                                 channel=stachan.split('.')[1])
            if chan:
                if len(chan[0].data) > min_length:
                    if abs(len(chan[0].data) - min_length) > 0.1 *\
                            chan[0].stats.sampling_rate:
                        raise IndexError('More than 0.1 s length '
                                         'difference, align and fix')
                    warnings.warn('Channels are not equal length, trimming')
                    chan[0].data = chan[0].data[0:min_length]
                if 'chan_mat' not in locals():
                    chan_mat = chan[0].data
                else:
                    chan_mat = np.vstack((chan_mat, chan[0].data))
        if not len(chan_mat.shape) > 1:
            warnings.warn('Matrix of traces is less than 2D for %s' % stachan)
            continue
        chan_mat = np.asarray(chan_mat)
        u, s, v = np.linalg.svd(chan_mat, full_matrices=full)
        svalues.append(s)
        svectors.append(v)
        uvectors.append(u)
        del(chan_mat)
    return svectors, svalues, uvectors, stachans


def empirical_SVD(stream_list, linear=True):
    """
    Empirical subspace detector generation function.

    Takes a list of \
    templates and computes the stack as the first order subspace detector, \
    and the differential of this as the second order subspace detector \
    following the empirical subspace method of
    `Barrett & Beroza, 2014 - SRL
    <http://srl.geoscienceworld.org/content/85/3/594.extract>`_.

    :type stream_list: list
    :param stream_list:
        list of streams to compute the subspace detectors from, where streams
        are :class:`obspy.core.stream.Stream` objects.
    :type linear: bool
    :param linear: Set to true by default to compute the linear stack as the \
        first subspace vector, False will use the phase-weighted stack as the \
        first subspace vector.

    :returns: list of two :class:`obspy.core.stream.Stream` s
    """
    from eqcorrscan.utils import stacking
    # Run a check to ensure all traces are the same length
    stachans = list(set([(tr.stats.station, tr.stats.channel)
                         for st in stream_list for tr in st]))
    for stachan in stachans:
        lengths = []
        for st in stream_list:
            lengths.append(len(st.select(station=stachan[0],
                                         channel=stachan[1])[0]))
        min_length = min(lengths)
        for st in stream_list:
            tr = st.select(station=stachan[0],
                           channel=stachan[1])[0]
            if len(tr.data) > min_length:
                if abs(len(tr.data) - min_length) > (0.1 *
                                                     tr.stats.sampling_rate):
                        raise IndexError('More than 0.1 s length '
                                         'difference, align and fix')
                warnings.warn(str(tr) + ' is not the same length as others, ' +
                              'trimming the end')
                tr.data = tr.data[0:min_length]
    if linear:
        first_subspace = stacking.linstack(stream_list)
    else:
        first_subspace = stacking.PWS_stack(streams=stream_list)
    second_subspace = first_subspace.copy()
    for i in range(len(second_subspace)):
        second_subspace[i].data = np.diff(second_subspace[i].data)
        second_subspace[i].stats.starttime += 0.5 * \
            second_subspace[i].stats.delta
    return [first_subspace, second_subspace]


def SVD_2_stream(SVectors, stachans, k, sampling_rate):
    """
    Convert the singular vectors output by SVD to streams.

    One stream will be generated for each singular vector level,
    for all channels.  Useful for plotting, and aiding seismologists thinking
    of waveforms!

    :type SVectors: list
    :param SVectors: List of :class:`numpy.ndarray` Singular vectors
    :type stachans: list
    :param stachans: List of station.channel Strings
    :type k: int
    :param k: Number of streams to return = number of SV's to include
    :type sampling_rate: float
    :param sampling_rate: Sampling rate in Hz

    :returns:
        SVstreams, List of :class:`obspy.core.stream.Stream`, with
        SVStreams[0] being composed of the highest rank singular vectors.
    """
    from obspy import Stream, Trace
    SVstreams = []
    for i in range(k):
        SVstream = []
        for j, stachan in enumerate(stachans):
            if len(SVectors[j]) <= k:
                warnings.warn('Too few traces at %s for a %02d dimensional '
                              'subspace. Detector streams will not include '
                              'this channel.' % (stachan, k))
            else:
                SVstream.append(Trace(SVectors[j][i],
                                      header={'station': stachan.split('.')[0],
                                              'channel': stachan.split('.')[1],
                                              'sampling_rate': sampling_rate}))
        SVstreams.append(Stream(SVstream))
    return SVstreams


def corr_cluster(trace_list, thresh=0.9):
    """
    Group traces based on correlations above threshold with the stack.

    Will run twice, once with a lower threshold to remove large outliers that
    would negatively affect the stack, then again with your threshold.

    :type trace_list: list
    :param trace_list:
        List of :class:`obspy.core.stream.Trace`s to compute similarity between
    :type thresh: float
    :param thresh: Correlation threshold between -1-1

    :returns:
        :class:`numpy.ndarray` of bool of whether that trace correlates well
        enough (above your given threshold) with the stack.

    .. note::
        We recommend that you align the data before computing the clustering,
        e.g., the P-arrival on all templates for the same channel should
        appear at the same time in the trace.  See the
        :func:`eqcorrscan.utils.stacking.align_traces` function for a way to do
        this.
    """
    stack = stacking.linstack([Stream(tr) for tr in trace_list])[0]
    output = np.array([False] * len(trace_list))
    group1 = []
    for i, tr in enumerate(trace_list):
        if normxcorr2(tr.data, stack.data)[0][0] > 0.6:
            output[i] = True
            group1.append(tr)
    if not group1:
        warnings.warn('Nothing made it past the first 0.6 threshold')
        return output
    stack = stacking.linstack([Stream(tr) for tr in group1])[0]
    group2 = []
    for i, tr in enumerate(trace_list):
        if normxcorr2(tr.data, stack.data)[0][0] > thresh:
            group2.append(tr)
            output[i] = True
        else:
            output[i] = False
    return output


def extract_detections(detections, templates, archive, arc_type,
                       extract_len=90.0, outdir=None, extract_Z=True,
                       additional_stations=[]):
    """
    Extract waveforms associated with detections

    Takes a list of detections for the template, template.  Waveforms will be
    returned as a list of :class:`obspy.core.stream.Stream` containing
    segments of extract_len.  They will also be saved if outdir is set.
    The default is unset.  The  default extract_len is 90 seconds per channel.

    :type detections: list
    :param detections: List of :class:`eqcorrscan.core.match_filter.DETECTION`.
    :type templates: list
    :param templates:
        A list of tuples of the template name and the template Stream used
        to detect detections.
    :type archive: str
    :param archive:
        Either name of archive or path to continuous data, see
        :func:`eqcorrscan.utils.archive_read` for details
    :type arc_type: str
    :param arc_type: Type of archive, either seishub, FDSN, day_vols
    :type extract_len: float
    :param extract_len:
        Length to extract around the detection (will be equally cut around
        the detection time) in seconds.  Default is 90.0.
    :type outdir: str
    :param outdir:
        Default is None, with None set, no files will be saved,
        if set each detection will be saved into this directory with files
        named according to the detection time, NOT than the waveform
        start time. Detections will be saved into template subdirectories.
    :type extract_Z: bool
    :param extract_Z:
        Set to True to also extract Z channels for detections delays will be
        the same as horizontal channels, only applies if only horizontal
        channels were used in the template.
    :type additional_stations: list
    :param additional_stations:
        List of tuples of (station, channel, network) to also extract data
        for using an average delay.

    :returns: list of :class:`obspy.core.streams.Stream`
    :rtype: list

    .. rubric: Example

    >>> from eqcorrscan.utils.clustering import extract_detections
    >>> from eqcorrscan.core.match_filter import DETECTION
    >>> from obspy import read, UTCDateTime
    >>> import os
    >>> # Use some dummy detections, you would use real one
    >>> detections = [DETECTION('temp1', UTCDateTime(2012, 3, 26, 9, 15), 2,
    ...                         ['WHYM', 'EORO'], 2, 1.2, 'corr'),
    ...               DETECTION('temp2',UTCDateTime(2012, 3, 26, 18, 5), 2,
    ...                         ['WHYM', 'EORO'], 2, 1.2, 'corr')]
    >>> path_to_templates = os.path.join('eqcorrscan', 'tests', 'test_data')
    >>> archive = os.path.join(path_to_templates, 'day_vols')
    >>> template_files = [os.path.join(path_to_templates, 'temp1.ms'),
    ...                   os.path.join(path_to_templates, 'temp2.ms')]
    >>> templates = [('temp' + str(i), read(filename))
    ...              for i, filename in enumerate(template_files)]
    >>> extracted = extract_detections(detections, templates,
    ...                                archive=archive, arc_type='day_vols')
    Working on detections for day: 2012-03-26T00:00:00.000000Z
    Cutting for detections at: 2012/03/26 09:15:00
    Cutting for detections at: 2012/03/26 18:05:00
    >>> print(extracted[0].sort())
    2 Trace(s) in Stream:
    AF.EORO..SHZ | 2012-03-26T09:14:15.000000Z - 2012-03-26T09:15:45.000000Z |\
 1.0 Hz, 91 samples
    AF.WHYM..SHZ | 2012-03-26T09:14:15.000000Z - 2012-03-26T09:15:45.000000Z |\
 1.0 Hz, 91 samples
    >>> print(extracted[1].sort())
    2 Trace(s) in Stream:
    AF.EORO..SHZ | 2012-03-26T18:04:15.000000Z - 2012-03-26T18:05:45.000000Z |\
 1.0 Hz, 91 samples
    AF.WHYM..SHZ | 2012-03-26T18:04:15.000000Z - 2012-03-26T18:05:45.000000Z |\
 1.0 Hz, 91 samples
    """
    from obspy import UTCDateTime
    import os
    from eqcorrscan.utils.archive_read import read_data
    # Sort the template according to start-times, needed so that stachan[i]
    # corresponds to delays[i]
    all_delays = []  # List of tuples of template name, delays
    all_stachans = []
    for template in templates:
        templatestream = template[1].sort(['starttime'])
        stachans = [(tr.stats.station, tr.stats.channel)
                    for tr in templatestream]
        mintime = templatestream[0].stats.starttime
        delays = [tr.stats.starttime - mintime for tr in templatestream]
        all_delays.append((template[0], delays))
        all_stachans.append((template[0], stachans))
    # Sort the detections and group by day
    detections.sort(key=lambda d: d.detect_time)
    detection_days = [detection.detect_time.date
                      for detection in detections]
    detection_days = list(set(detection_days))
    detection_days.sort()
    detection_days = [UTCDateTime(d) for d in detection_days]

    # Initialize output list
    detection_wavefiles = []

    # Also include Z channels when extracting detections
    if extract_Z:
        new_all_stachans = []
        new_all_delays = []
        for t, template in enumerate(all_stachans):
            stachans = template[1]
            delays = all_delays[t][1]
            new_stachans = []
            new_delays = []
            j = 0
            for i, stachan in enumerate(stachans):
                if j == 1:
                    new_stachans.append((stachan[0], stachan[1][0] + 'Z'))
                    new_delays.append(delays[i])
                    new_stachans.append(stachan)
                    new_delays.append(delays[i])
                    j = 0
                else:
                    new_stachans.append(stachan)
                    new_delays.append(delays[i])
                    j += 1
            new_all_stachans.append((template[0], new_stachans))
            new_all_delays.append((template[0], new_delays))
        all_delays = new_all_delays
        all_stachans = new_all_stachans
    if not len(additional_stations) == 0:
        print('Adding additional stations')
        for t, template in enumerate(all_stachans):
            av_delay = np.mean(all_delays[t][1])
            for sta in additional_stations:
                if sta not in template[1]:
                    print('Added station ' + '.'.join(sta))
                    template[1].append(sta)
                    all_delays[t][1].append(av_delay)
    del stachans
    # Loop through the days
    for detection_day in detection_days:
        print('Working on detections for day: ' + str(detection_day))
        stachans = list(set([stachans[1] for stachans in all_stachans][0]))
        # List of all unique stachans - read in all data
        st = read_data(archive=archive, arc_type=arc_type, day=detection_day,
                       stachans=stachans)
        st.merge(fill_value='interpolate')
        day_detections = [detection for detection in detections
                          if UTCDateTime(detection.detect_time.date) ==
                          detection_day]
        del stachans, delays
        for detection in day_detections:
            template = [t[1] for t in templates
                        if t[0] == detection.template_name]
            print('Cutting for detections at: ' +
                  detection.detect_time.strftime('%Y/%m/%d %H:%M:%S'))
            detect_wav = st.copy()
            for tr in detect_wav:
                tr.trim(starttime=UTCDateTime(detection.detect_time) -
                        extract_len / 2,
                        endtime=UTCDateTime(detection.detect_time) +
                        extract_len / 2)
            if outdir:
                if not os.path.isdir(os.path.join(outdir, template)):
                    os.makedirs(os.path.join(outdir, template))
                detect_wav.write(os.path.join(outdir, template,
                                              detection.detect_time.
                                              strftime('%Y-%m-%d_%H-%M-%S') +
                                              '.ms'),
                                 format='MSEED', encoding='STEIM2')
                print('Written file: %s' % os.path.join(outdir, template,
                                                        detection.detect_time.
                                                        strftime('%Y-%m-%d'
                                                                 '_%H-%M-%S') +
                                                        '.ms'))
            if not outdir:
                detection_wavefiles.append(detect_wav)
            del detect_wav
        del st
        if outdir:
            detection_wavefiles = []
    if not outdir:
        return detection_wavefiles
    else:
        return


def dist_mat_km(catalog):
    """
    Compute the distance matrix for all a catalog using epicentral separation.

    Will give physical distance in kilometers.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog for which to compute the distance matrix

    :returns: distance matrix
    :rtype: :class:`numpy.ndarray`
    """
    # Initialize square matrix
    dist_mat = np.array([np.array([0.0] * len(catalog))] *
                        len(catalog))
    # Calculate distance vector for each event
    for i, master in enumerate(catalog):
        mast_list = []
        if master.preferred_origin():
            master_ori = master.preferred_origin()
        else:
            master_ori = master.origins[0]
        master_tup = (master_ori.latitude,
                      master_ori.longitude,
                      master_ori.depth // 1000)
        for slave in catalog:
            if master.preferred_origin():
                slave_ori = slave.preferred_origin()
            else:
                slave_ori = slave.origins[0]
            slave_tup = (slave_ori.latitude,
                         slave_ori.longitude,
                         slave_ori.depth // 1000)
            mast_list.append(dist_calc(master_tup, slave_tup))
        # Sort the list into the dist_mat structure
        for j in range(i, len(catalog)):
            dist_mat[i, j] = mast_list[j]
    # Reshape the distance matrix
    for i in range(1, len(catalog)):
        for j in range(i):
            dist_mat[i, j] = dist_mat.T[i, j]
    return dist_mat


def space_cluster(catalog, d_thresh, show=True):
    """
    Cluster a catalog by distance only.

    Will compute the matrix of physical distances between events and utilize
    the :mod:`scipy.clustering.hierarchy` module to perform the clustering.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog of events to clustered
    :type d_thresh: float
    :param d_thresh: Maximum inter-event distance threshold

    :returns: list of :class:`obspy.core.event.Catalog` objects
    :rtype: list

    >>> from eqcorrscan.utils.clustering import space_cluster
    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> client = Client("NCEDC")
    >>> starttime = UTCDateTime("2002-01-01")
    >>> endtime = UTCDateTime("2002-02-01")
    >>> cat = client.get_events(starttime=starttime, endtime=endtime,
    ...                        minmagnitude=2)
    >>> groups = space_cluster(catalog=cat, d_thresh=2, show=False)

    >>> from eqcorrscan.utils.clustering import space_cluster
    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> client = Client("IRIS")
    >>> starttime = UTCDateTime("2002-01-01")
    >>> endtime = UTCDateTime("2002-02-01")
    >>> cat = client.get_events(starttime=starttime, endtime=endtime,
    ...                        minmagnitude=6, catalog="ISC")
    >>> groups = space_cluster(catalog=cat, d_thresh=1000, show=False)
    """
    # Compute the distance matrix and linkage
    dist_mat = dist_mat_km(catalog)
    dist_vec = squareform(dist_mat)
    Z = linkage(dist_vec, method='average')

    # Cluster the linkage using the given threshold as the cutoff
    indices = fcluster(Z, t=d_thresh, criterion='distance')
    group_ids = list(set(indices))
    indices = [(indices[i], i) for i in range(len(indices))]

    if show:
        # Plot the dendrogram...if it's not way too huge
        dendrogram(Z, color_threshold=d_thresh,
                   distance_sort='ascending')
        plt.show()

    # Sort by group id
    indices.sort(key=lambda tup: tup[0])
    groups = []
    for group_id in group_ids:
        group = Catalog()
        for ind in indices:
            if ind[0] == group_id:
                group.append(catalog[ind[1]])
            elif ind[0] > group_id:
                # Because we have sorted by group id, when the index is greater
                # than the group_id we can break the inner loop.
                # Patch applied by CJC 05/11/2015
                groups.append(group)
                break
    groups.append(group)
    return groups


def space_time_cluster(catalog, t_thresh, d_thresh):
    """
    Cluster detections in space and time.

    Use to separate repeaters from other events.  Clusters by distance
    first, then removes events in those groups that are at different times.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog of events to clustered
    :type t_thresh: float
    :param t_thresh: Maximum inter-event time threshold in seconds
    :type d_thresh: float
    :param d_thresh: Maximum inter-event distance in km

    :returns: list of :class:`obspy.core.event.Catalog` objects
    :rtype: list

    >>> from eqcorrscan.utils.clustering import space_time_cluster
    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> client = Client("IRIS")
    >>> starttime = UTCDateTime("2002-01-01")
    >>> endtime = UTCDateTime("2002-02-01")
    >>> cat = client.get_events(starttime=starttime, endtime=endtime,
    ...                         minmagnitude=6, catalog="ISC")
    >>> groups = space_time_cluster(catalog=cat, t_thresh=86400, d_thresh=1000)
    """
    initial_spatial_groups = space_cluster(catalog=catalog, d_thresh=d_thresh,
                                           show=False)
    # Need initial_spatial_groups to be lists at the moment
    initial_spatial_lists = []
    for group in initial_spatial_groups:
        initial_spatial_lists.append(list(group))
    # Check within these groups and throw them out if they are not close in
    # time.
    groups = []
    for group in initial_spatial_lists:
        for master in group:
            for event in group:
                if abs(event.preferred_origin().time -
                       master.preferred_origin().time) > t_thresh:
                    # If greater then just put event in on it's own
                    groups.append([event])
                    group.remove(event)
        groups.append(group)
    return [Catalog(group) for group in groups]


def re_thresh_csv(path, old_thresh, new_thresh, chan_thresh):
    """
    Remove detections by changing the threshold.

    Can only be done to remove detection by increasing threshold,
    threshold lowering will have no effect.

    :type path: str
    :param path: Path to the .csv detection file
    :type old_thresh: float
    :param old_thresh: Old threshold MAD multiplier
    :type new_thresh: float
    :param new_thresh: New threshold MAD multiplier
    :type chan_thresh: int
    :param chan_thresh: Minimum number of channels for a detection

    :returns: List of detections
    :rtype: list

    .. rubric:: Example

    >>> from eqcorrscan.utils.clustering import re_thresh_csv
    >>> import os
    >>> det_file = os.path.join('eqcorrscan', 'tests', 'test_data',
    ...                         'expected_tutorial_detections.txt')
    >>> detections = re_thresh_csv(path=det_file, old_thresh=8, new_thresh=10,
    ...                            chan_thresh=3)
    Read in 22 detections
    Left with 17 detections
    """
    from eqcorrscan.core.match_filter import read_detections
    old_detections = read_detections(path)
    old_thresh = float(old_thresh)
    new_thresh = float(new_thresh)
    # Be nice, ensure that the thresholds are float
    detections = []
    detections_in = 0
    detections_out = 0
    for detection in old_detections:
        detections_in += 1
        if abs(detection.detect_val) >=\
           (new_thresh / old_thresh) * detection.threshold and\
           detection.no_chans >= chan_thresh:
            detections_out += 1
            detections.append(detection)
    print('Read in %i detections' % detections_in)
    print('Left with %i detections' % detections_out)
    return detections


if __name__ == "__main__":
    import doctest
    doctest.testmod()
