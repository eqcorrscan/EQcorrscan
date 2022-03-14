"""
Functions to cluster seismograms by a range of constraints.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import os
import logging
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, Catalog, UTCDateTime, Trace
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

from eqcorrscan.utils import stacking
from eqcorrscan.utils.archive_read import read_data
from eqcorrscan.utils.correlate import (
    get_array_xcorr, get_stream_xcorr, CorrelationError)
from eqcorrscan.utils.pre_processing import _prep_data_for_correlation

Logger = logging.getLogger(__name__)


def cross_chan_correlation(
        st1, streams, shift_len=0.0, allow_individual_trace_shifts=True,
        xcorr_func='fftw', concurrency="concurrent", cores=1, **kwargs):
    """
    Calculate cross-channel correlation.

    Determine the cross-channel correlation between two streams of
    multichannel seismic data.

    :type st1: obspy.core.stream.Stream
    :param st1: Stream one
    :type streams: list
    :param streams: Streams to compare to.
    :type shift_len: float
    :param shift_len: How many seconds for templates to shift
    :type allow_individual_trace_shifts: bool
    :param allow_individual_trace_shifts:
        Controls whether templates are shifted by shift_len in relation to the
        picks as a whole, or whether each trace can be shifted individually.
        Defaults to True.
    :type xcorr_func: str, callable
    :param xcorr_func:
        The method for performing correlations. Accepts either a string or
        callable. See :func:`eqcorrscan.utils.correlate.register_array_xcorr`
        for more details
    :type concurrency: str
    :param concurrency: Concurrency for xcorr-func.
    :type cores: int
    :param cores: Number of threads to parallel over

    :returns:
        cross channel correlation, float - normalized by number of channels.
        locations of maximums
    :rtype: numpy.ndarray, numpy.ndarray

    .. Note::
        If no matching channels were found then the coherance and index for
        that stream will be nan.
    """
    # Cut all channels in stream-list to be the correct length (shorter than
    # st1 if stack = False by shift_len).
    allow_individual_trace_shifts = (
        allow_individual_trace_shifts and shift_len > 0)
    n_streams = len(streams)
    df = st1[0].stats.sampling_rate
    end_trim = int((shift_len * df) / 2)
    _streams = []
    if end_trim > 0:
        for stream in streams:
            _stream = stream.copy()  # Do not work on the users data
            for tr in _stream:
                tr.data = tr.data[end_trim: -end_trim]
                if tr.stats.sampling_rate != df:
                    raise NotImplementedError("Sampling rates differ")
            _streams.append(_stream)
        streams = _streams
    else:
        # _prep_data_for_correlation works in place on data.
        # We need to copy it first.
        streams = [stream.copy() for stream in streams]
    # Check which channels are in st1 and match those in the stream_list
    st1_preped, prep_streams, stream_indexes = _prep_data_for_correlation(
        stream=st1.copy(), templates=streams,
        template_names=list(range(len(streams))), force_stream_epoch=False)
    # Run the correlations
    multichannel_normxcorr = get_stream_xcorr(xcorr_func, concurrency)
    [cccsums, no_chans, _] = multichannel_normxcorr(
        templates=prep_streams, stream=st1_preped, cores=cores, stack=False,
        **kwargs)
    # Find maximas, sum and divide by no_chans
    if allow_individual_trace_shifts:
        coherances = cccsums.max(axis=-1).sum(axis=-1) / no_chans
    else:
        cccsums = cccsums.sum(axis=1)
        coherances = cccsums.max(axis=-1) / no_chans
    # Subtract half length of correlogram and convert positions to seconds
    positions = (cccsums.argmax(axis=-1) - end_trim) / df

    # This section re-orders the coherences to correspond to the order of the
    # input streams
    _coherances = np.empty(n_streams)
    if allow_individual_trace_shifts:
        n_max_traces = max([len(st) for st in prep_streams])
        # Set shifts for nan-traces to nan
        for i, tr in enumerate(st1_preped):
            if np.ma.is_masked(tr.data):
                positions[:, i] = np.nan
    else:
        positions = positions[:, np.newaxis]
        n_max_traces = 1
    n_shifts_per_stream = positions.shape[1]
    _positions = np.empty([n_streams, n_max_traces])

    _coherances.fill(np.nan)
    _positions.fill(np.nan)
    # Insert the correlations and shifts at the correct index for the templates
    _coherances[np.ix_(stream_indexes)] = coherances
    _positions[np.ix_(stream_indexes, range(n_shifts_per_stream))] = (
        positions)

    if not allow_individual_trace_shifts:  # remove empty third axis from array
        _positions = _positions[:, ]
    return _coherances, _positions


def handle_distmat_nans(dist_mat, replace_nan_distances_with=None):
    """
    Checks for nans and optionally fills missing correlations (nans) with a
    replacement value.

    :type dist_mat: np.ndarray
    :param dist_mag: Distance matrix.
    :type replace_nan_distances_with: None, 'mean', 'min', or float
    :param replace_nan_distances_with:
        Controls how the clustering handles nan-distances in the distance
        matrix. None/False only performs a check, while other choices (e.g.,
        1, 'mean', 'min' or float) replace nans in the distance matrix.

    :returns:
        Distance matrix.
        :type: np.ndarray
    """
    missing_corrs = ~np.isfinite(dist_mat)
    missing_vals = np.full_like(dist_mat, fill_value=1)
    if not replace_nan_distances_with:
        if np.isnan(dist_mat).any():
            raise CorrelationError(
                "The distance matrix contains nans, which may indicate that " +
                "some templates have no matching traces and can only be " +
                "indirectly linked via other templates. You should check " +
                " templates. Then you can try cluster() with "
                "replace_nan_distances_with=1")
    elif isinstance(replace_nan_distances_with, (int, float)):
        missing_vals = np.full_like(
            dist_mat, fill_value=replace_nan_distances_with)
    elif replace_nan_distances_with == 'mean':
        col_mean = np.nanmean(dist_mat, 0, keepdims=1)
        row_mean = np.nanmean(dist_mat, 1, keepdims=1)
        missing_vals = ((np.repeat(col_mean, len(row_mean), 0) +
                        np.repeat(row_mean, len(col_mean), 1)) / 2)
    elif replace_nan_distances_with == 'min':
        col_min = np.nanmin(dist_mat, 0, keepdims=1)
        row_min = np.nanmin(dist_mat, 1, keepdims=1)
        missing_vals = np.minimum(((np.repeat(col_min, len(row_min), 0),
                                    np.repeat(row_min, len(col_min), 1))))
    else:
        raise NotImplementedError(
            'replace_nan_distances_with={} is not supported'.format(
                replace_nan_distances_with))
    dist_mat = np.where(missing_corrs, missing_vals, dist_mat)
    assert np.allclose(dist_mat, dist_mat.T, atol=0.00001)
    # force perfect symmetry
    dist_mat = (dist_mat + dist_mat.T) / 2
    return dist_mat


def distance_matrix(stream_list, shift_len=0.0,
                    replace_nan_distances_with=None,
                    allow_individual_trace_shifts=True, cores=1):
    """
    Compute distance matrix for waveforms based on cross-correlations.

    Function to compute the distance matrix for all templates - will give
    distance as 1-abs(cccoh), e.g. a well correlated pair of templates will
    have small distances, and an equally well correlated reverse image will
    have the same distance as a positively correlated image - this is an issue.

    :type stream_list: list
    :param stream_list:
        List of the :class:`obspy.core.stream.Stream` to compute the distance
        matrix for
    :type shift_len: float
    :param shift_len: How many seconds for templates to shift
    :type allow_individual_trace_shifts: bool
    :param allow_individual_trace_shifts:
        Controls whether templates are shifted by shift_len in relation to the
        picks as a whole, or whether each trace can be shifted individually.
        Defaults to True.
    :type replace_nan_distances_with: None, 'mean', 'min', or float
    :param replace_nan_distances_with:
        Controls how the clustering handles nan-distances in the distance
        matrix. None/False only performs a check, while other choices (e.g.,
        1, 'mean', 'min' or float) replace nans in the distance matrix.
    :type cores: int
    :param cores: Number of cores to parallel process using, defaults to 1.

    :returns:
        - distance matrix (:py:class:`numpy.ndarray`) of size
          len(stream_list)**2
        - shift matrix (:py:class:`numpy.ndarray`) containing shifts between
          traces of the sorted streams. Size is len(stream_list)**2 * x, where
          x is 1 for shift_len=0 and/or allow_individual_trace_shifts=False.
          Missing correlations are indicated by nans.
        - shift dict (:py:class:`dict`):
          dictionary of (template_id: trace_dict) where trace_dict contains
          (trace.id: shift matrix (size `len(stream_list)**2`) for trace.id)

    .. warning::
        Because distance is given as :math:`1-abs(coherence)`, negatively
        correlated and positively correlated objects are given the same
        distance.

    .. note::
        Requires all traces to have the same sampling rate and same length.
    """
    n_streams = len(stream_list)
    # May have to allow duplicate channels for P- and S-picks at each station
    stream_list = [st.sort() for st in stream_list]
    uniq_traces = set([tr.id for st in stream_list for tr in st])
    n_uniq_traces = len(uniq_traces)
    # Initialize square matrix
    dist_mat = np.zeros([n_streams, n_streams])
    shift_mat = np.empty([n_streams, n_streams, n_uniq_traces])
    shift_mat[:] = np.nan
    shift_dict = dict()
    for i, master in enumerate(stream_list):
        dist_list, shift_list = cross_chan_correlation(
            st1=master, streams=stream_list, shift_len=shift_len,
            allow_individual_trace_shifts=allow_individual_trace_shifts,
            xcorr_func='fftw', cores=cores)
        dist_mat[i] = 1 - dist_list
        master_ids = [tr.id for tr in master]
        master_trace_indcs = [
            j for j, tr_id in enumerate(uniq_traces) if tr_id in master_ids]
        # Sort computed shifts into shift-matrix. shift_list could contain a
        # nan-column that needs to be ignored here (only when earliest trace is
        # missing)
        shift_mat[np.ix_([i], list(range(n_streams)), master_trace_indcs)] = (
            shift_list[:, ~np.all(np.isnan(shift_list), axis=0)])
        # Add trace-id with corresponding shift-matrix to shift-dictionary
        shift_mat_list = [shift_mat[:, :, mti] for mti in master_trace_indcs]
        trace_shift_dict = dict(zip(master_ids, shift_mat_list))
        shift_dict[i] = trace_shift_dict
    if shift_len == 0:
        dist_mat = handle_distmat_nans(
            dist_mat, replace_nan_distances_with=replace_nan_distances_with)
    else:
        # get the shortest distance for each correlation pair
        dist_mat_shortest = np.minimum(dist_mat, dist_mat.T)
        # Indicator says which matrix has shortest dist: value 0: mat2; 1: mat1
        mat_indicator = dist_mat_shortest == dist_mat
        mat_indicator = np.repeat(mat_indicator[:, :, np.newaxis],
                                  n_uniq_traces, axis=2)[:, :]
        # Get shift for the shortest distances
        shift_mat = (
            shift_mat * mat_indicator +
            np.transpose(shift_mat, [1, 0, 2]) * (1 - mat_indicator))
        dist_mat = dist_mat_shortest
    # Squeeze matrix to 2 axis (ignore nans) if 3rd dimension not needed
    if shift_len == 0 or allow_individual_trace_shifts is False:
        shift_mat = np.nanmean(shift_mat, axis=2)
    np.fill_diagonal(dist_mat, 0)
    return dist_mat, shift_mat.squeeze(), shift_dict


def cluster(template_list, show=True, corr_thresh=0.3, shift_len=0,
            allow_individual_trace_shifts=True, save_corrmat=False,
            replace_nan_distances_with=None, cores='all', **kwargs):
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

    When distance_matrix contains NaNs (event pairs that cannot be directly
    compared), then the mean correlation between templates is used instead of
    NaN (see https://github.com/eqcorrscan/EQcorrscan/issues/484).

    Will compute the distance matrix in parallel, using all available cores.
    The method, metric, and order to compute linkage from the distance matrix
    can be controled with parameters from scipy.cluster.hierarchy.linkage as
    kwargs.

    :type template_list: list
    :param template_list:
        List of tuples of the template (:class:`obspy.core.stream.Stream`)
        and the template id to compute clustering for
    :type show: bool
    :param show: plot linkage on screen if True, defaults to True
    :type corr_thresh: float
    :param corr_thresh: Cross-channel correlation threshold for grouping
    :type shift_len: float
    :param shift_len: How many seconds to allow the templates to shift
    :type allow_individual_trace_shifts: bool
    :param allow_individual_trace_shifts:
        Controls whether templates are shifted by shift_len in relation to the
        picks as a whole, or whether each trace can be shifted individually.
        Defaults to True.
    :type save_corrmat: bool
    :param save_corrmat:
        If True will save the distance matrix to dist_mat.npy in the local
        directory.
    :type replace_nan_distances_with: None, 'mean', 'min', or float
    :param replace_nan_distances_with:
        Controls how the clustering handles nan-distances in the distance
        matrix. None/False only performs a check, while other choices (e.g.,
        1, 'mean', 'min' or float) replace nans in the distance matrix.
    :type cores: int
    :param cores:
        number of cores to use when computing the distance matrix, defaults to
        'all' which will work out how many cpus are available and hog them.

    :returns:
        List of groups. Each group is a list of
        :class:`obspy.core.stream.Stream` making up that group.
    """
    if cores == 'all':
        num_cores = cpu_count()
    else:
        num_cores = cores
    # Extract only the Streams from stream_list
    stream_list = [x[0] for x in template_list]
    # Compute the distance matrix
    Logger.info('Computing the distance matrix using %i cores' % num_cores)
    dist_mat, shift_mat, shift_dict = distance_matrix(
        stream_list=stream_list, shift_len=shift_len, cores=num_cores,
        replace_nan_distances_with=replace_nan_distances_with,
        allow_individual_trace_shifts=allow_individual_trace_shifts)
    if save_corrmat:
        np.save('dist_mat.npy', dist_mat)
        Logger.info('Saved the distance matrix as dist_mat.npy')
    dist_mat = handle_distmat_nans(
        dist_mat, replace_nan_distances_with=replace_nan_distances_with)
    dist_vec = squareform(dist_mat)
    Logger.info('Computing linkage')
    Z = linkage(dist_vec, **kwargs)
    if show:
        Logger.info('Plotting the dendrogram')
        dendrogram(Z, color_threshold=1 - corr_thresh,
                   distance_sort='ascending')
        plt.show()
    # Get the indices of the groups
    Logger.info('Clustering')
    indices = fcluster(Z, t=1 - corr_thresh, criterion='distance')
    # Indices start at 1...
    group_ids = list(set(indices))  # Unique list of group ids
    Logger.info(' '.join(['Found', str(len(group_ids)), 'groups']))
    # Convert to tuple of (group id, stream id)
    indices = [(indices[i], i) for i in range(len(indices))]
    # Sort by group id
    indices.sort(key=lambda tup: tup[0])
    groups = []
    Logger.info('Extracting and grouping')
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
        Logger.info(msg)
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
                        shared_delays_master. \
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
    stachans = list(set([(tr.stats.station, tr.stats.channel)
                         for st in stream_list for tr in st]))
    stachans.sort()
    # Initialize a list for the output matrices, one matrix per-channel
    svalues = []
    svectors = []
    uvectors = []
    for stachan in stachans:
        lengths = []
        for st in stream_list:
            tr = st.select(station=stachan[0],
                           channel=stachan[1])
            if len(tr) > 0:
                tr = tr[0]
            else:
                Logger.warning('Stream does not contain %s'
                               % '.'.join(list(stachan)))
                continue
            lengths.append(len(tr.data))
        min_length = min(lengths)
        for stream in stream_list:
            chan = stream.select(station=stachan[0],
                                 channel=stachan[1])
            if chan:
                if len(chan[0].data) > min_length:
                    if abs(len(chan[0].data) - min_length) > 0.1 * \
                            chan[0].stats.sampling_rate:
                        raise IndexError('More than 0.1 s length '
                                         'difference, align and fix')
                    Logger.warning('Channels are not equal length, trimming')
                    chan[0].data = chan[0].data[0:min_length]
                if 'chan_mat' not in locals():
                    chan_mat = chan[0].data
                else:
                    chan_mat = np.vstack((chan_mat, chan[0].data))
        if not len(chan_mat.shape) > 1:
            Logger.warning('Matrix of traces is less than 2D for %s'
                           % '.'.join(list(stachan)))
            continue
        # Be sure to transpose chan_mat as waveforms must define columns
        chan_mat = np.asarray(chan_mat)
        u, s, v = np.linalg.svd(chan_mat.T, full_matrices=full)
        svalues.append(s)
        svectors.append(v)
        uvectors.append(u)
        del (chan_mat)
    return uvectors, svalues, svectors, stachans


def empirical_svd(stream_list, linear=True):
    """
    Empirical subspace detector generation function.

    Takes a list of templates and computes the stack as the first order
    subspace detector, and the differential of this as the second order
    subspace detector following the empirical subspace method of
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
                sr = tr.stats.sampling_rate
                if abs(len(tr.data) - min_length) > (0.1 * sr):
                    msg = 'More than 0.1 s length difference, align and fix'
                    raise IndexError(msg)
                msg = ' is not the same length as others, trimming the end'
                Logger.warning(str(tr) + msg)
                tr.data = tr.data[0:min_length]
    if linear:
        first_subspace = stacking.linstack(stream_list)
    else:
        first_subspace = stacking.PWS_stack(streams=stream_list)
    second_subspace = first_subspace.copy()
    for i in range(len(second_subspace)):
        second_subspace[i].data = np.diff(second_subspace[i].data)
        delta = second_subspace[i].stats.delta
        second_subspace[i].stats.starttime += 0.5 * delta

    return [first_subspace, second_subspace]


def svd_to_stream(uvectors, stachans, k, sampling_rate):
    """
    Convert the singular vectors output by SVD to streams.

    One stream will be generated for each singular vector level,
    for all channels.  Useful for plotting, and aiding seismologists thinking
    of waveforms!

    :type svectors: list
    :param svectors: List of :class:`numpy.ndarray` Singular vectors
    :type stachans: list
    :param stachans: List of station.channel Strings
    :type k: int
    :param k: Number of streams to return = number of SV's to include
    :type sampling_rate: float
    :param sampling_rate: Sampling rate in Hz

    :returns:
        svstreams, List of :class:`obspy.core.stream.Stream`, with
        svStreams[0] being composed of the highest rank singular vectors.
    """
    svstreams = []
    for i in range(k):
        svstream = []
        for j, stachan in enumerate(stachans):
            if len(uvectors[j]) <= k:
                Logger.warning('Too few traces at %s for a %02d dimensional '
                               'subspace. Detector streams will not include '
                               'this channel.' % ('.'.join(stachan[0],
                                                           stachan[1]), k))
            else:
                svstream.append(Trace(uvectors[j][i],
                                      header={'station': stachan[0],
                                              'channel': stachan[1],
                                              'sampling_rate': sampling_rate}))
        svstreams.append(Stream(svstream))
    return svstreams


def corr_cluster(trace_list, thresh=0.9):
    """
    Group traces based on correlations above threshold with the stack.

    Will run twice, once with 80% of threshold threshold to remove large
    outliers that would negatively affect the stack, then again with your
    threshold.

    :type trace_list: list
    :param trace_list:
        List of :class:`obspy.core.stream.Trace` to compute similarity between
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
    init_thresh = thresh * .8
    stack = stacking.linstack([Stream(tr) for tr in trace_list])[0]
    output = np.array([False] * len(trace_list))
    group1 = []
    array_xcorr = get_array_xcorr()
    for i, tr in enumerate(trace_list):
        cc = array_xcorr(np.array([tr.data]), stack.data, [0])[0][0][0]
        if cc > init_thresh:
            output[i] = True
            group1.append(tr)
    if len(group1) == 0:
        Logger.warning('Nothing made it past the first 80% threshold')
        return output
    stack = stacking.linstack([Stream(tr) for tr in group1])[0]
    group2 = []
    for i, tr in enumerate(trace_list):
        if array_xcorr(
                np.array([tr.data]), stack.data, [0])[0][0][0] > thresh:
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
    :param detections: List of :class:`eqcorrscan.core.match_filter.Detection`.
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
        Files written will be multiplexed miniseed files, the encoding will
        be chosen automatically and will likely be float.
    :type extract_Z: bool
    :param extract_Z:
        Set to True to also extract Z channels for detections delays will be
        the same as horizontal channels, only applies if only horizontal
        channels were used in the template.
    :type additional_stations: list
    :param additional_stations:
        List of tuples of (station, channel) to also extract data
        for using an average delay.

    :returns: list of :class:`obspy.core.streams.Stream`
    :rtype: list

    .. rubric: Example

    >>> from eqcorrscan.utils.clustering import extract_detections
    >>> from eqcorrscan.core.match_filter import Detection
    >>> from obspy import read, UTCDateTime
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> # Use some dummy detections, you would use real one
    >>> detections = [Detection(
    ...     template_name='temp1', detect_time=UTCDateTime(2012, 3, 26, 9, 15),
    ...     no_chans=2, chans=['WHYM', 'EORO'], detect_val=2, threshold=1.2,
    ...     typeofdet='corr', threshold_type='MAD', threshold_input=8.0),
    ...               Detection(
    ...     template_name='temp2', detect_time=UTCDateTime(2012, 3, 26, 18, 5),
    ...     no_chans=2, chans=['WHYM', 'EORO'], detect_val=2, threshold=1.2,
    ...     typeofdet='corr', threshold_type='MAD', threshold_input=8.0)]
    >>> archive = os.path.join(TEST_PATH, 'day_vols')
    >>> template_files = [os.path.join(TEST_PATH, 'temp1.ms'),
    ...                   os.path.join(TEST_PATH, 'temp2.ms')]
    >>> templates = [('temp' + str(i), read(filename))
    ...              for i, filename in enumerate(template_files)]
    >>> extracted = extract_detections(detections, templates,
    ...                                archive=archive, arc_type='day_vols')
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
    >>> # Extract from stations not included in the detections
    >>> extracted = extract_detections(
    ...    detections, templates, archive=archive, arc_type='day_vols',
    ...    additional_stations=[('GOVA', 'SHZ')])
    >>> print(extracted[0].sort())
    3 Trace(s) in Stream:
    AF.EORO..SHZ | 2012-03-26T09:14:15.000000Z - 2012-03-26T09:15:45.000000Z |\
 1.0 Hz, 91 samples
    AF.GOVA..SHZ | 2012-03-26T09:14:15.000000Z - 2012-03-26T09:15:45.000000Z |\
 1.0 Hz, 91 samples
    AF.WHYM..SHZ | 2012-03-26T09:14:15.000000Z - 2012-03-26T09:15:45.000000Z |\
 1.0 Hz, 91 samples
    >>> # The detections can be saved to a file:
    >>> extract_detections(detections, templates, archive=archive,
    ...                    arc_type='day_vols',
    ...                    additional_stations=[('GOVA', 'SHZ')], outdir='.')
    """
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
        Logger.info('Adding additional stations')
        for t, template in enumerate(all_stachans):
            av_delay = np.mean(all_delays[t][1])
            for sta in additional_stations:
                if sta not in template[1]:
                    Logger.info('Added station ' + '.'.join(sta))
                    template[1].append(sta)
                    all_delays[t][1].append(av_delay)
    del stachans
    # Loop through the days
    for detection_day in detection_days:
        Logger.info('Working on detections for day: ' + str(detection_day))
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
            Logger.info(
                'Cutting for detections at: ' +
                detection.detect_time.strftime('%Y/%m/%d %H:%M:%S'))
            detect_wav = st.copy()
            for tr in detect_wav:
                t1 = UTCDateTime(detection.detect_time) - extract_len / 2
                t2 = UTCDateTime(detection.detect_time) + extract_len / 2
                tr.trim(starttime=t1, endtime=t2)
            if outdir:
                if not os.path.isdir(os.path.join(outdir,
                                                  detection.template_name)):
                    os.makedirs(os.path.join(outdir, detection.template_name))
                detect_wav.write(os.path.join(outdir, detection.template_name,
                                              detection.detect_time.
                                              strftime('%Y-%m-%d_%H-%M-%S') +
                                              '.ms'),
                                 format='MSEED')
                Logger.info(
                    'Written file: %s' % '/'.join(
                        [outdir, detection.template_name,
                         detection.detect_time.strftime('%Y-%m-%d_%H-%M-%S')
                         + '.ms']))
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


def remove_unclustered(catalog, distance_cutoff, num_threads=None):
    """
    Remove events in catalog which do not have any other nearby events.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog for which to compute the distance matrix
    :type distance_cutoff: float
    :param distance_cutoff: Cutoff for considering events unclustered in km

    :returns: catalog
    :rtype: :class:`obspy.core.event.Catalog`
    """
    import ctypes
    from eqcorrscan.utils.libnames import _load_cdll
    from math import radians

    utilslib = _load_cdll('libutils')

    utilslib.remove_unclustered.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.uint8,
                               flags='C_CONTIGUOUS'),
        ctypes.c_float, ctypes.c_int]
    utilslib.remove_unclustered.restype = ctypes.c_int

    # Initialize square matrix
    mask = np.ascontiguousarray(np.zeros(len(catalog), dtype=np.uint8))
    latitudes, longitudes, depths = (
        np.empty(len(catalog), dtype=np.float32),
        np.empty(len(catalog), dtype=np.float32),
        np.empty(len(catalog), dtype=np.float32))
    for i, event in enumerate(catalog):
        origin = event.preferred_origin() or event.origins[0]
        latitudes[i] = radians(origin.latitude)
        longitudes[i] = radians(origin.longitude)
        depths[i] = origin.depth / 1000
    depths = np.ascontiguousarray(depths, dtype=np.float32)
    latitudes = np.ascontiguousarray(latitudes, dtype=np.float32)
    longitudes = np.ascontiguousarray(longitudes, dtype=np.float32)

    if num_threads is None:
        # Testing showed that 400 events per thread was best on the i7.
        num_threads = int(min(cpu_count(), len(catalog) // 400))
    if num_threads == 0:
        num_threads = 1

    ret = utilslib.remove_unclustered(
        latitudes, longitudes, depths, len(catalog), mask, distance_cutoff,
        num_threads)

    if ret != 0:  # pragma: no cover
        raise Exception("Internal error while computing distance matrix")

    _events = []
    for i, event in enumerate(catalog.events):
        if mask[i]:
            _events.append(event)
    catalog.events = _events
    return catalog


def dist_mat_km(catalog, num_threads=None):
    """
    Compute the distance matrix for a catalog using hypocentral separation.

    Will give physical distance in kilometers.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog for which to compute the distance matrix

    :returns: distance matrix
    :rtype: :class:`numpy.ndarray`
    """
    import ctypes
    from eqcorrscan.utils.libnames import _load_cdll

    utilslib = _load_cdll('libutils')

    utilslib.distance_matrix.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_long,
        np.ctypeslib.ndpointer(dtype=np.float32,
                               flags='C_CONTIGUOUS'),
        ctypes.c_int]
    utilslib.distance_matrix.restype = ctypes.c_int

    # Initialize square matrix
    dist_mat = np.zeros((len(catalog), len(catalog)), dtype=np.float32)
    latitudes, longitudes, depths = (
        np.empty(len(catalog)), np.empty(len(catalog)), np.empty(len(catalog)))
    for i, event in enumerate(catalog):
        origin = event.preferred_origin() or event.origins[0]
        latitudes[i] = origin.latitude
        longitudes[i] = origin.longitude
        depths[i] = origin.depth / 1000
    depths = np.ascontiguousarray(depths, dtype=np.float32)
    latitudes = np.ascontiguousarray(np.radians(latitudes), dtype=np.float32)
    longitudes = np.ascontiguousarray(np.radians(longitudes), dtype=np.float32)

    if num_threads is None:
        # Testing showed that 400 events per thread was best on the i7.
        num_threads = int(min(cpu_count(), len(catalog) // 400))
    if num_threads == 0:
        num_threads = 1

    ret = utilslib.distance_matrix(
        latitudes, longitudes, depths, len(catalog), dist_mat, num_threads)

    if ret != 0:  # pragma: no cover
        raise Exception("Internal error while computing distance matrix")
    # Fill distance matrix
    out = dist_mat.T + dist_mat
    return out


def dist_mat_time(catalog):
    """
    Compute the distance matrix for all a catalog using origin-time difference.

    Will give temporal separation in seconds.

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
            master_ori = master.origins[-1]
        for slave in catalog:
            if slave.preferred_origin():
                slave_ori = slave.preferred_origin()
            else:
                slave_ori = slave.origins[-1]
            mast_list.append(abs(master_ori.time - slave_ori.time))
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
    Cluster a catalog by distance only. DEPRECEATED - use catalog_cluster

    Will compute the matrix of physical distances between events and utilize
    the :mod:`scipy.clustering.hierarchy` module to perform the clustering.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog of events to clustered
    :type d_thresh: float
    :param d_thresh: Maximum inter-event distance threshold

    :returns: list of :class:`obspy.core.event.Catalog` objects
    :rtype: list
    """
    Logger.warning(
        "Depreciated, use `catalog_cluster` with `metric='distance'`")
    return catalog_cluster(catalog=catalog, thresh=d_thresh, metric="distance",
                           show=show)


def catalog_cluster(catalog, thresh, metric="distance", show=True):
    """
    Cluster a catalog by distance only.

    Will compute the matrix of physical distances between events and utilize
    the :mod:`scipy.clustering.hierarchy` module to perform the clustering.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog of events to clustered
    :type thresh: float
    :param thresh:
        Maximum separation, either in km (`metric="distance"`) or in seconds
        (`metric="time"`)
    :type metric: str
    :param metric: Either "distance" or "time"

    :returns: list of :class:`obspy.core.event.Catalog` objects
    :rtype: list

    >>> from eqcorrscan.utils.clustering import catalog_cluster
    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> client = Client("NCEDC")
    >>> starttime = UTCDateTime("2002-01-01")
    >>> endtime = UTCDateTime("2002-02-01")
    >>> cat = client.get_events(starttime=starttime, endtime=endtime,
    ...                         minmagnitude=2)
    >>> groups = catalog_cluster(catalog=cat, thresh=2, show=False)

    >>> from eqcorrscan.utils.clustering import catalog_cluster
    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> client = Client("https://earthquake.usgs.gov")
    >>> starttime = UTCDateTime("2002-01-01")
    >>> endtime = UTCDateTime("2002-02-01")
    >>> cat = client.get_events(starttime=starttime, endtime=endtime,
    ...                         minmagnitude=6)
    >>> groups = catalog_cluster(catalog=cat, thresh=1000, metric="time",
    ...     show=False)
    """
    # Compute the distance matrix and linkage
    if metric == "distance":
        dist_mat = dist_mat_km(catalog)
    elif metric == "time":
        dist_mat = dist_mat_time(catalog)
    else:
        raise NotImplementedError("Only supports distance and time metrics")
    dist_vec = squareform(dist_mat)
    Z = linkage(dist_vec, method='average')

    # Cluster the linkage using the given threshold as the cutoff
    indices = fcluster(Z, t=thresh, criterion='distance')
    group_ids = list(set(indices))
    indices = [(indices[i], i) for i in range(len(indices))]

    if show:
        # Plot the dendrogram...if it's not way too huge
        dendrogram(Z, color_threshold=thresh, distance_sort='ascending')
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
    >>> client = Client("https://earthquake.usgs.gov")
    >>> starttime = UTCDateTime("2002-01-01")
    >>> endtime = UTCDateTime("2002-02-01")
    >>> cat = client.get_events(starttime=starttime, endtime=endtime,
    ...                         minmagnitude=6)
    >>> groups = space_time_cluster(catalog=cat, t_thresh=86400, d_thresh=1000)
    """
    initial_spatial_groups = catalog_cluster(
        catalog=catalog, thresh=d_thresh, metric="distance", show=False)
    # Need initial_spatial_groups to be lists at the moment
    initial_spatial_lists = []
    for group in initial_spatial_groups:
        initial_spatial_lists.append(list(group))
    # Check within these groups and throw them out if they are not close in
    # time.
    groups = []
    for group in initial_spatial_lists:
        if len(group) > 1:
            sub_time_cluster = catalog_cluster(
                catalog=group, thresh=t_thresh, metric="time", show=False)
            groups.extend(sub_time_cluster)
        else:
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
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> det_file = os.path.join(TEST_PATH, 'expected_tutorial_detections.txt')
    >>> detections = re_thresh_csv(path=det_file, old_thresh=8, new_thresh=10,
    ...                            chan_thresh=3)
    >>> print(len(detections))
    17

    .. Note::
        This is a legacy function, and will read detections from all versions.

    .. Warning:: Only works if thresholding was done by MAD.
    """
    from eqcorrscan.core.match_filter import read_detections
    Logger.warning('Legacy function, please use '
                   'eqcorrscan.core.match_filter.Party.rethreshold.')
    old_detections = read_detections(path)
    old_thresh = float(old_thresh)
    new_thresh = float(new_thresh)
    # Be nice, ensure that the thresholds are float
    detections = []
    detections_in = 0
    detections_out = 0
    for detection in old_detections:
        detections_in += 1
        con1 = (new_thresh / old_thresh) * detection.threshold
        con2 = detection.no_chans >= chan_thresh
        required_thresh = (new_thresh / old_thresh) * detection.threshold
        con3 = abs(detection.detect_val) >= required_thresh
        if all([con1, con2, con3]):
            detections_out += 1
            detections.append(detection)
    Logger.info('Read in %i detections' % detections_in)
    Logger.info('Left with %i detections' % detections_out)
    return detections


if __name__ == "__main__":
    import doctest

    doctest.testmod()
