#!/usr/bin/python
"""
Functions to cluster seismograms by a range of constraints.

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
import warnings


def cross_chan_coherence(st1, st2, i=0):
    """
    Function to determine the cross-channel coherancy between two streams of \
    multichannel seismic data.

    :type st1: obspy Stream
    :param st1: Stream one
    :type st2: obspy Stream
    :param st2: Stream two
    :type i: int
    :param i: index used for parallel async processing, returned unaltered

    :returns: cross channel coherence, float - normalized by number of\
        channels, if i, returns tuple of (cccoh, i) where i is int, as intput.
    """

    from eqcorrscan.core.match_filter import normxcorr2
    cccoh = 0.0
    kchan = 0
    for tr in st1:
        tr1 = tr.data
        # Assume you only have one waveform for each channel
        tr2 = st2.select(station=tr.stats.station,
                         channel=tr.stats.channel)
        if tr2:
            cccoh += normxcorr2(tr1, tr2[0].data)[0][0]
            kchan += 1
    if kchan:
        cccoh = cccoh / kchan
        return (cccoh, i)
    else:
        warnings.warn('No matching channels')
        return (0, i)


def distance_matrix(stream_list, cores=1):
    """
    Function to compute the distance matrix for all templates - will give \
    distance as 1-abs(cccoh), e.g. a well correlated pair of templates will \
    have small distances, and an equally well correlated reverse image will \
    have the same distance as apositively correlated image - this is an issue.

    :type stream_list: List of obspy.Streams
    :param stream_list: List of the streams to compute the distance matrix for
    :type cores: int
    :param cores: Number of cores to parallel process using, defaults to 1.

    :returns: ndarray - distance matrix
    """
    from multiprocessing import Pool

    # Initialize square matrix
    dist_mat = np.array([np.array([0.0] * len(stream_list))] *
                        len(stream_list))
    for i, master in enumerate(stream_list):
        # Start a parallel processing pool
        pool = Pool(processes=cores)
        # Parallel processing
        results = [pool.apply_async(cross_chan_coherence, args=(master,
                                                                stream_list[j],
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


def cluster(stream_list, show=True, corr_thresh=0.3, save_corrmat=False,
            cores='all', debug=1):
    """
    Function to take a set of templates and cluster them, will return groups \
    as lists of streams.  Clustering is done by computing the cross-channel \
    correlation sum of each stream in stream_list with every other stream in \
    the list.  Scipy.cluster.hierachy functions are then used to compute the \
    complete distance matrix, where distance is 1 minus the normalised \
    cross-correlation sum such that larger distances are less similar events. \
    Groups are then created by clustering the distance matrix at distances \
    less than 1 - corr_thresh.

    Will compute the distance matrix in parallel, using all available cores

    :type stream_list: List of Obspy.Stream
    :param stream_list: List of templates to compute clustering for
    :type show: bool
    :param show: plot linkage on screen if True, defaults to True
    :type corr_thresh: float
    :param corr_thresh: Cross-channel correlation threshold for grouping
    :type save_corrmat: bool
    :param save_corrmat: If True will save the distance matrix to dist_mat.npy
    :type cores: int
    :param cores: numebr of cores to use when computing the distance matrix, \
        defaults to 'all' which will work out how many cpus are available \
        and hog them.
    :type debug: int
    :param debug: Level of debugging from 1-5, higher is more output, \
        currently only level 1 implimented.

    :returns: List of groups with each group a list of streams making up \
        that group.
    """
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    import matplotlib.pyplot as plt
    from multiprocessing import cpu_count
    if cores == 'all':
        num_cores = cpu_count()
    else:
        num_cores = cores
    # Compute the distance matrix
    if debug >= 1:
        msg = ' '.join(['Computing the distance matrix using',
                        str(num_cores), 'cores'])
        print msg
    dist_mat = distance_matrix(stream_list, cores=num_cores)
    if save_corrmat:
        np.save('dist_mat.npy', dist_mat)
        if debug >= 1:
            print 'Saved the distance matrix as dist_mat.npy'
    dist_vec = squareform(dist_mat)
    # plt.matshow(dist_mat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnB)
    if debug >= 1:
        print 'Computing linkage'
    Z = linkage(dist_vec)
    if show:
        if debug >= 1:
            print 'Plotting the dendrogram'
        dendrogram(Z, color_threshold=1 - corr_thresh,
                   distance_sort='ascending')
        plt.show()
    # Get the indeces of the groups
    if debug >= 1:
        print 'Clustering'
    indeces = fcluster(Z, 1 - corr_thresh, 'distance')
    group_ids = list(set(indeces))  # Unique list of group ids
    if debug >= 1:
        msg = ' '.join(['Found', str(len(group_ids)), 'groups'])
        print msg
    # Convert to tuple of (group id, stream id)
    indeces = [(indeces[i], i) for i in xrange(len(indeces))]
    # Sort by group id
    indeces.sort(key=lambda tup: tup[0])
    groups = []
    if debug >= 1:
        print 'Extracting and grouping'
    for group_id in group_ids:
        group = []
        for ind in indeces:
            if ind[0] == group_id:
                group.append(stream_list[ind[1]])
            elif ind[0] > group_id:
                # Because we have sorted by group id, when the index is greater
                # than the group_id we can break the inner loop.
                # Patch applied by CJC 05/11/2015
                groups.append(group)
                break
    return groups


def group_delays(stream_list):
    """
    Function to group template waveforms according to their delays.

    :type stream_list: list of obspy.Stream
    :param stream_list: List of the waveforms you want to group

    :returns: list of List of obspy.Streams where each initial list is a \
        group with the same delays.
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
        print msg
        # Calculate the delays
        starttimes = []
        chans = []
        for tr in st:
            starttimes.append(tr.stats.starttime)
            chans.append((tr.stats.station, tr.stats.channel))
        # This delays calculation will be an issue if we
        # have changes in channels
        delays = [starttimes[m] - min(starttimes)
                  for m in xrange(len(starttimes))]
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


def SVD(stream_list):
    """
    Function to compute the SVD of a number of templates and return the \
    singular vectors and singular values of the templates.

    :type stream_list: List of Obspy.Stream
    :param stream_list: List of the templates to be analysed

    :return: SVector(list of ndarray), SValues(list) for each channel, \
        Uvalues(list of ndarray) for each channel, \
        stachans, List of String (station.channel)

    .. note:: We recommend that you align the data before computing the \
        SVD, e.g., the P-arrival on all templates for the same channel \
        should appear at the same time in the trace.
    """
    # Convert templates into ndarrays for each channel
    # First find all unique channels:
    stachans = []
    for st in stream_list:
        for tr in st:
            stachans.append(tr.stats.station+'.'+tr.stats.channel)
    stachans = list(set(stachans))
    # Initialize a list for the output matrices, one matrix per-channel
    SValues = []
    SVectors = []
    Uvectors = []
    for stachan in stachans:
        chan_mat = [stream_list[i].select(station=stachan.split('.')[0],
                                          channel=stachan.
                                          split('.')[1])[0].data
                    for i in range(len(stream_list)) if
                    len(stream_list[i].select(station=stachan.split('.')[0],
                                              channel=stachan.split('.')[1]))
                    != 0]
        chan_mat = np.asarray(chan_mat)
        print chan_mat.shape
        U, s, V = np.linalg.svd(chan_mat, full_matrices=True)
        SValues.append(s)
        SVectors.append(V)
        Uvectors.append(U)
    return SVectors, SValues, Uvectors, stachans


def empirical_SVD(stream_list, linear=True):
    """
    Empirical subspace detector generation function.  Takes a list of \
    templates and computes the stack as the first order subspace detector, \
    and the differential of this as the second order subspace detector \
    following the emprical subspace method of Barrett & Beroza, 2014 - SRL.

    :type stream_list: list of stream
    :param stream_list: list of template streams to compute the subspace \
        detectors from
    :type linear: bool
    :param linear: Set to true by default to compute the linear stack as the \
        first subspace vector, False will use the phase-weighted stack as the \
        first subspace vector.

    :returns: list of two streams
    """
    from eqcorrscan.utils import stacking
    if linear:
        first_subspace = stacking.linstack(stream_list)
    second_subspace = first_subspace.copy()
    for i in range(len(second_subspace)):
        second_subspace[i].data = np.diff(second_subspace[i].data)
        second_subspace[i].stats.starttime += 0.5 * \
            second_subspace[i].stats.delta
    return [first_subspace, second_subspace]


def SVD_2_stream(SVectors, stachans, k, sampling_rate):
    """
    Function to convert the singular vectors output by SVD to streams, one \
    for each singular vector level, for all channels.

    :type SVectors: List of np.ndarray
    :param SVectors: Singular vectors
    :type stachans: List of Strings
    :param stachans: List of station.channel Strings
    :type k: int
    :param k: Number of streams to return = number of SV's to include
    :type sampling_rate: float
    :param sampling_rate: Sampling rate in Hz

    :returns: SVstreams, List of Obspy.Stream, with SVStreams[0] being \
        composed of the highest rank singular vectors.
    """
    from obspy import Stream, Trace
    SVstreams = []
    for i in range(k):
        SVstream = []
        for j, stachan in enumerate(stachans):
            SVstream.append(Trace(SVectors[j][i],
                                  header={'station': stachan.split('.')[0],
                                          'channel': stachan.split('.')[1],
                                          'sampling_rate': sampling_rate}))

        SVstreams.append(Stream(SVstream))
    return SVstreams


def corr_cluster(trace_list, thresh=0.9):
    """
    Group traces based on correlations above threshold with the stack - will \
    run twice, once with a lower threshold, then again with your threshold to \
    remove large outliers.

    :type trace_list: list of obspy.Trace
    :param trace_list: Traces to compute similarity between
    :type thresh: float
    :param thresh: Correlation threshold between -1-1

    :returns: np.ndarray of bool
    """
    from eqcorrscan.utils import stacking
    from obspy import Stream
    from core.match_filter import normxcorr2
    stack = stacking.linstack([Stream(tr) for tr in trace_list])[0]
    output = np.array([False]*len(trace_list))
    group1 = []
    for i, tr in enumerate(trace_list):
        if normxcorr2(tr.data, stack.data)[0][0] > 0.6:
            output[i] = True
            group1.append(tr)
    stack = stacking.linstack([Stream(tr) for tr in group1])[0]
    group2 = []
    for i, tr in enumerate(trace_list):
        if normxcorr2(tr.data, stack.data)[0][0] > thresh:
            group2.append(tr)
            output[i] = True
        else:
            output[i] = False
    return output


def extract_detections(detections, templates, contbase_list, extract_len=90.0,
                       outdir=None, extract_Z=True, additional_stations=[]):
    """
    Function to extract the waveforms associated with each detection in a \
    list of detections for the template, template.  Waveforms will be \
    returned as a list of obspy.Streams containing segments of extract_len.  \
    They will also be saved if outdir is set.  The default is unset.  The \
    default extract_len is 90 seconds per channel.

    :type detections: list tuple of of :class: datetime.datetime, string
    :param detections: List of datetime objects, and their associated \
        template name.
    :type templates: list of tuple of string and :class: obspy.Stream
    :param templates: A list of the tuples of the template name and the \
        template Stream used to detect detections.
    :type contbase_list: list of tuple of string
    :param contbase_list: List of tuples of the form \
        ['path', 'type', 'network']  Where path is the path to the continuous \
        database, type is the directory structure, which can be either \
        Yyyyy/Rjjj.01, which is the standard IRIS Year, julian day structure, \
        or, yyyymmdd which is a single directory for every day.
    :type extract_len: float
    :param extract_len: Length to extract around the detection (will be \
        equally cut around the detection time) in seconds.  Default is 90.0.
    :type outdir: bool or str
    :param outdir: Default is None, with None set, no files will be saved, \
        if set each detection will be saved into this directory with files \
        named according to the detection time, NOT than the waveform \
        start time. Detections will be saved into template subdirectories.
    :type extract_Z: bool
    :param extract_Z: Set to True to also extract Z channels for detections \
        delays will be the same as horizontal channels, only applies if \
        only horizontal channels were used in the template.
    :type additional_stations: list of tuple
    :param additional_stations: List of stations, chanels and networks to \
        also extract data for using an average delay.

    :returns: list of :class: obspy.Stream
    """
    from obspy import read
    from obspy import UTCDateTime
    import os
    # Sort the template according to starttimes, needed so that stachan[i]
    # corresponds to delays[i]
    all_delays = []  # List of tuples of template name, delays
    all_stachans = []
    for template in templates:
        templatestream = template[1].sort(['starttime'])
        stachans = [(tr.stats.station, tr.stats.channel, tr.stats.network)
                    for tr in templatestream]
        mintime = templatestream[0].stats.starttime
        delays = [tr.stats.starttime - mintime for tr in templatestream]
        all_delays.append((template[0], delays))
        all_stachans.append((template[0], stachans))
    # Sort the detections and group by day
    detections.sort()
    detection_days = [detection[0].date() for detection in detections]
    detection_days = list(set(detection_days))
    detection_days.sort()

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
                    new_stachans.append((stachan[0], stachan[1][0]+'Z',
                                         stachan[2]))
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
        for stachan in stachans:
            print('Extracting data for ' + '.'.join(stachan))
            contbase = [base for base in contbase_list
                        if base[2] == stachan[2]][0]
            if contbase[1] == 'yyyymmdd':
                dayfile = detection_day.strftime('%Y%m%d') + '/*' +\
                    stachan[0] + '.' + stachan[1][0] + '?' + stachan[1][-1] +\
                    '.*'
            elif contbase[1] == 'Yyyyy/Rjjj.01':
                dayfile = detection_day.strftime('Y%Y/R%j.01')+'/'+stachan[0] +\
                    '.*.'+stachan[1][0]+'?'+stachan[1][-1]+'.' +\
                    detection_day.strftime('%Y.%j')
            if 'st' not in locals():
                try:
                    st = read(contbase[0]+'/'+dayfile)
                except:
                    print 'No data for '+contbase[0]+'/'+dayfile
            else:
                try:
                    st += read(contbase[0]+'/'+dayfile)
                except:
                    print 'No data for '+contbase[0]+'/'+dayfile
        st.merge(fill_value='interpolate')
        day_detections = [detection for detection in detections
                          if detection[0].date() == detection_day]
        del stachans, delays
        for detection in day_detections:
            template = detection[1]
            t_stachans = [stachans[1] for stachans in all_stachans
                          if stachans[0] == template][0]
            t_delays = [delays[1] for delays in all_delays
                        if delays[0] == template][0]
            print 'Cutting for detections at: ' +\
                detection[0].strftime('%Y/%m/%d %H:%M:%S')
            detect_wav = st.copy()
            for tr in detect_wav:
                tr.trim(starttime=UTCDateTime(detection[0]) - extract_len / 2,
                        endtime=UTCDateTime(detection[0]) + extract_len / 2)
            if outdir:
                if not os.path.isdir(outdir+'/'+template):
                    os.makedirs(outdir+'/'+template)
                detect_wav.write(outdir+'/'+template+'/' +
                                 detection[0].strftime('%Y-%m-%d_%H-%M-%S') +
                                 '.ms', format='MSEED', encoding='STEIM2')
                print 'Written file: '+outdir+'/'+template+'/' +\
                    detection[0].strftime('%Y-%m-%d_%H-%M-%S')+'.ms'
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


def space_time_cluster(detections, t_thresh, d_thresh):
    """
    Function to cluster detections in space and time, use to seperate \
    repeaters from other events.

    :type detections: list
    :param detections: List of tuple of tuple of location (lat, lon, depth \
        (km)), and time as a datetime object
    :type t_thresh: float
    :param t_thresh: Maximum inter-event time threshold in seconds
    :type d_thresh: float
    :param d_thresh: Maximum inter-event distance in km

    :returns: List of tuple (detections, clustered) and list of indeces of \
        clustered detections
    """
    from eqcorrscan.utils.mag_calc import dist_calc
    # Ensure they are sorted by time first, not that we need it.
    detections.sort(key=lambda tup: tup[1])
    clustered = []
    clustered_indeces = []
    for master_ind, master in enumerate(detections):
        keep = False
        for slave in detections:
            if not master == slave and\
               abs((master[1] - slave[1]).total_seconds()) <= t_thresh and \
               dist_calc(master[0], slave[0]) <= d_thresh:
                # If the slave events is close in time and space to the master
                # keep it and break out of the loop.
                keep = True
                break
        if keep:
            clustered.append(master)
            clustered_indeces.append(master_ind)

    return clustered, clustered_indeces


def re_thresh_csv(path, old_thresh, new_thresh, chan_thresh):
    """
    Function to remove detections by changing the threshold, can only be done \
    to remove detection by increasing threshold, threshold lowering will have \
    no effect.

    :type path: str
    :param path: Path to the .csv detection file
    :type old_thresh: float
    :param old_thresh: Old threshold MAD multiplier
    :type new_thresh: float
    :param new_thresh: New threhsold MAD multiplier
    :type chan_thresh: int
    :param chan_thresh: Minimum number of channels for a detection

    returns: List of detections
    """
    f = open(path, 'r')
    old_thresh = float(old_thresh)
    new_thresh = float(new_thresh)
    # Be nice, ensure that the thresholds are float
    detections = []
    detections_in = 0
    detections_out = 0
    for line in f:
        if not line.split(', ')[0] == 'template' and len(line) > 2:
            detections_in += 1
            if abs(float(line.split(', ')[3])) >=\
               (new_thresh/old_thresh)*float(line.split(', ')[2]) and\
               int(line.split(', ')[4]) >= chan_thresh:
                detections_out += 1
                detections.append(line.split(', '))
    print 'Read in '+str(detections_in)+' detections'
    print 'Left with '+str(detections_out)+' detections'
    return detections


if __name__ == "__main__":
    import doctest
    doctest.testmod()
