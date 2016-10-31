#!/usr/bin/python
r"""
Code to determine the brightness function of seismic data according to a
three-dimensional travel-time grid.  This travel-time grid can be
generated using the grid2time function of the NonLinLoc package by Anthony
Lomax which can be found here: `NonLinLoc <http://alomax.free.fr/nlloc/>`_.
NonLinLoc is not distributed within this package but it is a very useful
stand-alone library for seismic event location.

This code is based on the method of
`Frank & Shapiro 2014
<http://gji.oxfordjournals.org/content/197/2/1215.short>`_.

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
import csv
import glob
import sys
import os

from obspy import Stream, Trace, read as obsread
from multiprocessing import Pool, cpu_count
from copy import deepcopy
from obspy.core.event import Catalog, Event, Pick, WaveformStreamID, Origin
from obspy.core.event import EventDescription, CreationInfo, Comment

from eqcorrscan.core.match_filter import DETECTION, normxcorr2
from eqcorrscan.utils import findpeaks
from eqcorrscan.core.template_gen import _template_gen


def _read_tt(path, stations, phase, phaseout='S', ps_ratio=1.68,
             lags_switch=True):
    """
    Read in .csv files of slowness generated from Grid2Time.

    Converts these data to a useful format here.

    It should be noted that this can read either P or S travel-time grids, not
    both at the moment.

    :type path: str
    :param path: The path to the .csv Grid2Time outputs
    :type stations: list
    :param stations: List of station names to read slowness files for.
    :type phase: str
    :param phase: Input phase type.
    :type phaseout: str
    :param phaseout: What phase to return the lagtimes in.
    :type ps_ratio: float
    :param ps_ratio: p to s ratio for conversion
    :type lags_switch: bool
    :param lags_switch: Return lags or raw travel-times, if set to true will \
        return lags.

    :returns: Stations
    :rtype: list
    :returns: List of lists of tuples of node locations
    :rtype: list
    :returns: Array of lags.
    :rtype: :class:`numpy.ndarray`

    .. note::
        **Output:**

        station[1] refers to nodes[1] and lags[1] nodes[1][1] refers
        to station[1] and lags[1][1] nodes[n][n] is a tuple of latitude,
        longitude and depth.

    .. note:: This function currently needs comma separated grid files in \
        NonLinLoc format.  Only certain versions of NonLinLoc write these csv \
        files, however it should be possible to read the binary files \
        directly.  If you find you need this capability let us know and we \
        can try and implement it.
    """
    # Locate the slowness file information
    gridfiles = []
    stations_out = []
    for station in stations:
        gridfiles += (glob.glob(path + '*.' + phase + '.' + station +
                      '.time.csv'))
        if glob.glob(path + '*.' + phase + '.' + station + '*.csv'):
            stations_out += [station]
    if not stations_out:
        raise IOError('No slowness files found')
    # Read the files
    allnodes = []
    for gridfile in gridfiles:
        print('     Reading slowness from: ' + gridfile)
        f = open(gridfile, 'r')
        grid = csv.reader(f, delimiter=str(' '))
        traveltime = []
        nodes = []
        for row in grid:
            nodes.append((float(row[0]), float(row[1]), float(row[2])))
            traveltime.append(float(row[3]))
        traveltime = np.array(traveltime)
        if not phase == phaseout:
            if phase == 'S':
                traveltime = traveltime / ps_ratio
            else:
                traveltime = traveltime * ps_ratio
        if lags_switch:
            lags = traveltime - min(traveltime)
        else:
            lags = traveltime
        if 'alllags' not in locals():
            alllags = [lags]
        else:
            alllags = np.concatenate((alllags, [lags]), axis=0)
        allnodes = nodes
        # each element of allnodes should be the same as the
        # other one, e.g. for each station the grid must be the
        # same, hence allnodes=nodes
        f.close()
    alllags = np.array(alllags)
    return stations_out, allnodes, alllags


def _resample_grid(stations, nodes, lags, mindepth, maxdepth, corners):
    """
    Resample the lagtime grid to a given volume.

    For use if the grid from Grid2Time is too large or you want to run a
    faster, downsampled scan.

    :type stations: list
    :param stations:
        List of station names from in the form where stations[i] refers to
        nodes[i][:] and lags[i][:]
    :type nodes: list
    :param nodes:
        List of node points where nodes[i] referes to stations[i] and
        nodes[:][:][0] is latitude in degrees, nodes[:][:][1] is longitude in
        degrees, nodes[:][:][2] is depth in km.
    :type lags: numpy.ndarray
    :param lags:
        Array of arrays where lags[i][:] refers to stations[i]. lags[i][j]
        should be the delay to the nodes[i][j] for stations[i] in seconds.
    :type mindepth: float
    :param mindepth: Upper limit of volume
    :type maxdepth: float
    :param maxdepth: Lower limit of volume
    :type corners: matplotlib.path.Path
    :param corners:
        matplotlib Path of the corners for the 2D polygon to cut to in lat and
        lon.

    :returns: Stations
    :rtype: list
    :returns: List of lists of tuples of node locations
    :rtype: list
    :returns: Array of lags.
    :rtype: :class:`numpy.ndarray`

    .. note::
        **Output:**

        station[1] refers to nodes[1] and lags[1] nodes[1][1] refers
        to station[1] and lags[1][1] nodes[n][n] is a tuple of latitude,
        longitude and depth.
    """
    resamp_nodes = []
    resamp_lags = []
    # Cut the volume
    for i, node in enumerate(nodes):
        # If the node is within the volume range, keep it
        if mindepth < float(node[2]) < maxdepth and\
           corners.contains_point(node[0:2]):
                resamp_nodes.append(node)
                resamp_lags.append([lags[:, i]])
    # Reshape the lags
    print(np.shape(resamp_lags))
    resamp_lags = np.reshape(resamp_lags, (len(resamp_lags), len(stations))).T
    # Resample the nodes - they are sorted in order of size with largest long
    # then largest lat, then depth.
    print(' '.join(['Grid now has ', str(len(resamp_nodes)), 'nodes']))
    return stations, resamp_nodes, resamp_lags


def _rm_similarlags(stations, nodes, lags, threshold):
    """
    Remove nodes that have a very similar network moveout to another node.

    This function will, for each node, calculate the difference in lagtime
    at each station at every node, then sum these for each node to get a
    cumulative difference in network moveout.  This will result in an
    array of arrays with zeros on the diagonal.

    :type stations: list
    :param stations:
        List of station names from in the form where stations[i] refers to
        nodes[i][:] and lags[i][:]
    :type nodes: list
    :param nodes:
        List of node points where nodes[i] referes to stations[i] and
        nodes[:][:][0] is latitude in degrees, nodes[:][:][1] is longitude in
        degrees, nodes[:][:][2] is depth in km.
    :type lags: numpy.ndarray
    :param lags:
        Array of arrays where lags[i][:] refers to stations[i]. lags[i][j]
        should be the delay to the nodes[i][j] for stations[i] in seconds.
    :type threshold: float
    :param threshold: Threshold for removal in seconds

    :returns: Stations
    :rtype: list
    :returns: List of lists of tuples of node locations
    :rtype: list
    :returns: Array of lags.
    :rtype: :class:`numpy.ndarray`

    .. note::
        **Output:**

        station[1] refers to nodes[1] and lags[1] nodes[1][1] refers
        to station[1] and lags[1][1] nodes[n][n] is a tuple of latitude,
        longitude and depth.
    """
    netdif = abs((lags.T -
                  lags.T[0]).sum(axis=1).reshape(1, len(nodes))) > threshold
    for i in range(len(nodes)):
        _netdif = abs((lags.T -
                       lags.T[i]).sum(axis=1).reshape(1, len(nodes)))\
            > threshold
        netdif = np.concatenate((netdif, _netdif), axis=0)
        sys.stdout.write("\r" + str(float(i) // len(nodes) * 100) + "% \r")
        sys.stdout.flush()
    nodes_out = [nodes[0]]
    node_indeces = [0]
    print("\n")
    print(len(nodes))
    for i in range(1, len(nodes)):
        if np.all(netdif[i][node_indeces]):
            node_indeces.append(i)
            nodes_out.append(nodes[i])
    lags_out = lags.T[node_indeces].T
    print("Removed " + str(len(nodes) - len(nodes_out)) + " duplicate nodes")
    return stations, nodes_out, lags_out


def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))


def _node_loop(stations, lags, stream, clip_level,
               i=0, mem_issue=False, instance=0, plot=False):
    """
    Internal function to allow for brightness to be paralleled.

    :type stations: list
    :param stations: List of stations to use.
    :type lags: numpy.ndarray
    :param lags: List of lags where lags[i[:]] are the lags for stations[i].
    :type stream: obspy.core.stream.Stream
    :param stream: Data stream to find the brightness for.
    :type clip_level: float
    :param clip_level: Upper limit for energy as a multiplier to the mean \
        energy.
    :type i: int
    :param i: Index of loop for parallelisation.
    :type mem_issue: bool
    :param mem_issue: If True will write to disk rather than storing data in \
        RAM.
    :type instance: int
    :param instance: instance for bulk parallelisation, only used if \
        mem_issue=true.
    :type plot: bool
    :param plot: Turn plotting on or off, defaults to False.

    :returns: index
    :rtype: int
    :returns: network response
    :rtype: numpy.ndarray
    """
    if plot:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(stream) + 1, 1, sharex=True)
        axes = axes.ravel()

    for l, tr in enumerate(stream):
        j = [k for k in range(len(stations))
             if stations[k] == tr.stats.station]
        # Check that there is only one matching station
        if len(j) > 1:
            warnings.warn('Too many stations')
            j = [j[0]]
        if len(j) == 0:
            warnings.warn('No station match')
            continue
        lag = lags[j[0]]
        pad = np.zeros(int(round(lag * tr.stats.sampling_rate)))
        lagged_energy = np.square(np.concatenate((tr.data, pad)))[len(pad):]
        # Clip energy
        lagged_energy = np.clip(lagged_energy, 0,
                                clip_level * np.mean(lagged_energy))
        if 'energy' not in locals():
            energy = (lagged_energy /
                      _rms(lagged_energy)).reshape(1, len(lagged_energy))
            # Cope with zeros enountered
            energy = np.nan_to_num(energy)
            # This is now an array of floats - we can convert this to int16
            # normalize to have max at max of int16 range
            if not max(energy[0]) == 0.0:
                scalor = 1 / max(energy[0])
                energy = (500 * (energy * scalor)).astype(np.int16)
            else:
                energy = energy.astype(np.int16)
        else:
            norm_energy = (lagged_energy /
                           _rms(lagged_energy)).reshape(1, len(lagged_energy))
            norm_energy = np.nan_to_num(norm_energy)
            # Convert to int16
            if not max(norm_energy[0]) == 0.0:
                scalor = 1 / max(norm_energy[0])
                norm_energy = (500 * (norm_energy * scalor)).astype(np.int16)
            else:
                norm_energy = norm_energy.astype(np.int16)
            # Apply lag to data and add it to energy - normalize the data here
            energy = np.concatenate((energy, norm_energy), axis=0)
        if plot:
            axes[l].plot(lagged_energy * 200, 'r')
            axes[l].plot(tr.data, 'k')
            # energy_tr=Trace(energy[l])
            energy_tr = Trace(lagged_energy)
            # energy_tr=Trace(tr.data)
            energy_tr.stats.station = tr.stats.station
            energy_tr.stats.sampling_rate = tr.stats.sampling_rate
            if 'energy_stream' not in locals():
                energy_stream = Stream(energy_tr)
            else:
                energy_stream += energy_tr
    energy = np.sum(energy, axis=0).reshape(1, len(lagged_energy))
    energy = energy.astype(np.uint16)
    # Convert any nans to zeros
    energy = np.nan_to_num(energy)
    if plot:
        energy_tr = Trace(energy[0])
        energy_tr.stats.station = 'Stack'
        energy_tr.stats.sampling_rate = tr.stats.sampling_rate
        energy_tr.stats.network = 'Energy'
        energy_stream += energy_tr
        axes[-1].plot(energy[0])
        plt.subplots_adjust(hspace=0)
        plt.show()
        # energy_stream.plot(equal_scale=False, size=(800,600))
    if not mem_issue:
        return i, energy
    else:
        if not os.path.isdir('tmp' + str(instance)):
            os.makedirs('tmp' + str(instance))
        np.save('tmp' + str(instance) + '/node_' + str(i), energy)
        return i, str('tmp' + str(instance) + '/node_' + str(i))


def _cum_net_resp(node_lis, instance=0):
    """
    Compute the cumulative network response by reading saved energy .npy files.

    :type node_lis: numpy.ndarray
    :param node_lis: List of nodes (ints) to read from
    :type instance: int
    :param instance: Instance flag for parallel workflows, defaults to 0.

    :returns: cumulative network response
    :rtype: numpy.ndarray
    :returns: node indeces for each sample of the cumulative network response.
    :rtype: list
    """
    cum_net_resp = np.load('tmp' + str(instance) +
                           '/node_' + str(node_lis[0]) + '.npy')[0]
    os.remove('tmp' + str(instance) + '/node_' + str(node_lis[0]) + '.npy')
    indeces = np.ones(len(cum_net_resp)) * node_lis[0]
    for i in node_lis[1:]:
        node_energy = np.load('tmp' + str(instance) + '/node_' +
                              str(i) + '.npy')[0]
        updated_indeces = np.argmax([cum_net_resp, node_energy], axis=0)
        temp = np.array([cum_net_resp, node_energy])
        cum_net_resp = np.array([temp[updated_indeces[j]][j]
                                 for j in range(len(updated_indeces))])
        del temp, node_energy
        updated_indeces[updated_indeces == 1] = i
        indeces = updated_indeces
        os.remove('tmp' + str(instance) + '/node_' + str(i) + '.npy')
    return cum_net_resp, indeces


def _find_detections(cum_net_resp, nodes, threshold, thresh_type,
                     samp_rate, realstations, length):
    """
    Find detections within the cumulative network response.

    :type cum_net_resp: numpy.ndarray
    :param cum_net_resp: Array of cumulative network response for nodes
    :type nodes: list
    :param nodes: Nodes associated with the source of energy in the \
        cum_net_resp
    :type threshold: float
    :param threshold: Threshold value
    :type thresh_type: str
    :param thresh_type: Either MAD (Median Absolute Deviation) or abs \
        (absolute) or RMS (Root Mean Squared)
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz
    :type realstations: list
    :param realstations:
        List of stations used to make the cumulative network response, will be
        reported in the :class:`eqcorrscan.core.match_filter.DETECTION`
    :type length: float
    :param length: Maximum length of peak to look for in seconds

    :returns:
        Detections as :class:`eqcorrscan.core.match_filter.DETECTION` objects.
    :rtype: list
    """
    cum_net_resp = np.nan_to_num(cum_net_resp)  # Force no NaNs
    if np.isnan(cum_net_resp).any():
        raise ValueError("Nans present")
    print('Mean of data is: ' + str(np.median(cum_net_resp)))
    print('RMS of data is: ' + str(np.sqrt(np.mean(np.square(cum_net_resp)))))
    print('MAD of data is: ' + str(np.median(np.abs(cum_net_resp))))
    if thresh_type == 'MAD':
        thresh = (np.median(np.abs(cum_net_resp)) * threshold)
    elif thresh_type == 'abs':
        thresh = threshold
    elif thresh_type == 'RMS':
        thresh = _rms(cum_net_resp) * threshold
    print('Threshold is set to: ' + str(thresh))
    print('Max of data is: ' + str(max(cum_net_resp)))
    peaks = findpeaks.find_peaks2_short(cum_net_resp, thresh,
                                        length * samp_rate, debug=0)
    detections = []
    if peaks:
        for peak in peaks:
            node = nodes[peak[1]]
            detections.append(DETECTION(template_name=str(node[0]) +
                                        '_' + str(node[1]) + '_' +
                                        str(node[2]),
                                        detect_time=peak[1] / samp_rate,
                                        no_chans=len(realstations),
                                        detect_val=peak[0], threshold=thresh,
                                        typeofdet='brightness',
                                        chans=realstations,
                                        id=str(node[0]) + '_' + str(node[1]) +
                                        '_' + str(node[2]) +
                                        str(peak[1] / samp_rate)))
    else:
        detections = []
    print('I have found ' + str(len(peaks)) + ' possible detections')
    return detections


def coherence(stream_in, stations=['all'], clip=False):
    """
    Determine the average network coherence of a given template or detection.

    You will want your stream to contain only signal as noise will reduce the
    coherence (assuming it is incoherent random noise).

    :type stream_in: obspy.core.stream.Stream
    :param stream_in: The stream of seismic data you want to calculate the \
            coherence for.
    :type stations: list
    :param stations: List of stations to use for coherence, default is all
    :type clip: tuple
    :param clip: Default is to use all the data given (`False`) - \
            tuple of start and end in seconds from start of trace

    :return: tuple of coherence and number of channels used.
    :rtype: tuple
    """
    stream = stream_in.copy()  # Copy the data before we remove stations
    # First check that all channels in stream have data of the same length
    maxlen = np.max([len(tr.data) for tr in stream])
    if maxlen == 0:
        warnings.warn('template without data')
        return 0.0, len(stream)
    if not stations[0] == 'all':
        for tr in stream:
            if tr.stats.station not in stations:
                stream.remove(tr)  # Remove stations we don't want to use
    for tr in stream:
        if not len(tr.data) == maxlen and not len(tr.data) == 0:
            warnings.warn(tr.stats.station + '.' + tr.stats.channel +
                          ' is not the same length, padding \n' +
                          'Length is ' + str(len(tr.data)) + ' samples')
            pad = np.zeros(maxlen - len(tr.data))
            if tr.stats.starttime.hour == 0:
                tr.data = np.concatenate((pad, tr.data), axis=0)
            else:
                tr.data = np.concatenate((tr.data, pad), axis=0)
        elif len(tr.data) == 0:
            tr.data = np.zeros(maxlen)
    # Clip the data to the set length
    if clip:
        for tr in stream:
            tr.trim(tr.stats.starttime + clip[0], tr.stats.starttime + clip[1])
    _coherence = 0.0
    # Loop through channels and generate a correlation value for each
    # unique cross-channel pairing
    for i in range(len(stream)):
        for j in range(i + 1, len(stream)):
            _coherence += np.abs(normxcorr2(stream[i].data,
                                            stream[j].data))[0][0]
    _coherence = 2 * _coherence / (len(stream) * (len(stream) - 1))
    return _coherence, len(stream)


def brightness(stations, nodes, lags, stream, threshold, thresh_type,
               template_length, template_saveloc, coherence_thresh,
               coherence_stations=['all'], coherence_clip=False,
               gap=2.0, clip_level=100, instance=0, pre_pick=0.2,
               plotsave=True, cores=1, debug=0):
    """
    Calculate the brightness function for a single day.

    Written to calculate the brightness function for a single day of data,
    using moveouts from a 3D travel-time grid.

    .. Note::
        Data in stream must be all of the same length and have the same
        sampling rates, see :func:`eqcorrscan.utils.pre_processing.dayproc`

    :type stations: list
    :param stations:
        List of station names from in the form where stations[i] refers to
        nodes[i][:] and lags[i][:]
    :type nodes: list
    :param nodes:
        List of node points where nodes[i] refers to stations[i] and
        nodes[:][:][0] is latitude in degrees, nodes[:][:][1] is longitude in
        degrees, nodes[:][:][2] is depth in km.
    :type lags: numpy.ndarray
    :param lags:
        Array of arrays where lags[i][:] refers to stations[i]. lags[i][j]
        should be the delay to the nodes[i][j] for stations[i] in seconds.
    :type stream: obspy.core.stream.Stream
    :param stream: Data through which to look for detections.
    :type threshold: float
    :param threshold:
        Threshold value for detection of template within the brightness
        function.
    :type thresh_type: str
    :param thresh_type:
        Either MAD or abs where MAD is the Median Absolute Deviation and abs
        is an absolute brightness.
    :type template_length: float
    :param template_length: Length of template to extract in seconds
    :type template_saveloc: str
    :param template_saveloc: Path of where to save the templates.
    :type coherence_thresh: tuple
    :param coherence_thresh:
            Threshold for removing incoherent peaks in the network response,
            those below this will not be used as templates. Must be in the
            form of (a,b) where the coherence is given by: :math:`a-kchan/b`
            where kchan is the number of channels used to compute the
            coherence.
    :type coherence_stations: list
    :param coherence_stations:
        List of stations to use in the coherence thresholding - defaults to
        `all` which uses all the stations.
    :type coherence_clip: tuple
    :param coherence_clip:
        Start and end in seconds of data to window around, defaults to False,
        which uses all the data given.
    :type gap: float
    :param gap: Minimum inter-event time in seconds for detections.
    :type clip_level: float
    :param clip_level:
        Multiplier applied to the mean deviation of the energy as an upper
        limit, used to remove spikes (earthquakes, lightning, electrical
        spikes) from the energy stack.
    :type instance: int
    :param instance:
        Optional, used for tracking when using a distributed computing system.
    :type pre_pick: float
    :param pre_pick: Seconds before the detection time to include in template
    :type plotsave: bool
    :param plotsave:
        Save or show plots, if `False` will try and show the plots on screen -
        as this is designed for bulk use this is set to `True` to save any
        plots rather than show them if you create them - changes the backend
        of matplotlib, so if is set to `False` you will see NO PLOTS!
    :type cores: int
    :param cores: Number of cores to use, defaults to 1.
    :type debug: int
    :param debug: Debug level from 0-5, higher is more output.

    :return: list of templates as :class:`obspy.core.stream.Stream` objects
    :rtype: list
    """
    if plotsave:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.ioff()
    from eqcorrscan.utils import plotting
    # Check that we actually have the correct stations
    realstations = []
    for station in stations:
        st = stream.select(station=station)
        if st:
            realstations += station
    del st
    stream_copy = stream.copy()
    # Force convert to int16
    for tr in stream_copy:
        # int16 max range is +/- 32767
        if max(abs(tr.data)) > 32767:
            tr.data = 32767 * (tr.data / max(abs(tr.data)))
            # Make sure that the data aren't clipped it they are high gain
            # scale the data
        tr.data = tr.data.astype(np.int16)
    # The internal _node_loop converts energy to int16 too to converse memory,
    # to do this it forces the maximum of a single energy trace to be 500 and
    # normalises to this level - this only works for fewer than 65 channels of
    # data
    if len(stream_copy) > 130:
        raise OverflowError('Too many streams, either re-code and cope with' +
                            'either more memory usage, or less precision, or' +
                            'reduce data volume')
    detections = []
    detect_lags = []
    parallel = True
    plotvar = False
    mem_issue = False
    # Loop through each node in the input
    # Linear run
    print('Computing the energy stacks')
    if not parallel:
        for i in range(0, len(nodes)):
            print(i)
            if not mem_issue:
                j, a = _node_loop(stations, lags[:, i], stream, plot=True)
                if 'energy' not in locals():
                    energy = a
                else:
                    energy = np.concatenate((energy, a), axis=0)
                print('energy: ' + str(np.shape(energy)))
            else:
                j, filename = _node_loop(stations, lags[:, i], stream, i,
                                         mem_issue)
        energy = np.array(energy)
        print(np.shape(energy))
    else:
        # Parallel run
        num_cores = cores
        if num_cores > len(nodes):
            num_cores = len(nodes)
        if num_cores > cpu_count():
            num_cores = cpu_count()
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(_node_loop, args=(stations, lags[:, i],
                                                      stream, i, clip_level,
                                                      mem_issue, instance))
                   for i in range(len(nodes))]
        pool.close()
        if not mem_issue:
            print('Computing the cumulative network response from memory')
            energy = [p.get() for p in results]
            pool.join()
            energy.sort(key=lambda tup: tup[0])
            energy = [node[1] for node in energy]
            energy = np.concatenate(energy, axis=0)
            print(energy.shape)
        else:
            pool.join()
    # Now compute the cumulative network response and then detect possible
    # events
    if not mem_issue:
        print(energy.shape)
        indeces = np.argmax(energy, axis=0)  # Indeces of maximum energy
        print(indeces.shape)
        cum_net_resp = np.array([np.nan] * len(indeces))
        cum_net_resp[0] = energy[indeces[0]][0]
        peak_nodes = [nodes[indeces[0]]]
        for i in range(1, len(indeces)):
            cum_net_resp[i] = energy[indeces[i]][i]
            peak_nodes.append(nodes[indeces[i]])
        del energy, indeces
    else:
        print('Reading the temp files and computing network response')
        node_splits = int(len(nodes) // num_cores)
        indeces = [range(node_splits)]
        for i in range(1, num_cores - 1):
            indeces.append(range(node_splits * i, node_splits * (i + 1)))
        indeces.append(range(node_splits * (i + 1), len(nodes)))
        pool = Pool(processes=num_cores)
        results = [pool.apply_async(_cum_net_resp, args=(indeces[i], instance))
                   for i in range(num_cores)]
        pool.close()
        results = [p.get() for p in results]
        pool.join()
        responses = [result[0] for result in results]
        print(np.shape(responses))
        node_indeces = [result[1] for result in results]
        cum_net_resp = np.array(responses)
        indeces = np.argmax(cum_net_resp, axis=0)
        print(indeces.shape)
        print(cum_net_resp.shape)
        cum_net_resp = np.array([cum_net_resp[indeces[i]][i]
                                 for i in range(len(indeces))])
        peak_nodes = [nodes[node_indeces[indeces[i]][i]]
                      for i in range(len(indeces))]
        del indeces, node_indeces
    if plotvar:
        cum_net_trace = deepcopy(stream[0])
        cum_net_trace.data = cum_net_resp
        cum_net_trace.stats.station = 'NR'
        cum_net_trace.stats.channel = ''
        cum_net_trace.stats.network = 'Z'
        cum_net_trace.stats.location = ''
        cum_net_trace.stats.starttime = stream[0].stats.starttime
        cum_net_trace = Stream(cum_net_trace)
        cum_net_trace += stream.select(channel='*N')
        cum_net_trace += stream.select(channel='*1')
        cum_net_trace.sort(['network', 'station', 'channel'])
        # np.save('cum_net_resp.npy',cum_net_resp)
        #     cum_net_trace.plot(size=(800,600), equal_scale=False,\
        #                        outfile='NR_timeseries.eps')

    # Find detection within this network response
    print('Finding detections in the cumulatve network response')
    detections = _find_detections(cum_net_resp, peak_nodes, threshold,
                                  thresh_type, stream[0].stats.sampling_rate,
                                  realstations, gap)
    del cum_net_resp
    templates = []
    nodesout = []
    good_detections = []
    if detections:
        print('Converting detections into templates')
        # Generate a catalog of detections
        detections_cat = Catalog()
        for j, detection in enumerate(detections):
            if debug > 3:
                print('Converting for detection ' + str(j) + ' of ' +
                      str(len(detections)))
            # Create an event for each detection
            event = Event()
            # Set up some header info for the event
            event.event_descriptions.append(EventDescription())
            event.event_descriptions[0].text = 'Brightness detection'
            event.creation_info = CreationInfo(agency_id='EQcorrscan')
            copy_of_stream = deepcopy(stream_copy)
            # Convert detections to obspy.core.event type -
            # name of detection template is the node.
            node = (detection.template_name.split('_')[0],
                    detection.template_name.split('_')[1],
                    detection.template_name.split('_')[2])
            # Look up node in nodes and find the associated lags
            index = nodes.index((float(node[0]), float(node[1]),
                                 float(node[2])))
            detect_lags = lags[:, index]
            ksta = Comment(text='Number of stations=' + str(len(detect_lags)))
            event.origins.append(Origin())
            event.origins[0].comments.append(ksta)
            event.origins[0].time = copy_of_stream[0].stats.starttime +\
                detect_lags[0] + detection.detect_time
            event.origins[0].latitude = node[0]
            event.origins[0].longitude = node[1]
            event.origins[0].depth = node[2]
            for i, detect_lag in enumerate(detect_lags):
                station = stations[i]
                st = copy_of_stream.select(station=station)
                if len(st) != 0:
                    for tr in st:
                        _waveform_id = WaveformStreamID(station_code=tr.stats.
                                                        station,
                                                        channel_code=tr.stats.
                                                        channel,
                                                        network_code='NA')
                        event.picks.append(Pick(waveform_id=_waveform_id,
                                                time=tr.stats.starttime +
                                                detect_lag +
                                                detection.detect_time +
                                                pre_pick,
                                                onset='emergent',
                                                evalutation_mode='automatic'))
            if debug > 0:
                print('Generating template for detection: ' + str(j))
            template = (_template_gen(event.picks, copy_of_stream,
                        template_length, 'all'))
            template_name = template_saveloc + '/' +\
                str(template[0].stats.starttime) + '.ms'
            # In the interests of RAM conservation we write then read
            # Check coherency here!
            temp_coher, kchan = coherence(template, coherence_stations,
                                          coherence_clip)
            coh_thresh = float(coherence_thresh[0]) - kchan / \
                float(coherence_thresh[1])
            coherent = False
            if temp_coher > coh_thresh:
                template.write(template_name, format="MSEED")
                print('Written template as: ' + template_name)
                print('---------------------------------coherence LEVEL: ' +
                      str(temp_coher))
                coherent = True
            elif debug > 0:
                print('Template was incoherant, coherence level: ' +
                      str(temp_coher))
                coherent = False
            del copy_of_stream, tr, template
            if coherent:
                templates.append(obsread(template_name))
                nodesout += [node]
                good_detections.append(detection)
            elif debug > 0:
                print('No template for you')
    if plotvar:
        all_detections = [(cum_net_trace[-1].stats.starttime +
                           detection.detect_time).datetime
                          for detection in detections]
        good_detections = [(cum_net_trace[-1].stats.starttime +
                            detection.detect_time).datetime
                           for detection in good_detections]
        if not plotsave:
            plotting.NR_plot(cum_net_trace[0:-1],
                             Stream(cum_net_trace[-1]),
                             detections=good_detections,
                             size=(18.5, 10),
                             title='Network response')
            # cum_net_trace.plot(size=(800,600), equal_scale=False)
        else:
            savefile = 'plots/' +\
                cum_net_trace[0].stats.starttime.datetime.strftime('%Y%m%d') +\
                '_NR_timeseries.pdf'
            plotting.NR_plot(cum_net_trace[0:-1],
                             Stream(cum_net_trace[-1]),
                             detections=good_detections,
                             size=(18.5, 10), save=True, savefile=savefile,
                             title='Network response')
    nodesout = list(set(nodesout))
    return templates, nodesout


if __name__ == "__main__":
    import doctest
    doctest.testmod()
