#!/usr/bin/python
"""
Code to compute the linkage between seismograms and cluster them accordingly

Written by Calum Chamberlain, in alpha stages of development as of 24/06/2015

Implimented to streamline templates after template detection in beamforming
methods, employed by implimentation of Frank et al. code.

As such this code is designed to work only for templates with the same channels
"""

# I want to make sure that the lags for those that I cluster are similar as well
# as making sure that they are similarly correlated - will cluster based on
# cross-channel correlation sum

def cross_chan_coherance(st1, st2):
    """
    Function to determine the cross-channel coherancy between two streams of
    multichannel seismic data.

    :type st1: obspy Stream
    :type st2: obspy Stream

    :returns: cross channel coherance, float - normalized by number of channels
    """
    from core.match_filter import normxcorr2
    cccoh=0.0
    kchan=0
    for tr in st1:
        tr1=tr.data
        # Assume you only have one waveform for each channel
        tr2=st2.select(station=tr.stats.station, \
                       channel=tr.stats.channel)
        if tr2:
            cccoh+=normxcorr2(tr1,tr2[0].data)[0][0]
            kchan+=1
    cccoh=cccoh/kchan
    return cccoh

def distance_matrix(templates):
    """
    Function to compute the distance matrix for all templates - will give
    distance as 1-abs(cccoh), e.g. a well correlated pair of templates will
    have small distances, and an equally well correlated reverse image will
    have the same distance as apositively correlated image - this is an issue

    :type templates: List of obspy.Streams

    :returns: ndarray - distance matrix
    """
    import numpy as np
    # Initialize square matrix
    dist_mat=np.array([np.array([0.0]*len(templates))]*len(templates))
    for i in xrange(len(templates)):
        for j in xrange(i,len(templates)):
            if i==j:
                dist_mat[i,j]=0.0
            else:
                dist_mat[i,j]=1-np.abs(cross_chan_coherance(templates[i],templates[j]))
    for i in xrange(1,len(templates)):
        for j in xrange(i):
            dist_mat[i,j]=dist_mat.T[i,j]
    return dist_mat

def cluster(templates, show=True):
    """
    Function to take a set of templates and cluster them, will return clustered
    templates

    :type template: List of Obspy.Stream

    :returns: List of cluster groups, array of length len(templates), with
                each number relating to a cluster
    """
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, dendrogram
    import matplotlib.pyplot as plt
    dist_mat=distance_matrix(templates)
    dist_vec=squareform(dist_mat)
    # plt.matshow(dist_mat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnB)
    Z = linkage(dist_vec)
    D = dendrogram(Z)
    if show:
        plt.show()
    return Z

def group_delays(templates):
    """
    Function to group template waveforms according to their delays

    :type templates: List of obspy.Stream
    :param templates: List of the waveforms you want to group

    :returns: List of List of obspy.Streams where each initial list is a group
            with the same delays
    """
    groups=[]
    group_delays=[]
    group_chans=[]
    for i in xrange(1,len(templates)):
        print 'Working on waveform '+str(i)+' of '+str(len(templates))
        # Calculate the delays
        starttimes=[]
        chans=[]
        for tr in templates[i]:
            starttimes.append(tr.stats.starttime)
            chans.append((tr.stats.station, tr.stats.channel))
        delays=[starttimes[m]-min(starttimes) for m in xrange(len(starttimes))]
        if len(groups)==0:
            groups.append([templates[i]])
            group_delays.append(delays)
            group_chans.append(chans)
        else:
            j=0
            match=False
            while not match:
                kmatch=0
                for k in xrange(len(chans)):
                    # Check if the channel and delay match another group
                    if chans[k] in group_chans[j] and \
                       delays[k]==group_delays[j][group_chans[j].index(chans[k])]:
                        kmatch+=1 # increase the match index
                if kmatch==len(chans): # If all the channels match, add it to the group
                    groups[j].append(templates[i])
                    match=True
                elif j<len(groups)-1:
                    j+=1
                else:
                    # Create a new group and break the loop
                    groups.append([templates[i]])
                    group_delays.append(delays)
                    group_chans.append(chans)
                    match=True # Use this to break the loop
    return groups

# def generate_families(templates):
    # """
    # Function to take a load of templates, group them according to their
    # delays, then within these groups group them according to similarity
    # and stack the results into tempates

    # :type templates: List of obspy.Stream
    # :param templates: List of all possible templates to use

    # :returns: List of templates (stacked)
    # """
    # import stacking
    # groups=group_delays(templates)
    # for group in groups:
        # Z=cluster(group, show=False)
    # #### UNFINISHED!!!
    # return stacked_templates

def allign_traces(traces):
    """
    Function to allign all traces in a stream based on their cross-correlations

    :type traces: List of :class: obspy.Trace
    :param traces: List of of single channels unalligned - should contain\
            data from multiple earthquakes for a single site, channel and phase\
            I will compute the linear stack of all traces and allign the traces\
            relative to this.

    :returns: List of ndarray
    """
    # Set the master trace as the linear stack of all traces in the stream
    from stacking import linstack
    from obspy import Stream, Trace
    from obspy.signal.cross_correlation import xcorr
    import numpy as np
    import warnings, scipy
    # Check the type, linstack needs a list of streams
    for i in xrange(len(traces)):
        if type(traces[i])==Trace:
            traces[i]=Stream(traces[i])
        elif type(traces[i]) != Stream:
            raise IOError('You have not passed me an obspy stream')
    master=linstack(traces)
    alligned_data=[]
    for trace in traces:
        shift_len=int(0.25*trace[0].stats.npts)
        # Mostly copied from xcorrPickCorrection
        _cc_shift, cc_max, cc = xcorr(master[0].data, trace[0].data, shift_len,\
                                      full_xcorr=True)
        print cc_max, _cc_shift
        # make array with time shifts in seconds corresponding to cc function
        correction = _cc_shift
        print correction
        if not abs(correction) < 1.0:
            continue
        pad=np.zeros(abs(correction)*trace[0].stats.sampling_rate)
        if len(pad) < len(trace[0].data):
            if correction > 0:
                alligned_data.append(np.append(pad, trace[0].data[0:-len(pad)]))
            elif correction < 0:
                alligned_data.append(np.append(trace[0].data[len(pad),:], pad))
            else:
                alligned_data.append(trace[0].data)
        else:
            alligned_data.append(trace[0].data)
    return alligned_data

def SVD_testing(templates):
    """
    Function to compute the SVD of a number of templates and return the singular
    vectors and singular values of the templates.

    :type templates: List of Obspy.Stream
    :param templates: List of the templates to be analysed

    :return: SVector(ndarray), SValues(ndarray) for each channel, stachans, List
            of String (station.channel)

    .. rubric:: Note

    **IN ALPHA, not working as expected**

    It is recommended that you align the data before computing the SVD, .e.g.,
    the P-arrival on all templates for the same channel should appear at the same
    time in the trace.
    """
    import numpy as np
    # Convert templates into ndarrays for each channel
    # First find all unique channels:
    stachans=[]
    for template in templates:
        for tr in template:
            stachans.append(tr.stats.station+'.'+tr.stats.channel)
    stachans=list(set(stachans))
    # Initialize a list for the output matrices, one matrix per-channel
    SValues=[]
    SVectors=[]
    for stachan in stachans:
        chan_mat=[templates[i].select(station=stachan.split('.')[0], \
                                  channel=stachan.split('.')[1])[0].data \
                  for i in xrange(len(templates)) if \
                  len(templates[i].select(station=stachan.split('.')[0], \
                                  channel=stachan.split('.')[1])) != 0]
        chan_mat=[chan_mat[i]/np.max(chan_mat[i]) for i in xrange(len(chan_mat))]
        chan_mat=np.asarray(chan_mat).T
        print chan_mat.shape
        print stachan
        U, s, V = np.linalg.svd(chan_mat, full_matrices=False)
        SValues.append(s)
        SVectors.append(U.T)
    return SVectors, SValues, stachans

def SVD_2_stream_testing(SVectors, stachans, k, sampling_rate):
    """
    Function to convert the singular vectors output by SVD to streams, one for
    each singular vector level, for all channels.

    :type SVectors: ndarray
    :param SVectors: Singular vectors
    :type stachans: List of Strings
    :param stachans: List of station.channel Strings
    :type k: Int
    :param k: Number of streams to return = number of SV's to include
    :sampling_rate: Float

    :returns: SVstreams, List of Obspy.Stream, with SVStreams[0] being
            composed of the highest rank singular vectors.

    .. rubric:: Note

    **IN ALPHA, not working as expected**
    """
    from obspy import Stream, Trace
    SVstreams=[]
    for i in xrange(k):
        SVstream=[]
        for j in xrange(len(stachans)):
            if len(SVectors[i]) > j:
                SVstream.append(Trace(SVectors[i][j], \
                                        header={'station': stachans[j].split('.')[0],
                                                'channel': stachans[j].split('.')[1],
                                                'sampling_rate': sampling_rate}))
        SVstreams.append(Stream(SVstream))
    return SVstreams
