#!/usr/bin/python
"""
Code to compute the linkage between seismograms and cluster them accordingly

Written by Calum Chamberlain, in alpha stages of development as of 24/06/2015

Implimented to streamline templates after template detection in beamforming
methods, employed by implimentation of Frank et al. code.

As such this code is designed to work only for templates with the same channels

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

# I want to make sure that the lags for those that I cluster are similar as well
# as making sure that they are similarly correlated - will cluster based on
# cross-channel correlation sum
import numpy as np

def cross_chan_coherence(st1, st2):
    """
    Function to determine the cross-channel coherancy between two streams of
    multichannel seismic data.

    :type st1: obspy Stream
    :param st1: Stream one
    :type st2: obspy Stream
    :param st2: Stream two

    :returns: cross channel coherence, float - normalized by number of channels
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
    :param template: List of the streams to compute the distance matrix for

    :returns: ndarray - distance matrix
    """
    # Initialize square matrix
    dist_mat=np.array([np.array([0.0]*len(templates))]*len(templates))
    for i in xrange(len(templates)):
        for j in xrange(i,len(templates)):
            if i==j:
                dist_mat[i,j]=0.0
            else:
                dist_mat[i,j]=1-np.abs(cross_chan_coherence(templates[i],templates[j]))
    for i in xrange(1,len(templates)):
        for j in xrange(i):
            dist_mat[i,j]=dist_mat.T[i,j]
    return dist_mat

def cluster(templates, show=True):
    """
    Function to take a set of templates and cluster them, will return distance\
            matrix of grouped templates

    :type templates: List of Obspy.Stream
    :param templates: List of templates to compute clustering for
    :type show: bool
    :param show: plot linkage on screen if True, defaults to True

    :returns: List of cluster groups, array of length len(templates), with
                each number relating to a cluster
    .. rubric:: Note:
        Not fully feautured yet, returns the Z matrix, but doesn't tell you what\
        can be clustered.
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

    :returns: List of List of obspy.Streams where each initial list is a group\
            with the same delays
    """
    groups=[]
    group_delays=[]
    group_chans=[]
    # Sort templates by number of channels
    templates=[(template, len(template)) for template in templates]
    templates.sort(key=lambda tup:tup[1])
    templates=[template[0] for template in templates]
    for i in xrange(1,len(templates)):
        print 'Working on waveform '+str(i)+' of '+str(len(templates))
        # Calculate the delays
        starttimes=[]
        chans=[]
        for tr in templates[i]:
            starttimes.append(tr.stats.starttime)
            chans.append((tr.stats.station, tr.stats.channel))
        # This delays calculation will be an issue if we have changes in channels
        delays=[starttimes[m]-min(starttimes) for m in xrange(len(starttimes))]
        delays=[round(d,2) for d in delays]
        if len(groups)==0:
            groups.append([templates[i]])
            group_delays.append(delays)
            group_chans.append(chans)
        else:
            j=0
            match=False
            while not match:
                kmatch=0
                # Find the set of shared stations and channels
                shared_chans=[]
                shared_delays_slave=[]
                shared_delays_master=[]
                k=0
                for chan in chans:
                    if chan in group_chans[j]:
                        shared_chans.append(chan)
                        shared_delays_slave.append(delays[k])
                        shared_delays_master.append(group_delays[j][group_chans[j].index(chan)])
                    k+=1
                # Normalize master and slave delay times
                shared_delays_slave=[delay-min(shared_delays_slave)\
                                     for delay in shared_delays_slave]
                shared_delays_master=[delay-min(shared_delays_master)\
                                      for delay in shared_delays_master]
                for k in xrange(len(shared_chans)):
                    # Check if the channel and delay match another group
                    if shared_delays_slave[k]==shared_delays_master[k]:
                        kmatch+=1 # increase the match index
                if kmatch==len(shared_chans): # If all the channels match, add it to the group
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

def SVD(templates):
    """
    Function to compute the SVD of a number of templates and return the singular
    vectors and singular values of the templates.

    :type templates: List of Obspy.Stream
    :param templates: List of the templates to be analysed

    :return: SVector(list of ndarray), SValues(list) for each channel, \
            Uvalues(list of ndarray) for each channel, \
            stachans, List of String (station.channel)

    .. rubric:: Note

    It is recommended that you align the data before computing the SVD, e.g.,
    the P-arrival on all templates for the same channel should appear at the same
    time in the trace.
    """
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
    Uvectors=[]
    for stachan in stachans:
        chan_mat=[templates[i].select(station=stachan.split('.')[0], \
                                  channel=stachan.split('.')[1])[0].data \
                  for i in xrange(len(templates)) if \
                  len(templates[i].select(station=stachan.split('.')[0], \
                                  channel=stachan.split('.')[1])) != 0]
        # chan_mat=[chan_mat[i]/np.max(chan_mat[i]) for i in xrange(len(chan_mat))]
        chan_mat=np.asarray(chan_mat)
        print chan_mat.shape
        print stachan
        U, s, V = np.linalg.svd(chan_mat, full_matrices=False)
        SValues.append(s)
        SVectors.append(V)
        Uvectors.append(U)
    return SVectors, SValues, Uvectors, stachans

def empirical_SVD(templates, linear=True):
    """
    Empirical subspace detector generation function.  Takes a list of templates
    and computes the stack as the first order subspace detector, and the
    differential of this as the second order subspace detector following
    the emprical subspace method of Barrett & Beroza, 2014 - SRL.

    :type templates: list of stream
    :param templates: list of template streams to compute the subspace detectors\
        from
    :type linear: Bool
    :param linear: Set to true by default to compute the linear stack as the\
        first subspace vector, False will use the phase-weighted stack as the\
        first subspace vector.

    :returns: list of two streams
    """
    import stacking
    if linear:
        first_subspace=stacking.linstack(templates)
    second_subspace=first_subspace.copy()
    for i in xrange(len(second_subspace)):
        second_subspace[i].data=np.diff(second_subspace[i].data)
        second_subspace[i].stats.starttime+=0.5*second_subspace[i].stats.delta
    return [first_subspace, second_subspace]

def SVD_2_stream(SVectors, stachans, k, sampling_rate):
    """
    Function to convert the singular vectors output by SVD to streams, one for
    each singular vector level, for all channels.

    :type SVectors: List of np.ndarray
    :param SVectors: Singular vectors
    :type stachans: List of Strings
    :param stachans: List of station.channel Strings
    :type k: int
    :param k: Number of streams to return = number of SV's to include
    :type sampling_rate: float
    :param sampling_rate: Sampling rate in Hz

    :returns: SVstreams, List of Obspy.Stream, with SVStreams[0] being
            composed of the highest rank singular vectors.
    """
    from obspy import Stream, Trace
    SVstreams=[]
    for i in xrange(k):
        SVstream=[]
        for j in xrange(len(stachans)):
            SVstream.append(Trace(SVectors[j][i], \
                                    header={'station': stachans[j].split('.')[0],
                                            'channel': stachans[j].split('.')[1],
                                            'sampling_rate': sampling_rate}))

        SVstreams.append(Stream(SVstream))
    return SVstreams

def corr_cluster(traces, thresh=0.9):
    """
    Group traces based on correlations above threshold with the stack - will
    run twice, once with a lower threshold, then again with your threshold to
    remove large outliers

    :type traces: List of :class:obspy.Trace
    :param traces: Traces to compute similarity between
    :type thresh: float
    :param thrsh: Correlation threshold between -1-1

    :returns: np.ndarray of bool
    """
    import stacking
    from obspy import Stream
    from core.match_filter import normxcorr2
    stack=stacking.linstack([Stream(tr) for tr in traces])[0]
    output=np.array([False]*len(traces))
    group1=[]
    for i in xrange(len(traces)):
        if normxcorr2(traces[i].data,stack.data)[0][0] > 0.6:
            output[i]=True
            group1.append(traces[i])
    stack=stacking.linstack([Stream(tr) for tr in group1])[0]
    group2=[]
    for i in xrange(len(traces)):
        if normxcorr2(traces[i].data,stack.data)[0][0] > thresh:
            group2.append(tr)
            output[i]=True
        else:
            output[i]=False
    return output

def extract_detections(detections, templates, extract_len=90.0, \
                        outdir=None, extract_Z=True, additional_stations=[]):
    """
    Function to extract the waveforms associated with each detection in a list
    of detections for the template, template.  Waveforms will be returned as
    a list of obspy.Streams containing segments of extract_len.  They will also
    be saved if outdir is set.  The default is unset.  The default extract_len
    is 90 seconds per channel.

    :type detections: List tuple of of :class: datetime.datetime, string
    :param detections: List of datetime objects, and their associated template\
            name
    :type templates: List of tuple of string and :class: obspy.Stream
    :param templates: A list of the tuples of the template name and the template\
            Stream used to detect detections.
    :type extract_len: float
    :param extract_len: Length to extract around the detection (will be equally\
            cut around the detection time) in seconds.  Default is 90.0.
    :type outdir: Bool or String
    :param outdir: Default is None, with None set, no files will be saved,\
            if set each detection will be saved into this directory with files\
            named according to the detection time, NOT than the waveform\
            start time. Detections will be saved into template subdirectories.
    :type extract_Z: Bool
    :param extract_Z: Set to True to also extract Z channels for detections\
            delays will be the same as horizontal channels, only applies if\
            only horizontal channels were used in the template.
    :type additional_stations: List of tuple
    :param additional_stations: List of stations, chanels and networks to also\
            extract data for using an average delay.

    :returns: List of :class: obspy.Stream
    """
    from obspy import read, UTCDateTime, Stream
    from utils import pre_processing
    import datetime as dt
    from par import match_filter_par as matchdef
    from par import template_gen_par as templatedef
    import os
    from joblib import Parallel, delayed
    # Sort the template according to starttimes, needed so that stachan[i]
    # corresponds to delays[i]
    all_delays=[] # List of tuples of template name, delays
    all_stachans=[]
    for template in templates:
        templatestream=template[1].sort(['starttime'])
        stachans=[(tr.stats.station,tr.stats.channel,tr.stats.network) \
                  for tr in templatestream]
        mintime=templatestream[0].stats.starttime
        delays=[tr.stats.starttime-mintime for tr in templatestream]
        all_delays.append((template[0], delays))
        all_stachans.append((template[0], stachans))
    # Sort the detections and group by day
    detections.sort()
    detection_days=[detection[0].date() for detection in detections]
    detection_days=list(set(detection_days))
    detection_days.sort()

    # Initialize output list
    detection_wavefiles=[]

    # Also include Z channels when extracting detections
    if extract_Z:
        new_all_stachans=[]
        new_all_delays=[]
        t=0
        for template in all_stachans:
            stachans=template[1]
            delays=all_delays[t][1]
            new_stachans=[]
            new_delays=[]
            j=0
            for i in xrange(len(stachans)):
                if j==1:
                    new_stachans.append((stachans[i][0], stachans[i][1][0]+'Z',\
                                         stachans[i][2]))
                    new_delays.append(delays[i])
                    new_stachans.append(stachans[i])
                    new_delays.append(delays[i])
                    j=0
                else:
                    new_stachans.append(stachans[i])
                    new_delays.append(delays[i])
                    j+=1
            new_all_stachans.append((template[0], new_stachans))
            new_all_delays.append((template[0], new_delays))
            t+=1
        all_delays=new_all_delays
        all_stachans=new_all_stachans
    if not len(additional_stations)==0:
        t=0
        for template in all_stachans:
            av_delay=np.mean(all_delays[t][1])
            for sta in additional_stations:
                if not sta in template[1]:
                    template[1].append(sta)
                    all_delays[t][1].append(av_delay)
            t+=1

    # Loop through the days
    for detection_day in detection_days:
        print 'Working on detections for day: '+str(detection_day)
        stachans=list(set([stachans[1] for stachans in all_stachans][0]))
        # List of all unique stachans - read in all data
        for stachan in stachans:
            contbase=[base for base in matchdef.contbase\
                      if base[2]==stachan[2]][0]
            if contbase[1]=='yyyymmdd':
                dayfile=detection_day.strftime('%Y%m%d')+'/*'+stachan[0]+\
                        '.'+stachan[1][0]+'?'+stachan[1][1]+'.*'
            elif contbase[1]=='Yyyyy/Rjjj.01':
                dayfile=detection_day.strftime('Y%Y/R%j.01')+'/'+stachan[0]+\
                        '.*.'+stachan[1][0]+'?'+stachan[1][1]+'.'+detection_day.strftime('%Y.%j')
            if not 'st' in locals():
                try:
                    st=read(contbase[0]+'/'+dayfile)
                except:
                    print 'No data for '+contbase[0]+'/'+dayfile
            else:
                try:
                    st+=read(contbase[0]+'/'+dayfile)
                except:
                    print 'No data for '+contbase[0]+'/'+dayfile
        st.merge(fill_value='interpolate')
        # We now have a stream of day long data, we should process it!
        # st=Parallel(n_jobs=3)(delayed(pre_processing.dayproc)(tr, templatedef.lowcut,\
                                                               # templatedef.highcut,\
                                                               # templatedef.filter_order,\
                                                               # templatedef.samp_rate,\
                                                               # matchdef.debug, detection_day)\
                                # for tr in st)
        # st=Stream(st)
        # for tr in st:
            # Convert to int32 for STEIM2 format
            # tr.data=tr.data.astype(np.int32)

        day_detections=[detection for detection in detections\
                        if detection[0].date() == detection_day]
        for detection in day_detections:
            template=detection[1]
            t_stachans=[stachans[1] for stachans in all_stachans \
                      if stachans[0] == template][0]
            t_delays=[delays[1] for delays in all_delays\
                    if delays[0] == template][0]
            print 'Cutting for detections at: '+detection[0].strftime('%Y/%m/%d %H:%M:%S')
            detect_wav=st.copy()
            for tr in detect_wav:
                tr.trim(starttime=UTCDateTime(detection[0])-extract_len/2,\
                            endtime=UTCDateTime(detection[0])+extract_len/2)
            if outdir:
                if not os.path.isdir(outdir+'/'+template):
                    os.makedirs(outdir+'/'+template)
                detect_wav.write(outdir+'/'+template+'/'+\
                                  detection[0].strftime('%Y-%m-%d_%H-%M-%S')+\
                                  '.ms', format='MSEED', encoding='STEIM2')
                                  # '.ms', format='MSEED', encoding='STEIM2')
                print 'Written file: '+outdir+'/'+template+'/'+\
                         detection[0].strftime('%Y-%m-%d_%H-%M-%S')+'.ms'
            if not outdir:
                detection_wavefiles.append(detect_wav)
            del detect_wav
        del st
        if outdir:
            detection_wavefiles=[]
    if not outdir:
        return detection_wavefiles
    else:
        return
