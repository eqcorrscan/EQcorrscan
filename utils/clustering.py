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
    kchan=len(st1)
    for tr in st1:
        tr1=tr.data
        # Assume you only have one waveform for each channel
        tr2=st2.select(station=tr.stats.station, \
                       channel=tr.stats.channel)[0].data
        cccoh+=normxcorr2(tr1,tr2)[0][0]
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
        for j in xrange(len(templates)):
            if i==j:
                dist_mat[i,j]=0.0
            else:
                dist_mat[i,j]=1-np.abs(cross_chan_coherance(templates[i],templates[j]))
    return dist_mat

def cluster(templates):
    """
    Function to take a set of templates and cluster them, will return clustered
    templates

    :type template: List of Obspy.Stream

    :returns: List of cluster groups, array of length len(templates), with
                each number relating to a cluster
    """
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage
    dist_mat=distance_matrix(templates)
    dist_vec=squareform(dist_mat)
    Z = linkage(dist_vec)
