#!/usr/bin/python
"""
Utility module of the EQcorrscan package to allow for different methods of
stacking of seismic signal in one place.

In alpha stages and only with linear stacking implimented thusfar

Calum Chamberlain 24/06/2015

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
def linstack(streams):
    """
    Function to compute the linear stack of a series of seismic streams of
    multiplexed data

    :type streams: List of Streams
    :param stream: List of streams to stack

    :returns: stack - Stream
    """
    # import matplotlib.pyplot as plt
    stack=streams[np.argmax([len(stream) for stream in streams])].copy()
    for tr in stack:
        tr.data=tr.data/np.sqrt(np.mean(np.square(tr.data)))
        tr.data=np.nan_to_num(tr.data)
    for i in xrange(1,len(streams)):
        # print "Stacking stream "+str(i)
        for tr in stack:
            # print tr.stats.station+'.'+tr.stats.channel
            matchtr=streams[i].select(station=tr.stats.station,\
                                       channel=tr.stats.channel)
            if matchtr:
                norm=matchtr[0].data/np.sqrt(np.mean(np.square(matchtr[0].data)))
                norm=np.nan_to_num(norm)
                tr.data=np.sum((norm,\
                               tr.data), axis=0)
    return stack

def PWS_stack(streams, weight=2):
    """
    Function to compute the phase weighted stack of a series of streams.
    Recommend aligning the traces before stacking.

    :type streams: list of obspy.Stream
    :param streams: List of Stream to stack
    :type weight: float
    :param weight: Exponent to the phase stack used for weighting.

    :return: obspy.Stream
    """
    from scipy.signal import hilbert
    # First get the linear stack which we will weight by the phase stack
    Linstack=linstack(streams)
    # Compute the instantaneous phase
    instaphases=[]
    print "Computing instantaneous phase"
    for stream in streams:
        instaphase=stream.copy()
        for tr in instaphase:
            analytic=hilbert(tr.data)
            envelope=np.sqrt(np.sum((np.square(analytic),\
                                     np.square(tr.data)), axis=0))
            tr.data=analytic/envelope
        instaphases.append(instaphase)
    # Compute the phase stack
    print "Computing the phase stack"
    Phasestack=linstack(instaphases)
    # print type(Phasestack)
    # Compute the phase-weighted stack
    for tr in Phasestack:
        tr.data=Linstack.select(station=tr.stats.station)[0].data*\
                np.abs(tr.data**weight)
    return Phasestack

def align_traces(trace_list, shift_len):
    """
    Function to allign traces relative to each other based on their
    cross-correlation value

    :type trace_list: List of Traces
    :param trace_list: List of traces to allign
    :type shift_len: int
    :param shift_len: Length to allow shifting within in samples

    :returns: list of shifts for best allignment in seconds
    """
    from obspy.signal.cross_correlation import xcorr
    from copy import deepcopy
    traces=deepcopy(trace_list)
    # Use trace with largest MAD amplitude as master
    master=traces[0]
    MAD_master=np.median(np.abs(master.data))
    master_no=0
    for i in xrange(1,len(traces)):
        if np.median(np.abs(traces[i])) > MAD_master:
            master=traces[i]
            MAD_master=np.median(np.abs(master.data))
            master_no=i
    shifts=[]
    for i in xrange(len(traces)):
        if not master.stats.sampling_rate == traces[i].stats.sampling_rate:
            raise ValueError('Sampling rates not the same')
        shift, cc=xcorr(master, traces[i], shift_len)
        shifts.append(shift/master.stats.sampling_rate)
    return shifts
