#!/usr/bin/python
"""
Part of the EQcorrscan package written by Calum Chamberlain of Victoria
University of Wellington.

Functions designed to test the reosultion of templat matching cross-correlation
function inspired by Emily Warren-Smith.
"""

import sys, os
bob=os.path.realpath(__file__)
bob=bob.split('/')
path='/'
for i in xrange(len(bob)-2):
    path+=bob[i]+'/'
print path
sys.path.insert(0,path)
import numpy as np
from match_filter import normxcorr2
from utils import findpeaks

def ccc_gen(image, template):
    """
    Function to test if a detection is possible at this moveout

    :type image: obspy.Stream
    :type template: obspy.Stream
    :type delays: list
    :type threhsold: float

    :returns: ccc, a matrix of correlation values
    """
    ccc=np.array([np.array([0]*(len(image[0].data)-len(template[0].data)+1))]*len(image))
    print 'Correlation matrix is shaped: '+str(ccc.shape)
    for i in xrange(len(image)):
        templatetr=template.select(station=image[i].stats.station,\
                                   channel=image[i].stats.channel)
        ccc[i]=normxcorr2(templatetr[0].data, image[i].data)[0]
    return ccc

def _node_loop(freq, node, ccc, lags, threshold, thresh_type):
    """
    """
    for j in xrange(len(ccc)):
        pad=np.array([0]*int(round(lags[j,0]*freq)))
        if not 'cccsum' in locals():
            cccsum=np.append(ccc[j],pad)[len(pad):]
        else:
            bob=np.append(ccc[j],pad)[len(pad):]
            cccsum=np.sum([cccsum, bob], axis=0)
    if thresh_type=='MAD':
        rawthresh=threshold*np.median(np.abs(cccsum))
    elif thresh_type=='absolute':
        rawthresh=threshold
    elif thresh_type=='av_chan_corr':
        rawthresh=threshold*(cccsum/len(image))
    else:
        rawthresh=threshold*np.median(np.abs(cccsum))
    peaks = findpeaks.find_peaks2(cccsum, rawthresh,\
                                  freq, 0)
    if not len(peaks) == 0:
        return node
    else:
        return

def moveout_check(template, nodes, lags, threshold, thresh_type, lowcut,\
                  highcut, filt_order):
    """
    Function to check different moveouts for detectability
    """
    # Generate random noise and seed with the template
    from copy import deepcopy
    from obspy.signal.filter import bandpass
    from joblib import Parallel, delayed
    parallel=True
    image=deepcopy(template)
    for i in xrange(len(image)):
        image[i].data=np.random.randn(len(image[i].data)+\
                                      (86400*image[i].stats.sampling_rate)\
                                      -len(image[i].data))
        image[i].data[(len(image[i].data)-len(template[i].data))/2:\
                      (len(image[i].data)-len(template[i].data))/2+len(template[i].data)]=\
                image[i].data[(len(image[i].data)-len(template[i].data))/2:\
                      (len(image[i].data)-len(template[i].data))/2+len(template[i].data)]+\
                template[i].data/np.mean(template[i].data**2)**0.5
        image[i].data=bandpass(image[i].data, lowcut, highcut,\
                               image[i].stats.sampling_rate, filt_order)
    ccc=ccc_gen(image, template)
    # Lags
    possible_locations=[]
    freq=image[0].stats.sampling_rate
    if not parallel:
        for i in xrange(len(nodes)):
            possible_locations+=_node_loop(freq, nodes[i], ccc, lags[:,[i]],\
                                           threshold, thresh_type)
    else:
        possible_locations = Parallel(n_jobs=10, verbose=5)(delayed(_node_loop)\
                                                            (freq, nodes[i], ccc,\
                                                             lags[:,[i]], threshold,\
                                                             thresh_type)\
                                                            for i in xrange(len(nodes)))
    return possible_locations

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        raise IOError("Needs one argument, the template to check")
    templatename=str(sys.argv[1])
    from par import match_filter_par as defaults
    from par import template_gen_par as tempdef
    from par import bright_lights_par as brightdef
    from obspy import read
    from core.bright_lights import _read_tt, _resample_grid
    template=read(templatename)
    stations_in=[]
    for tr in template:
        stations_in+=[tr.stats.station]
    stations_in=list(set(stations_in))
    stations, nodes, lags = _read_tt(brightdef.nllpath,stations_in,\
                                    brightdef.phase, phaseout='S', \
                                    ps_ratio=brightdef.ps_ratio)
    stations, nodes, lags = _resample_grid(stations, nodes, lags, brightdef.volume,\
                                           brightdef.resolution)
    # print np.shape(lags)
    for station in stations:
        if not 'template_dummy' in locals():
            template_dummy=template.select(station=station)
        else:
            template_dummy+=template.select(station=station)
    template=template_dummy
    for tr in template:
        for i in xrange(len(stations)):
            if tr.stats.station == stations[i]:
                if not 'alllags' in locals():
                    alllags=[lags[i]]
                else:
                    # print stations[i]
                    alllags=np.concatenate((alllags, [lags[i]]), axis=0)
                    # print np.shape(alllags)
    lags=alllags
    print 'Lags is shaped: '+str(np.shape(lags))
    print 'I have '+str(len(template))+' channels of data'
# Indexing will be an issue, currently don't check that stations match between data and lags
    possible_locations = moveout_check(template, nodes, lags, defaults.threshold,\
                                       defaults.threshtype, tempdef.lowcut,\
                                       tempdef.highcut, tempdef.filter_order)
    from utils import EQcorrscan_plotting as plotting
    if not len(possible_locations) == 0:
        plotting.threeD_gridplot(possible_locations)
    else:
        raise ValueError("No possible location found")
