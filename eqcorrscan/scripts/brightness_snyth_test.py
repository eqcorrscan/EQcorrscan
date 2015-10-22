#!/usr/bin/env python
"""
Script to test the brightness function with synthetic data
"""
import pickle
import sys, glob, datetime as dt
import matplotlib.pyplot as plt
sys.path.append('/Volumes/GeoPhysics_09/users-data/chambeca/my_programs/Building/EQcorrscan')
from par import LFE_template_gen_par as templatedef
from par import match_filter_par_LFEs as matchdef
from par import bright_lights_par as brightdef
from obspy import UTCDateTime, Stream, Trace, read as obsread
# First generate the templates
from core import bright_lights, match_filter
from utils import pre_processing
from utils import EQcorrscan_plotting as plotting
from obspy.signal.filter import bandpass
from joblib import Parallel, delayed
import warnings, pickle, copy
import numpy as np

# Read in the nodes
with open('Nodes_in.pkl','rb') as pickle_load:
    nodes = pickle.load(pickle_load)
with open('Lags_in.pkl','rb') as pickle_load:
    lags = pickle.load(pickle_load)
with open('Stations_in.pkl','rb') as pickle_load:
    stations = pickle.load(pickle_load)
print 'Read in '+str(len(nodes))+' nodes'

# Generate the synthetic data with a spike at nodes[0]
flat=np.zeros(1500)
samp_rate=20 # 4 Hz
offset=int(10*samp_rate)
length=int(0.5*samp_rate)
ksta=0
knode=0
realstr=True
starttime=UTCDateTime('20140926')
for station in stations:
    if not 'stream' in locals():
        spiked=copy.deepcopy(flat)
        spiked[offset+(samp_rate*lags[ksta][knode]):\
               offset+length+(samp_rate*lags[ksta][knode])]=1
        tr=Trace(spiked)
        tr.stats.station=station
        tr.stats.channel='S1'
        tr.stats.network='SYN'
        tr.stats.sampling_rate=samp_rate
        tr.stats.starttime=starttime
        stream=Stream(tr)
    else:
        spiked=copy.deepcopy(flat)
        spiked[offset+(samp_rate*lags[ksta][knode]):\
               offset+length+(samp_rate*lags[ksta][knode])]=1
        tr=Trace(spiked)
        tr.stats.station=station
        tr.stats.channel='S1'
        tr.stats.network='SYN'
        tr.stats.sampling_rate=samp_rate
        tr.stats.starttime=starttime
        stream+=tr
    ksta+=1
if realstr:
    # stream=obsread('scripts/brightness_test.ms')
    # stream.detrend('demean')
    # stream=obsread('/Volumes/GeoPhysics_09/users-data/chambeca/SAMBA_archive/day_volumes_S/'+\
                # 'Y2011/R247.01/*N.2011.247')
    # stream.detrend('demean')
    # stream.resample(samp_rate)
    # stream.write('scripts/brightness_test_daylong.ms',format='MSEED')
    stream=obsread('scripts/brightness_test_daylong.ms')
    stream.trim(starttime=UTCDateTime('2011-09-04 17:05:00'),\
                endtime=UTCDateTime('2011-09-04 17:15:00'))#, pad=True,\
               # fill_value=0)
    # for tr in stream:
        # if tr.stats.station=='WVZ':
            # stream.remove(tr)
stream.filter('bandpass',freqmin=4.0, freqmax=8.0)
# stream.trim(stream[0].stats.starttime+90, stream[0].stats.endtime)
stream.trim(stream[0].stats.starttime, stream[0].stats.endtime, pad=True, fill_value=0)
stream.plot(size=(800,600),equal_scale=False)

instance=0

# Cut the nodes...
cutnodes=[nodes[0]]+[nodes[116]]
cutlags=np.array([lags[:,0]]+[lags[:,116]]).T
detect_templates, detect_nodes=bright_lights.brightness(stations, \
                        nodes, lags, stream,
                        brightdef.threshold, brightdef.thresh_type,\
                        brightdef.coherance, instance, matchdef, templatedef)
plotting.threeD_gridplot(detect_nodes)
# detect_templates, detect_nodes=bright_lights.brightness(stations, \
                        # cutnodes, cutlags, stream,
                        # brightdef.threshold, brightdef.thresh_type,\
                        # brightdef.coherance, instance, matchdef, templatedef)
