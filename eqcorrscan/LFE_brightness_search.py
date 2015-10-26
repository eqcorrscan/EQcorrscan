#!/usr/bin/env python

#------------------------------------------------------------------------------
#   Purpose:    Script to call all elements of EQcorrscan module to search
#               continuous data for likely LFE repeats
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
LFE_brightness_search - Script to generate templates using the birghtness function
then seaach for repeats of them in contnuous data.

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

import sys, glob
import datetime as dt
instance=0
Split=False
startdate=False
sys.path.insert(0,"/home/processor/Desktop/EQcorrscan")
parallel=True
oldnodes=False
if len(sys.argv) == 2:
    flag=str(sys.argv[1])
    if flag == '--debug':
        Test=True
        Prep=False
    elif flag == '--debug-prep':
        Test=False
        Prep=True
    elif flag == '--old-nodes':
        Test=False
        Prep=False
        oldnodes=True
    else:
        raise IOError("I don't recognise your arguments I know --debug and --debug-prep, and --instance, --splits, --startdate, --enddate")
elif len(sys.argv) == 5:
    args=sys.argv[1:len(sys.argv)]
    if args[0] == '--instance' or args[2]=='--instance':
        # Arguments to allow the code to be run in multiple instances
        Split=True
        Test=False
        Prep=False
        for i in xrange(len(args)):
            if args[i] == '--instance':
                instance=int(args[i+1])
                print 'I will run this for instance '+str(instance)
            elif args[i] == '--splits':
                splits=int(args[i+1])
                print 'I will divide the days into '+str(splits)+' chunks'
    elif args[0] == '--startdate' or args[2] == '--startdate':
        Split=False
        Test=False
        Prep=False
        for i in xrange(len(args)):
            if args[i]== '--startdate':
                startdate=dt.datetime.strptime(str(args[i+1]),'%Y/%m/%d')
                print 'Start date is: '+dt.datetime.strftime(startdate, '%Y/%m/%d')
            elif args[i] == '--enddate':
                enddate=dt.datetime.strptime(str(args[i+1]), '%Y/%m/%d')
                print 'End date is: '+dt.datetime.strftime(enddate, '%Y/%m/%d')
    else:
        raise IOError("I don't recognise your arguments,  I know --debug and"+\
                      "--debug-prep, and --instance, --splits, --startdate, --enddate")
elif len(sys.argv) == 6:
    args=sys.argv[1:len(sys.argv)]
    if args[0] == '--instance' or args[1]=='--instance' or args[2]=='--instance':
        # Arguments to allow the code to be run in multiple instances
        Split=True
        Test=False
        Prep=False
        for i in xrange(len(args)):
            if args[i] == '--instance':
                instance=int(args[i+1])
                print 'I will run this for instance '+str(instance)
            elif args[i] == '--splits':
                splits=int(args[i+1])
                print 'I will divide the days into '+str(splits)+' chunks'
            elif args[i] == '--old-nodes':
                oldnodes=True
    elif args[0] == '--startdate' or argc[1] == '--startdate' or args[2] == '--startdate':
        Split=False
        Test=False
        Prep=False
        for i in xrange(len(args)):
            if args[i]== '--startdate':
                startdate=dt.datetime.strptime(str(args[i+1]),'%Y/%m/%d')
                print 'Start date is: '+dt.datetime.strftime(startdate, '%Y/%m/%d')
            elif args[i] == '--enddate':
                enddate=dt.datetime.strptime(str(args[i+1]), '%Y/%m/%d')
                print 'End date is: '+dt.datetime.strftime(enddate, '%Y/%m/%d')
            elif args[i] == '--old-nodes':
                oldnodes=True
    else:
        raise IOError("I don't recognise your arguments,  I know --debug and"+\
                      "--debug-prep, and --instance, --splits, --startdate, --enddate")

elif len(sys.argv) == 7:
    args=sys.argv[1:len(sys.argv)]
    Split=False
    Test=False
    Prep=False
    for i in xrange(len(args)):
        if args[i]== '--startdate':
            startdate=dt.datetime.strptime(str(args[i+1]),'%Y/%m/%d')
            print 'Start date is: '+dt.datetime.strftime(startdate, '%Y/%m/%d')
        elif args[i] == '--enddate':
            enddate=dt.datetime.strptime(str(args[i+1]), '%Y/%m/%d')
            print 'End date is: '+dt.datetime.strftime(enddate, '%Y/%m/%d')
        elif args[i] == '--instance':
            instance=int(args[i+1])
            print 'This instance is given the flag '+str(instance)
elif not len(sys.argv) == 1:
    raise ValueError("I only take one argument, no arguments, or two flags with arguments")
else:
    Test=False
    Prep=False
    Split=False

from par import LFE_template_gen_par as templatedef
from par import match_filter_par_LFEs as matchdef
from par import bright_lights_par as brightdef
from utils import seismo_logs
if brightdef.plotsave:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.ioff()
#from par import lagcalc as lagdef
from obspy import UTCDateTime, Stream, read as obsread
# First generate the templates
from core import bright_lights, match_filter
from utils import pre_processing
from utils import EQcorrscan_plotting as plotting
from obspy.signal.filter import bandpass
from joblib import Parallel, delayed
import warnings, pickle

templates=[]
delays=[]
stations=[]
print 'Template generation parameters are:'
print '     sfilebase: '+templatedef.sfilebase
print '     samp_rate: '+str(templatedef.samp_rate)+' Hz'
print '     lowcut: '+str(templatedef.lowcut)+' Hz'
print '     highcut: '+str(templatedef.highcut)+' Hz'
print '     length: '+str(templatedef.length)+' s'
print '     swin: '+templatedef.swin+'\n'

if not oldnodes:
    # Use the brightness function to search for possible templates
    # First read in the travel times
    print 'Reading in the original grids'
    stations, allnodes, alllags = \
            bright_lights._read_tt(brightdef.nllpath,brightdef.stations,\
                                    brightdef.phase, phaseout='S', \
                                    ps_ratio=brightdef.ps_ratio)
    print 'I have read in '+str(len(allnodes))+' nodes'
    # Resample the grid to allow us to run it quickly!
    print 'Cutting the grid'
    stations, nodes, lags = bright_lights._resample_grid(stations, allnodes,
                                                         alllags,
                                                         brightdef.mindepth,
                                                         brightdef.maxdepth,
                                                         brightdef.corners,
                                                         brightdef.resolution)
    del allnodes, alllags
    # Check that we still have a grid!
    if len(nodes) == 0:
        raise IOError("You have no nodes left")
    # Remove lags that have a similar network moveout, e.g. the sum of the
    # differences in moveouts is small.
    print "Removing simlar lags"
    stations, nodes, lags = bright_lights._rm_similarlags(stations, nodes, lags,
                                                          brightdef.nodesimthresh)
   # print "Plotting new grid"
   # plotting.threeD_gridplot(nodes, save=brightdef.plotsave, savefile='Nodes_in.png')
    # Call the main function!

    fnodesin=open('Nodes_in.csv','w')
    for node in nodes:
        fnodesin.write(node[0]+','+node[1]+','+node[2]+'\n')
    fnodesin.close()
    with open('Nodes_in.pkl','wb') as pickle_file:
        pickle.dump(nodes, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
    with open('Lags_in.pkl','wb') as pickle_file:
        pickle.dump(lags, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
    with open('Stations_in.pkl','wb') as pickle_file:
        pickle.dump(stations, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
else:
    with open('Nodes_in.pkl','rb') as pickle_load:
        nodes = pickle.load(pickle_load)
    with open('Lags_in.pkl','rb') as pickle_load:
        lags = pickle.load(pickle_load)
    with open('Stations_in.pkl','rb') as pickle_load:
        stations = pickle.load(pickle_load)
    print 'Read in '+str(len(nodes))+' nodes'

templates=[]
nodesout=[]

if startdate:
    dates=[UTCDateTime(startdate)+i for i in xrange(0, int(UTCDateTime(enddate) -\
                                                           UTCDateTime(startdate)),\
                                                    86400)]
else:
    dates=brightdef.dates

ndays=len(dates)
print 'Will loop through '+str(ndays)+' days'
if Split:
    if instance==splits-1:
        ndays=ndays-(ndays/splits)*(splits-1)
        dates=dates[-ndays:]
    else:
        ndays=ndays/splits
        dates=dates[ndays*instance:(ndays*instance)+ndays]
    print 'This instance will run for '+str(ndays)+' days'
    print 'This instance will run from '+str(min(dates))


for day in dates: #Loop through dates
    # Set up where data are going to be read in from
    if 'stream' in locals():
        del stream
    print 'Working on day '+str(day)
    if not Test:
                # Read in the data
        # if 'stream' in locals():
            # del stream
        for station in stations:
            # Check logs for this day
            if station == 'WHAT2':
                sta='WHAT'
                useful_chans=['2']
            elif station == 'POCR2':
                sta='POCR'
                useful_chans=['2']
            elif station == 'FRAN':
                sta=station
                useful_chans=['2']
            else:
                sta=station
                useful_chans=['N','2']
            rawdir='/projects/nesi00219/logfiles/Volumes/Taranaki_01/data/boseca/SAMBA_mar09/'+sta+'/'+\
                    str(day.year)+str(day.julday).zfill(3)
            errors, full = seismo_logs.check_all_logs(rawdir, \
                                                      1.0/templatedef.samp_rate)
            if len(errors) > 1:
                warnings.warn('Found GPS errors exceeding one sample for station '+station)
                continue
            else:
                print 'No GPS errors found, continuing'
            if len(station) == 3:
                netcode='NZ'
            else:
                netcode='AF'
            for base in matchdef.contbase:
                if base[2]==netcode:
                    contbase=base
            if not 'contbase' in locals():
                raise NameError('contbase not defined for netcode '+netcode)
            if contbase[1]=='yyyymmdd':
                daydir=str(day.year)+str(day.month).zfill(2)+\
                        str(day.day).zfill(2)
            elif contbase[1]=='Yyyyy/Rjjj.01':
                daydir='Y'+str(day.year)+'/R'+str(day.julday).zfill(3)+'.01'
            print '     Reading data from: '
            for chan in useful_chans: # only take N horizontal components
                if glob.glob(contbase[0]+'/'+daydir+'/*'+station+'.*'+chan+'.*'):
                    print contbase[0]+'/'+daydir+'/*'+station+'.*'+chan+'.*'
                    if not 'stream' in locals():
                        stream=obsread(contbase[0]+'/'+daydir+'/*'+station+'.*'+chan+'.*')
                    else:
                        stream+=obsread(contbase[0]+'/'+daydir+'/*'+station+'.*'+chan+'.*')
    else:
        for station in stations:
            fname='test_data/'+station+'-*-'+str(day.year)+\
                               '-'+str(day.month).zfill(2)+\
                               '-'+str(day.day).zfill(2)+'-processed.ms'
            if glob.glob(fname):
                if not 'stream' in locals():
                    stream=obsread(fname)
                else:
                    stream+=obsread(fname)
    # Process the stream
    if not Test:
        print 'Processing the data'
        stream=stream.merge(fill_value='interpolate')
                        # Merge stream so that each trace is a single channel to
                       # send to pre-processing
        if not parallel:
            for tr in stream:
                # tr.plot()
                tr=pre_processing.dayproc(tr, brightedef.lowcut, brightdef.highcut,\
                                            brightdef.filter_order, brightdef.samp_rate,\
                                            templatedef.debug, day)
        else:
            stream=Parallel(n_jobs=10)(delayed(pre_processing.dayproc)(tr, brightdef.lowcut,\
                                                                   brightdef.highcut,\
                                                                   brightdef.filter_order,\
                                                                   brightdef.samp_rate,\
                                                                   templatedef.debug, day)\
                                for tr in stream)
            stream=Stream(stream)
            print stream
    if not Prep:
        #stream_copy=stream.copy() # Keep the stream safe
        print "Running the detection routine"
        # Check that the data are okay
        detect_templates, detect_nodes=bright_lights.brightness(stations, \
                        nodes, lags, stream,
                        brightdef.threshold, brightdef.thresh_type,\
                        brightdef.coherance, instance, matchdef, templatedef)
        del detect_templates#, stream # Delete templates from memory to conserve RAM!
        #stream=stream_copy
        nodesout+=detect_nodes
        if Split:
            plotting.threeD_gridplot(nodesout, save=brightdef.plotsave,\
                                 savefile='Detected_nodes_'+str(instance)+'.png')
        else:
            plotting.threeD_gridplot(nodesout, save=brightdef.plotsave,\
                                 savefile='Detected_nodes.png')

    else:
        for tr in stream:
            print "Writing data as: test_data/"+tr.stats.station+'-'+tr.stats.channel+\
                    '-'+str(tr.stats.starttime.year)+\
                    '-'+str(tr.stats.starttime.month).zfill(2)+\
                    '-'+str(tr.stats.starttime.day).zfill(2)+\
                    '-processed.ms'
            tr.write('test_data/'+tr.stats.station+'-'+tr.stats.channel+\
                    '-'+str(tr.stats.starttime.year)+\
                    '-'+str(tr.stats.starttime.month).zfill(2)+\
                    '-'+str(tr.stats.starttime.day).zfill(2)+\
                    '-processed.ms', format='MSEED')
        sys.exit()

import numpy as np
# Now do the detections with these templates!
# station=[]
# print np.shape(templates)
# for template in templates:
    # # Calculate the delays for each template, do this only once so that we
    # # don't have to do it heaps!
    # # Get minimum start time
    # mintime=UTCDateTime(3000,1,1,0,0)
    # for tr in template:
        # if tr.stats.starttime < mintime:
            # mintime=tr.stats.starttime
    # delay=[]
    # # Generate list of delays
    # for tr in template:
        # delay.append(tr.stats.starttime-mintime)
    # delays.append(delay)
    # # Generate list of stations in templates
    # for tr in template:
        # station.append(tr.stats.station+'.'+tr.stats.channel[0]+\
                        # '*'+tr.stats.channel[1])
        # stations.append(tr.stats.station+'.'+tr.stats.channel[0]+\
                        # '*'+tr.stats.channel[1])
    # Save a miniseed copy of the templates
    #print 'Writing template as: '+str(template[0].stats.starttime)+'.ms'
    #template.write(templatedef.saveloc+'/'+str(template[0].stats.starttime)+'.ms',format="MSEED")

# Sort stations into a unique list - this list will ensure we only read in data
# we need, which is VERY important as I/O is very costly and will eat memory
stations=list(set(stations))

if Split:
    node_file=open('Nodes_detected_'+str(instance)+'.csv','w')
else:
    node_file=open('Nodes_detected.csv','w')
for node in nodesout:
    node_file.write(node[0]+','+node[1]+','+node[2]+'\n')
node_file.close()

# print len(stations)
# for station in stations:
    # print station

# Now run the match filter routine

# # Loop over days
# ndays=int(matchdef.enddate-matchdef.startdate+1)
# newsfiles=[]
# for i in range(0,ndays):
    # # Set up where data are going to be read in from
    # day=matchdef.startdate+(i*86400)
    # if matchdef.baseformat=='yyyy/mm/dd':
        # daydir=str(day.year)+'/'+str(day.month).zfill(2)+'/'+\
                # str(day.day).zfill(2)
    # elif matchdef.baseformat=='Yyyyy/Rjjj.01':
        # daydir='Y'+str(day.year)+'/R'+str(day.julday).zfill(3)+'.01'
    # # Read in data using obspy's reading routines, data format will be worked
    # # out by the obspy module
    # # Note you might have to change this bit to match your naming structure
    # for stachan in stations:
        # # station is of the form STA.CHAN, to allow these to be in an odd
        # # arrangements we can seperate them
        # try:
            # station=stachan.split('.')[0]
            # channel=stachan.split('.')[1]
        # except:
            # print 'Issues with this station name: '+stachan
            # sys.exit()
        # if glob.glob(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*'):
            # #print 'Reading data from: '+\
            # #    glob.glob(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*')[0]
            # if not 'st' in locals():
                # st=obsread(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*')
            # else:
                # st+=obsread(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*')
        # else:
            # print 'No data for '+stachan+' for day '+daydir+' in '+matchdef.contbase
    # if not 'st' in locals():
        # print 'No data found for day: '+str(day)
        # break
    # # Process data
    # print 'Processing the data'
    # for tr in st:
        # tr=pre_processing.dayproc(tr, templatedef.lowcut, templatedef.highcut,\
                                            # templatedef.filter_order, templatedef.samp_rate,\
                                            # matchdef.debug, day)

    # # Call the match_filter module - returns detections, a list of detections
    # # containted within the detection class with elements, time, template,
    # # number of channels used and cross-channel correlation sum.
    # print 'Running the detection routine'
    # detections=match_filter.match_filter(templatedef.sfiles, templates, delays, st,
                                         # matchdef.threshold, matchdef.threshtype,
                                         # matchdef.trig_int*st[0].stats.sampling_rate,
                                         # matchdef.plot)

    # for detection in detections:
        # print 'template: '+detection.template_name+' detection at: '\
            # +str(detection.detect_time)+' with a cccsum of: '+str(detection.detect_val)
    # # Call the lag generation routine - returns list of s-files generated
    # #newsfiles+=lagcalc(detections, st)
