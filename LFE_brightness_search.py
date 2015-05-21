#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Script to call all elements of EQcorrscan module to search
#               continuous data for likely LFE repeats
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
LFE_brightness_search - Script to generate templates using the birghtness function
then seaach for repeats of them in contnuous data.
"""

import sys, glob
sys.path.insert(0,"/home/processor/Desktop/EQcorrscan")

if len(sys.argv) == 2:
    flag=str(sys.argv[1])
    if flag == '--debug':
        Test=True
        Prep=False
    elif flag == '--debug-prep':
        Test=False
        Prep=True
    else:
        raise ValueError("I don't recognise the argument, I only know --debug and --debug-prep")
elif not len(sys.argv) == 1:
    raise ValueError("I only take one argument, no arguments, or two flags with arguments")
else:
    Test=False
    Prep=False
    Split=False

from par import template_gen_par as templatedef
from par import match_filter_par as matchdef
from par import bright_lights_par as brightdef
#from par import lagcalc as lagdef
from obspy import UTCDateTime, read as obsread
# First generate the templates
from core import bright_lights, match_filter
from utils import pre_processing
from obspy.signal.filter import bandpass

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


# Use the brightness function to search for possible templates
# First read in the travel times
print 'Reading in the original grids'
stations, allnodes, alllags = \
            bright_lights._read_tt(brightdef.nllpath,brightdef.stations,\
                                    brightdef.phase)
# Resample the grid to allow us to run it quickly!
print 'Cutting the grid'
stations, nodes, lags = bright_lights._resample_grid(stations, allnodes,
                                                     alllags,
                                                     brightdef.volume,
                                                     brightdef.resolution)
# Remove lags that have a similar network moveout, e.g. the sum of the
# differences in moveouts is small.
print "Removing simlar lags"
stations, nodes, lags = bright_lights._rm_similarlags(stations, nodes, lags,
                                                      brightdef.nodesimthresh)
print "Plotting new grid"
lats=[]
longs=[]
depths=[]
for node in nodes:
    lats.append(node[0])
    longs.append(node[1])
    depths.append(node[2])
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(lats, longs, depths)
plt.show()
# Call the main function!
templates=[]
for i in xrange(0,int(brightdef.enddate-brightdef.startdate),86400): #Loop through dates
    # Set up where data are going to be read in from
    day=brightdef.startdate+i
    print 'Working on day '+str(day)
    if not Test:
                # Read in the data
        if 'stream' in locals():
            del stream
        for station in stations:
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
            print '     Reading data from: '+contbase[0]+'/'+daydir+'/*'+station+'*'
            for chan in ['N','E','1','2']: # only take horizontal components
                if glob.glob(contbase[0]+'/'+daydir+'/*'+station+'*'+chan+'*'):
                    if not 'stream' in locals():
                        stream=obsread(contbase[0]+'/'+daydir+'/*'+station+'*'+chan+'*')
                    else:
                        stream+=obsread(contbase[0]+'/'+daydir+'/*'+station+'*'+chan+'*')
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
    print 'Processing the data'
    if not Test:
        for tr in stream:
            tr=pre_processing.dayproc(tr, templatedef.lowcut, templatedef.highcut,\
                                        templatedef.filter_order, templatedef.samp_rate,\
                                        matchdef.debug, day)
    if not Prep:
        templates+=bright_lights.brightness(stations, nodes, lags, stream,\
                        brightdef.threshold, brightdef.thresh_type)
    else:
        for tr in stream:
            tr.write('test_data/'+tr.stats.station+'-'+tr.stats.channel+\
                    '-'+str(tr.stats.starttime.year)+\
                    '-'+str(tr.stats.starttime.month).zfill(2)+\
                    '-'+str(tr.stats.starttime.day).zfill(2)+\
                    '-processed.ms', format='MSEED')
        sys.exit()

import numpy as np
# Now do the detections with these templates!
station=[]
print np.shape(templates)
for template in templates:
    # Calculate the delays for each template, do this only once so that we
    # don't have to do it heaps!
    # Get minimum start time
    mintime=UTCDateTime(3000,1,1,0,0)
    for tr in template:
        if tr.stats.starttime < mintime:
            mintime=tr.stats.starttime
    delay=[]
    # Generate list of delays
    for tr in template:
        delay.append(tr.stats.starttime-mintime)
    delays.append(delay)
    # Generate list of stations in templates
    for tr in template:
        station.append(tr.stats.station+'.'+tr.stats.channel[0]+\
                        '*'+tr.stats.channel[1])
        stations.append(tr.stats.station+'.'+tr.stats.channel[0]+\
                        '*'+tr.stats.channel[1])
    # Save a miniseed copy of the templates
    #print 'Writing template as: '+str(template[0].stats.starttime)+'.ms'
    #template.write(templatedef.saveloc+'/'+str(template[0].stats.starttime)+'.ms',format="MSEED")

# Sort stations into a unique list - this list will ensure we only read in data
# we need, which is VERY important as I/O is very costly and will eat memory
stations=list(set(stations))
print len(stations)
for station in stations:
    print station

# Now run the match filter routine


# Loop over days
ndays=int(matchdef.enddate-matchdef.startdate+1)
newsfiles=[]
for i in range(0,ndays):
    # Set up where data are going to be read in from
    day=matchdef.startdate+(i*86400)
    if matchdef.baseformat=='yyyy/mm/dd':
        daydir=str(day.year)+'/'+str(day.month).zfill(2)+'/'+\
                str(day.day).zfill(2)
    elif matchdef.baseformat=='Yyyyy/Rjjj.01':
        daydir='Y'+str(day.year)+'/R'+str(day.julday).zfill(3)+'.01'
    # Read in data using obspy's reading routines, data format will be worked
    # out by the obspy module
    # Note you might have to change this bit to match your naming structure
    for stachan in stations:
        # station is of the form STA.CHAN, to allow these to be in an odd
        # arrangements we can seperate them
        try:
            station=stachan.split('.')[0]
            channel=stachan.split('.')[1]
        except:
            print 'Issues with this station name: '+stachan
            sys.exit()
        if glob.glob(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*'):
            #print 'Reading data from: '+\
            #    glob.glob(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*')[0]
            if not 'st' in locals():
                st=obsread(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*')
            else:
                st+=obsread(matchdef.contbase+'/'+daydir+'/'+station+'*.'+channel+'.*')
        else:
            print 'No data for '+stachan+' for day '+daydir+' in '+matchdef.contbase
    if not 'st' in locals():
        print 'No data found for day: '+str(day)
        break
    # Process data
    print 'Processing the data'
    for tr in st:
        if tr.stats.sampling_rate != templatedef.samp_rate:
            #print 'Reasmpling '+tr.stats.station+'.'+tr.stats.channel
            tr.resample(templatedef.samp_rate)
        tr.detrend('linear')    # Detrend data before filtering
        #print 'Filtering '+tr.stats.station+'.'+tr.stats.channel
        tr.data=bandpass(tr.data, templatedef.lowcut, templatedef.highcut,
                    tr.stats.sampling_rate, 4, True)
        # Account for two letter channel names in s-files and therefore templates
        tr.stats.channel=tr.stats.channel[0]+tr.stats.channel[2]
        # Sanity check to ensure files are daylong
        if int(tr.stats.npts/tr.stats.sampling_rate) != 86400:
            print 'Data for '+tr.stats.station+'.'+tr.stats.channel+' is not of daylong length, will zero pad'
            # Use obspy's trim function with zero padding
            tr.trim(day,day+86400,pad=True,fill_value=0)


    # Call the match_filter module - returns detections, a list of detections
    # containted within the detection class with elements, time, template,
    # number of channels used and cross-channel correlation sum.
    print 'Running the detection routine'
    detections=match_filter.match_filter(templatedef.sfiles, templates, delays, st,
                                         matchdef.threshold, matchdef.threshtype,
                                         matchdef.trig_int*st[0].stats.sampling_rate,
                                         matchdef.plot)

    for detection in detections:
        print 'template: '+detection.template_name+' detection at: '\
            +str(detection.detect_time)+' with a cccsum of: '+str(detection.detect_val)
    # Call the lag generation routine - returns list of s-files generated
    #newsfiles+=lagcalc(detections, st)
