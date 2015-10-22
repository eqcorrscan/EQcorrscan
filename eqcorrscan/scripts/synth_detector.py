#!/usr/bin/env python
"""
Script to utilise the match_filter routines of EQcorrscan and the synth_seis
routines in the utils of EQcorrscan.  This script will read in a grid of
travel-times and generate a series of synthetic templates for this grid.
It will then use this selection of synthetic templates in a match-filter
routine to detect similar, real-earthquakes in coninuous seismic data.

These detected events can then be stacked and used as further templates...
ad infinitium...
"""
import os, sys, glob, datetime as dt, numpy as np
sys.path.append('/projects/nesi00219/EQcorrscan')
instance=0
Split=False
startdate=False
parallel=True
if len(sys.argv) == 2:
    flag=str(sys.argv[1])
    if flag == '--debug':
        Test=True
        Prep=False
    elif flag == '--debug-prep':
        Test=False
        Prep=True
    else:
        raise IOError("I don't recognise your arguments I know --debug and "+\
        "--debug-prep, and --instance, --splits, --startdate, --enddate")
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
        raise IOError("I don't recognise your arguments,  I know --debug and --debug-prep, and --instance, --splits, --startdate, --enddate")
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

######### END of arguments section #######

from par import template_gen_par as templatedef
from par import match_filter_par as matchdef
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
from utils import pre_processing, synth_seis
from utils import EQcorrscan_plotting as plotting
from obspy.signal.filter import bandpass
from joblib import Parallel, delayed
import warnings

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
#
#
## Use the brightness function to search for possible templates
## First read in the travel times
#print 'Reading in the original grids'
#stations, allnodes, alltravel_times = \
#            bright_lights._read_tt(brightdef.nllpath,brightdef.stations,\
#                                    brightdef.phase, phaseout='S', \
#                                    ps_ratio=brightdef.ps_ratio, lags_switch=False)
#print 'I have read in '+str(len(allnodes))+' nodes'
#
## We now have a grid of travel-times which can then be used to generate synthetic\
## seismograms using the utils.synth_seis functions.
#
## We should trim the grid to the area we want to work in
#print 'Cutting the grid'
#stations, nodes, travel_times = bright_lights._resample_grid(stations, allnodes,
#                                                     alltravel_times,
#                                                     brightdef.mindepth,
#                                                     brightdef.maxdepth,
#                                                     brightdef.corners,
#                                                     brightdef.resolution)
#del allnodes, alltravel_times
## Check that we still have a grid!
#if len(nodes) == 0:
#    raise IOError("You have no nodes left")
#
## Call the template generation function
#synth_templates=synth_seis.template_grid(stations, nodes, travel_times, 'S', \
#        PS_ratio=brightdef.ps_ratio, samp_rate=templatedef.samp_rate,\
#        flength=int(templatedef.samp_rate*templatedef.length))
#print 'We have '+str(len(synth_templates))+' synthetic templates'
#
## Write out the synthetics!
#i=0
#template_names=[] # List of the template names, which will be the node location
#templates=[]
#for synth in synth_templates:
#    # We need the data to be in int32
#    stations=[tr.stats.station for tr in synth if tr.stats.station not in ['WHAT','POCR']]
#    if len(list(set(stations))) < 5:
#        # Only write and use templates with at least five stations
#        i+=1
#        continue
#    for tr in synth:
#        tr.data=(tr.data*1000).astype(np.int32)
#        tr.filter('bandpass', freqmin=templatedef.lowcut,\
#                    freqmax=templatedef.highcut)
#        # Name the channels so that they can be found!
#        if tr.stats.station in ['FOZ','JCZ','LBZ','WVZ']:
#            tr.stats.channel='HHN'
#            tr.stats.network='NZ'
#        elif tr.stats.station == 'GCSZ':
#            tr.stats.channel='EH1'
#            tr.stats.network='NZ'
#        elif tr.stats.station == 'RPZ':
#            tr.stats.channel='HH1'
#            tr.stats.network='NZ'
#        elif tr.stats.station in ['EORO','WHYM','COSA','GOVA','LABE','MTFO',\
#                                    'COVA','POCR','SOLU','WHAT',\
#                                    'MTBA','SOLU']:
#            tr.stats.channel='SHN'
#            tr.stats.network='AF'
#        elif tr.stats.station in ['FRAN','POCR2','WHAT2']:
#            tr.stats.channel='SH2'
#            tr.stats.network='AF'
#    synth.write('templates/synthetics/'+str(nodes[i][0])+'_'+str(nodes[i][1])+\
#                '_'+str(nodes[i][2])+'_template.ms', format='MSEED')#,\
#                #encoding='STEIM2', reclen=512)
#    template_names.append(str(nodes[i][0])+'_'+str(nodes[i][1])+\
#                '_'+str(nodes[i][2]))
#    templates.append(synth)
#    i+=1

#del nodes, travel_time

template_names=glob.glob('templates/synthetics/*_template.ms')
templates=[obsread(tfile) for tfile in template_names]
template_names=[t.split('/')[-1].split('_template.ms')[0] \
                for t in template_names]

print 'We have '+str(len(templates))+' templates with at least five stations'
print 'Working out what stations we have'

stations=[]
for template in templates:
    # Calculate the delays for each template, do this only once so that we
    # don't have to do it heaps!
    # Check that all templates are the correct length
    for tr in template:
        if not templatedef.samp_rate*templatedef.length == tr.stats.npts:
            raise ValueError('Template for '+tr.stats.station+'.'+\
                             tr.stats.channel+' is not the correct length, recut.'+\
                             ' It is: '+str(tr.stats.npts)+' and should be '+
                             str(templatedef.samp_rate*templatedef.length))
    # Generate list of stations in templates
    for tr in template:
        # Correct FOZ channels
        if tr.stats.station=='FOZ' and len(tr.stats.channel)==3:
            tr.stats.channel='HH'+tr.stats.channel[2]
        if len(tr.stats.channel)==3:
            stations.append(tr.stats.station+'.'+tr.stats.channel[0]+\
                            '*'+tr.stats.channel[2]+'.'+tr.stats.network)
            tr.stats.channel=tr.stats.channel[0]+tr.stats.channel[2]
        elif len(tr.stats.channel)==2:
            stations.append(tr.stats.station+'.'+tr.stats.channel[0]+\
                            '*'+tr.stats.channel[1]+'.'+tr.stats.network)
        else:
            raise ValueError('Channels are not named with either three or two characters')
stations=list(set(stations))
# Use the templates to scan through the datas!
# Work out what days are to be scanned through
if startdate:
    dates=[UTCDateTime(startdate)+i for i in xrange(0, int(UTCDateTime(enddate) -\
                                                           UTCDateTime(startdate)),\
                                                    86400)]
else:
    dates=matchdef.dates

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


f=open('detections/synth_'+str(min(dates).year)+\
       str(min(dates).month).zfill(2)+\
       str(min(dates).day).zfill(2)+'_'+\
       str(max(dates).year)+str(max(dates).month).zfill(2)+\
       str(max(dates).day).zfill(2),'w')
f.write('template, detect-time, cccsum, threshold, number of channels\n')

all_templates=templates
all_template_names=template_names

for day in dates:
    if 'st' in locals():
        del st
    # Read in data using obspy's reading routines, data format will be worked
    # out by the obspy module
    # Note you might have to change this bit to match your naming structure
    actual_stations=[] # List of the actual stations used
    for stachan in stations:
        # station is of the form STA.CHAN, to allow these to be in an odd
        # arrangements we can seperate them
        station=stachan.split('.')[0]
        channel=stachan.split('.')[1]
        netcode=stachan.split('.')[2]
        rawdir='/projects/nesi00219/logfiles/Volumes/Taranaki_01/data/boseca/SAMBA_mar09/'+station+'/'+\
                    str(day.year)+str(day.julday).zfill(3)
        #rawdir='/Volumes/Taranaki_01/data/boseca/SAMBA_mar09/'+station+'/'+\
        #            str(day.year)+str(day.julday).zfill(3)
        errors, full = seismo_logs.check_all_logs(rawdir, \
                                                  1.0/templatedef.samp_rate)
        if len(errors) > 1:
            continue
        if not Test:
            # Set up the base directory format
            for base in matchdef.contbase:
                if base[2]==netcode:
                    contbase=base
            if not 'contbase' in locals():
                raise NameError('contbase is not defined for '+netcode)
            baseformat=contbase[1]
            if baseformat=='yyyy/mm/dd':
                daydir=str(day.year)+'/'+str(day.month).zfill(2)+'/'+\
                        str(day.day).zfill(2)
            elif baseformat=='Yyyyy/Rjjj.01':
                daydir='Y'+str(day.year)+'/R'+str(day.julday).zfill(3)+'.01'
            elif baseformat=='yyyymmdd':
                daydir=str(day.year)+str(day.month).zfill(2)+str(day.day).zfill(2)

            # Try and find the appropriate files
            if baseformat=='Yyyyy/Rjjj.01':
                if glob.glob(contbase[0]+'/'+daydir+'/'+station+'.*.'+channel+\
                             '.'+str(day.year)+'.'+str(day.julday).zfill(3)):
                    chan_available=True
                else:
                    chan_available=False
            else:
                if glob.glob(contbase[0]+'/'+daydir+'/*'+station+'.'+channel+'.*'):
                    chan_available=True
                else:
                    chan_available=False
            if chan_available:
                if not 'st' in locals():
                    if baseformat=='Yyyyy/Rjjj.01':
                        st=obsread(contbase[0]+'/'+daydir+'/*'+station+'.*.'+\
                                   channel+'.'+str(day.year)+'.'+\
                                   str(day.julday).zfill(3))
                    else:
                        st=obsread(contbase[0]+'/'+daydir+'/*'+station+'.'+\
                                   channel+'.*')
                else:
                    if baseformat=='Yyyyy/Rjjj.01':
                        st+=obsread(contbase[0]+'/'+daydir+'/*'+station+'.*.'+\
                                    channel+'.'+str(day.year)+'.'+\
                                    str(day.julday).zfill(3))
                    else:
                        st+=obsread(contbase[0]+'/'+daydir+'/*'+station+'.'+\
                                    channel+'.*')
                actual_stations.append(station) # Add to this list only if we have the data
            else:
                print 'No data for '+stachan+' for day '+daydir+' in '\
                        +contbase[0]
        else:
            fname='test_data/'+station+'-'+channel+'-'+str(day.year)+\
                           '-'+str(day.month).zfill(2)+\
                           '-'+str(day.day).zfill(2)+'-processed.ms'
            if glob.glob(fname):
                if not 'st' in locals():
                    st=obsread(fname)
                else:
                    st+=obsread(fname)
                actual_stations.append(station)
    actual_stations=list(set(actual_stations))

    st=st.merge(fill_value='interpolate') # Enforce trace continuity
    if not 'st' in locals():
        print 'No data found for day: '+str(day)
    elif len(actual_stations) < matchdef.minsta:
        print 'Data from fewer than '+str(matchdef.minsta)+' stations found, will not detect'
    else:
        if not Test:
            # Process data
            print 'Processing the data for day '+daydir
            if matchdef.debug >= 4:
                for tr in st:
                    tr=pre_processing.dayproc(tr, templatedef.lowcut, templatedef.highcut,\
                                            templatedef.filter_order, templatedef.samp_rate,\
                                            matchdef.debug, day)
            else:
                st=Parallel(n_jobs=10)(delayed(pre_processing.dayproc)(tr, templatedef.lowcut,\
                                                                   templatedef.highcut,\
                                                                   templatedef.filter_order,\
                                                                   templatedef.samp_rate,\
                                                                   matchdef.debug, day)\
                                for tr in st)
        if not Prep:
            # For some reason st is now a list rather than a stream
            if 'stream_st' in locals():
                del stream_st
            for tr in st:
                if 'stream_st' in locals():
                    stream_st+=tr
                else:
                    stream_st=Stream(tr)
            st=stream_st
            # Call the match_filter module - returns detections, a list of detections
            # containted within the detection class with elements, time, template,
            # number of channels used and cross-channel correlation sum.
            print 'Running the detection routine'
            if not os.path.isdir('temp_'+str(instance)):
                os.makedirs('temp_'+str(instance))
            groups=0
            detections=[]
            # Cope with having heaps of templates
            if len(all_templates) > 100:
                groups=int(len(all_templates)/100)
            for i in xrange(groups):
                if i==groups:
                    templates=all_templates[i*100:]
                    template_names=all_template_names[i*100:]
                else:
                    templates=all_templates[i*100:(i+1)*100]
                    template_names=all_template_names[i*100:(i+1)*100]
                detections+=match_filter.match_filter(template_names, templates, st,
                                                 matchdef.threshold, matchdef.threshtype,
                                                 matchdef.trig_int,  matchdef.plot,
                                                 'temp_'+str(instance))

            for detection in detections:
                # output detections to file
                f.write(detection.template_name+', '+str(detection.detect_time)+\
                        ', '+str(detection.detect_val)+', '+str(detection.threshold)+\
                        ', '+str(detection.no_chans)+'\n')
                print 'template: '+detection.template_name+' detection at: '\
                    +str(detection.detect_time)+' with a cccsum of: '+str(detection.detect_val)
            if detections:
                f.write('\n')
        else:
            for tr in st:
                tr.write('test_data/'+tr.stats.station+'-'+tr.stats.channel+\
                         '-'+str(tr.stats.starttime.year)+\
                         '-'+str(tr.stats.starttime.month).zfill(2)+\
                         '-'+str(tr.stats.starttime.day).zfill(2)+\
                         '-processed.ms', format='MSEED')
f.close()
