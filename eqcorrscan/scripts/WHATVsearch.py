#!/usr/bin/env python

#------------------------------------------------------------------------------
#   Purpose:    Script to call all elements of EQcorrscan module to search
#               continuous data for likely LFE repeats
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
LFEsearch - Script to generate templates from previously picked EQs and then
search for repeats of them in contnuous data.

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

import sys, os, glob
bob=os.path.realpath(__file__)
bob=bob.split('/')
path='/'
for i in xrange(len(bob)):
    path+=bob[i]+'/'
print path
sys.path.insert(0,path)


from par import template_gen_par as templatedef
from par import match_filter_par as matchdef
#from par import lagcalc as lagdef
from obspy import UTCDateTime, Stream, read as obsread
# First generate the templates
from core import template_gen
from utils import seismo_logs

Split=False
instance=False
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
elif len(sys.argv) == 5:
    # Arguments to allow the code to be run in multiple instances
    Split=True
    Test=False
    Prep=False
    args=sys.argv[1:len(sys.argv)]
    for i in xrange(len(args)):
        if args[i] == '--instance':
            instance=int(args[i+1])
            print 'I will run this for instance '+str(instance)
        elif args[i] == '--splits':
            splits=int(args[i+1])
            print 'I will divide the days into '+str(splits)+' chunks'

elif not len(sys.argv) == 1:
    raise ValueError("I only take one argument, no arguments, or two flags with arguments")
else:
    Test=False
    Prep=False
    Split=False

templates=[]
delays=[]
stations=[]
print 'Template generation parameters are:'
print 'sfilebase: '+templatedef.sfilebase
print 'samp_rate: '+str(templatedef.samp_rate)+' Hz'
print 'lowcut: '+str(templatedef.lowcut)+' Hz'
print 'highcut: '+str(templatedef.highcut)+' Hz'
print 'length: '+str(templatedef.length)+' s'
print 'swin: '+templatedef.swin+'\n'
for sfile in templatedef.sfiles:
    print 'Working on: '+sfile+'\r'
    if not os.path.isfile(templatedef.saveloc+'/'+sfile+'_template.ms'):
        print sfile
        template=template_gen.from_contbase(sfile)

        print 'saving template as: '+templatedef.saveloc+'/'+\
                str(template[0].stats.starttime)+'.ms'
        template.write(templatedef.saveloc+'/'+\
                   sfile+'_template.ms',format="MSEED")
    else:
        template=obsread(templatedef.saveloc+'/'+sfile+'_template.ms')
    templates+=[template]
    # Will read in seisan s-file and generate a template from this,
    # returned name will be the template name, used for parsing to the later
    # functions

# for tfile in templatedef.tfiles:
    # # Loop through pre-existing template files
    # sys.stdout.write("\rReading in pre-existing template: "+tfile+"\r")
    # sys.stdout.flush()
    # templates.append(obsread(tfile))

templates=[obsread(tfile) for tfile in templatedef.tfiles]

print 'Read in '+str(len(templates))+' templates'

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
            raise ValueError('Channels are not named with either three or two charectars')


# Template generation and processing is over, now to the match-filtering

# Sort stations into a unique list - this list will ensure we only read in data
# we need, which is VERY important as I/O is very costly and will eat memory
stations=list(set(stations))

# Now run the match filter routine
from core import match_filter
from obspy import read as obsread
# from obspy.signal.filter import bandpass
# from obspy import Stream, Trace
# import numpy as np
from utils import pre_processing
from joblib import Parallel, delayed

# Loop over days
ndays=len(matchdef.dates)
newsfiles=[]
# f=open('detections/run_start_'+str(UTCDateTime().year)+\
       #str(UTCDateTime().month).zfill(2)+\
       #str(UTCDateTime().day).zfill(2)+'T'+\
       #str(UTCDateTime().hour).zfill(2)+str(UTCDateTime().minute).zfill(2),'w')
print 'Will loop through '+str(ndays)+' days'
if Split:
    if instance==splits-1:
        ndays=ndays-(ndays/splits)*(splits-1)
        dates=matchdef.dates[-ndays:]
    else:
        ndays=ndays/splits
        dates=matchdef.dates[ndays*instance:(ndays*instance)+ndays]
    print 'This instance will run for '+str(ndays)+' days'
    print 'This instance will run from '+str(min(dates))
else:
    dates=matchdef.dates

f=open('detections/'+str(min(dates).year)+\
       str(min(dates).month).zfill(2)+\
       str(min(dates).day).zfill(2)+'_'+\
       str(max(dates).year)+str(max(dates).month).zfill(2)+\
       str(max(dates).day).zfill(2),'w')
f.write('template, detect-time, cccsum, threshold, number of channels\n')

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
        rawdir='/Volumes/Taranaki_01/data/boseca/SAMBA_mar09/'+station+'/'+\
                    str(day.year)+str(day.julday).zfill(3)
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
            template_names=[tfile.split('/')[-1] for tfile in templatedef.tfiles]
            template_names.append(templatedef.sfiles)
            if not os.path.isdir('temp_'+str(instance)):
                os.makedirs('temp_'+str(instance))
            detections=match_filter.match_filter(template_names, templates, st,
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
