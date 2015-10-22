#!/usr/bin/env python
"""
Tutorial

    This script is designed as a tutorial to highlight how to call the main
    functions within the EQcorrscan module.  In this tutorial we will see how
    to generate a template and run this through the match-filter routine.

    The template will be generated from a pre-picked earthquake, however there
    are other ways to generate templates, for example this package also contains
    a simple brightness function that is designed to scan continuous seismic
    data and look for impulsive energy originating from a discrete point in a
    seismic velocity model.  The use of this brightness function is not
    included in this tutorial script yet because it is still in beta.

This package is dstributed under the LGPL v3.0, by using this script and the
functions contained within the EQcorrscan package you implicitly accept the
licence.  For the full wording of the licence refer to the licence.txt file.

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

# First we import the required modules:
from obspy import read, Stream
from core import template_gen, match_filter
from par import match_filter_par as matchdef
from utils import pre_processing, Sfile_util
import glob

# Now we find the s-file we want to use to generate a template from
sfiles=glob.glob('test_data/tutorial_data/*L.S*')

# Generate the template from these sfiles:
templates=[] # Open a list to be filled - only applies for multiple templates
template_names=[] # List of template names for later ID
i=0 # Template name iterator
for sfile in sfiles:
    # Read in the picks from the S-file, note, in the full case one fo the main\
            # functions in template_gen would be used rather than this, but for\
            # the tutorial we will read in the data here - also note that this\
            # template generation is inefficient for multiple templates, if using\
            # daylong data for multiple templates you would want to only read\
            # the seismic data once and cut it multiple times.
    picks=Sfile_util.readpicks(sfile)
    for pick in picks:
        if not 'wavefiles' in locals():
            wavefiles=glob.glob('test_data/tutorial_data/'+\
                                   pick.station+'.*')
        else:
            wavefiles+=glob.glob('test_data/tutorial_data/'+\
                                 pick.station+'.*')
    wavefiles=list(set(wavefiles))
    for wavefile in wavefiles:
        print 'Reading data from '+wavefile
        if not 'st' in locals():
            st=read(wavefile)
        else:
            st+=read(wavefile)
    st=st.merge(fill_value='interpolate')
    day=st[0].stats.starttime.date
    for tr in st:
        tr=pre_processing.dayproc(tr, 1.0, 20.0, 3, 100.0,\
                                  matchdef.debug, day)
    # Apply a small amoutn of delay before the pick
    for pick in picks:
        pick.time=pick.time-0.1
    template=template_gen._template_gen(picks, st, 1.0, 'all')
    # This will generate an obspy.Stream object
    # Append this Stream to the list of templates
    templates+=[template]
    template_names.append('tutorial_'+str(i))
    # Plot the template just to check that all is well!
    template.plot(size=(800,600), equal_scale=False)
    # Save template for later
    template.write('test_data/tutorial_data/'+template_names[i]+'_template.ms',\
                   format='MSEED')
    i+=1
    del template, st

# Extract the stations from the templates
for template in templates:
    if not 'stachans' in locals():
        stachans=[(tr.stats.station, tr.stats.channel) for tr in template]
    else:
        stachans+=[(tr.stats.station, tr.stats.channel) for tr in template]

# Make this a unique list
stachans=list(set(stachans))

# Read in the continuous data for these station, channel combinations
for stachan in stachans:
    print 'Reading data from: test_data/tutorial_data/'+stachan[0]+'.*..*'+stachan[1][-1]+'.*'
    # Generate a new stream object and add to it
    if not 'st' in locals():
        st=read('test_data/tutorial_data/'+stachan[0]+'.*..*'+stachan[1][-1]+'.*')
    else:
        st+=read('test_data/tutorial_data/'+stachan[0]+'.*..*'+stachan[1][-1]+'.*')

# Merge the data to account for miniseed files being written in chunks
# We need continuous day-long data, so data are padded if there are gaps
st=st.merge(fill_value='interpolate')

# Work out what day we are working on, required as we will pad the data to be daylong
day=st[0].stats.starttime.date

# Process the data in the same way as the template
for tr in st:
    tr=pre_processing.dayproc(tr, 1.0, 20.0, 3, 100.0,\
                              matchdef.debug, day)

# Compute detections
detections=match_filter.match_filter(template_names, templates, st,\
                                     matchdef.threshold, matchdef.threshtype,\
                                     matchdef.trig_int, True,\
                                     'temp_0')

# We now have a list of detections! We can output these to a file to check later
f=open('tutorial_detections.csv','w')
for detection in detections:
    f.write(detection.template_name+', '+str(detection.detect_time)+\
            ', '+str(detection.detect_val)+', '+str(detection.threshold)+\
            ', '+str(detection.no_chans)+'\n')
f.close()
