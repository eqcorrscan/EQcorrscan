r"""Tutorial  This script is designed as a tutorial to highlight how to \
call the main functions within the EQcorrscan module.

In this tutorial
we will see how to generate a template and run this through the
matched-filter routine.
The template will be generated from a pre-picked earthquake, however there
are other ways to generate templates, for example this package also contains
a simple brightness function that is designed to scan continuous seismic
data and look for impulsive energy originating from a discrete point in a
seismic velocity model.

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
import os
from obspy import read
from eqcorrscan.core import template_gen, match_filter
# Before calling these module imports for parameter files you should insert
# your own path into sys.path so that we find your parameter files.
from eqcorrscan.utils import pre_processing, Sfile_util
import glob

# Set up the default parameters - these used to be stored in parameter files
debug = 3  # High debug level should output lots to keep you informed
threshold = 8.0  # Threshold level as MAD multiplier
threshtype = 'MAD'  # Threshold type, in this case Median Absolute Deviation
trig_int = 6.0  # Minimum trigger interval for one template in seconds

# Now we find the s-file we want to use to generate a template from
data_directory = os.path.join('test_data', 'tutorial_data')
sfiles = glob.glob(os.path.join(data_directory, '*L.S*'))
print sfiles

templates = []
template_names = []
for i, sfile in enumerate(sfiles):
    # Read in the picks from the S-file, note, in the full case one fo the main
    # functions in template_gen would be used rather than this, but for
    # the tutorial we will read in the data here - also note that this
    # template generation is inefficient for multiple templates, if using
    # daylong data for multiple templates you would want to only read
    # the seismic data once and cut it multiple times.
    event = Sfile_util.readpicks(sfile)
    for pick in event.picks:
        print pick
        if 'wavefiles' not in locals():
            wavefiles = glob.glob(os.path.join(data_directory,
                                               '.'.join([pick.waveform_id.
                                                         station_code, '*'])))
        else:
            wavefiles += glob.glob(os.path.join(data_directory,
                                                '.'.join([pick.waveform_id.
                                                          station_code, '*'])))
    wavefiles = list(set(wavefiles))
    for wavefile in wavefiles:
        print ' '.join(['Reading data from', wavefile])
        if 'st' not in locals():
            st = read(wavefile)
        else:
            st += read(wavefile)

    st = st.merge(fill_value='interpolate')
    day = st[0].stats.starttime.date

    # Process the data with our required parameters
    for tr in st:
        tr = pre_processing.dayproc(tr, 1.0, 20.0, 3, 100.0,
                                    debug, day)

    # Use the template generation function to cut our templates
    template = template_gen._template_gen(event.picks, st, length=1.0,



                                          swin='all', prepick=0.1, plot=True)
    # This will generate an obspy.Stream object
    # Append this Stream to the list of templates
    templates += [template]
    template_names.append('_'.join(['tutorial', str(i)]))

    # Save template for later
    template.write(os.path.join(data_directory, '_'.join([template_names[i],
                                                          'template.ms'])),
                   format='MSEED')
    # Delete excess information from memory If you are re-using this script
    # with the same templates you should be able to comment out this loop
    # once you have generated your templates once.
    del template, st

# Extract the stations from the templates
for template in templates:
    del tr
    if 'stachans' not in locals():
        stachans = [(tr.stats.station, tr.stats.channel) for tr in template]
    else:
        stachans += [(tr.stats.station, tr.stats.channel) for tr in template]

# Make this a unique list
stachans = list(set(stachans))

# Read in the continuous data for these station, channel combinations
for stachan in stachans:
    data_file = ''.join([stachan[0], '.*..*', stachan[1][-1], '.*'])
    data_file = os.path.join(data_directory, data_file)
    print ' '.join(['Reading data from:', data_file])
    # Generate a new stream object and add to it
    if 'st' not in locals():
        st = read(data_file)
    else:
        st += read(data_file)

# Merge the data to account for miniseed files being written in chunks
# We need continuous day-long data, so data are padded if there are gaps
st = st.merge(fill_value='interpolate')

# Work out what day we are working on, required as
# we will pad the data to be daylong
day = st[0].stats.starttime.date

# Process the data in the same way as the template
for tr in st:
    tr = pre_processing.dayproc(tr, 1.0, 20.0, 3, 100.0,
                                debug, day)

# Compute detections
detections = match_filter.match_filter(template_names, templates, st,
                                       threshold, threshtype, trig_int,
                                       plotvar=True, cores=2, tempdir=False,
                                       debug=debug, plot_format='pdf')

# We now have a list of detections! We can output these to
# a file to check later
f = open('tutorial_detections.csv', 'w')
for detection in detections:
    line = ', '.join([detection.template_name, str(detection.detect_time),
                      str(detection.detect_val), str(detection.threshold),
                      str(detection.no_chans)])
    f.write(line)
    print line
    f.write(os.linesep)
f.close()
