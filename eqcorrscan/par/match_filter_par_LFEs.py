#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for match filter python code
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for match filter python code

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

from obspy import UTCDateTime

# dates=[UTCDateTime('2010-12-03')+i\
       # for i in xrange(0,int(UTCDateTime('2011-03-15') - UTCDateTime('2010-12-03')),\
                       # 86400)]
# dates+=[UTCDateTime('2013-09-05')+i\
       # for i in xrange(0,int(UTCDateTime('2013-12-17') - UTCDateTime('2013-09-05')),\
                       # 86400)]


# dates=[UTCDateTime('2009-03-26')+i \
 # for i in xrange(0, int(UTCDateTime('2015-03-09') - UTCDateTime('2009-03-26')),\
                 # 86400)] # Example of a generator expression for all-time
dates=[UTCDateTime('2009-05-12'), UTCDateTime('2009-07-14'),\
       UTCDateTime('2009-07-15'), UTCDateTime('2009-11-15'),\
       UTCDateTime('2010-07-05'), UTCDateTime('2010-07-14'),\
       UTCDateTime('2010-08-20'), UTCDateTime('2010-08-31'),\
       UTCDateTime('2010-10-05'), UTCDateTime('2011-08-03'),\
       UTCDateTime('2011-09-02'), UTCDateTime('2011-09-04'),\
       UTCDateTime('2013-03-28')] #tremor days
                                        # List of dates to run through, can be
                                        # made in any pythonic way, but must be
                                        # a list of obspy.UTCDateTime objects

threshold=8.0                       # Threshold to declare a detection at, type
                                    # of threshold is set below

threshtype='MAD'                    # Type of threshold, can be either MAD or
                                    # absolute or av_chan_corr:
                                    # absolute is an absolute correlation
                                    # threshold, which will be applied to the
                                    # network correlation sum independent on
                                    # the number of stations used
                                    # MAD is a Median Absolute Deviation based
                                    # threshold, in that the threshold is set
                                    # at 'threshold'*MAD, where MAD is the Median
                                    # Absolute Deviation in the summed
                                    # correlation sum - this method is
                                    # recommended over ABS.
                                    # av_chan_corr will detect based on the
                                    # average single channel correlation value
                                    # assuming all stations are present.

trig_int=6                          # Trigger interval, determines how often to
                                    # make detections.  this value should be in
                                    # seconds.

minsta=3                            # Minimum number of stations to run the
                                    # detection routine for
contbase=[('/Volumes/GeoPhysics_09/users-data/chambeca/SAMBA_archive/day_volumes_S',\
          'Yyyyy/Rjjj.01','AF'),
          ('/Volumes/GeoPhysics_09/users-data/chambeca/GeoNet_archive/day_volumes_S',\
         'Yyyyy/Rjjj.01','NZ')]
          # ('/Volumes/GeoPhysics_09/users-data/chambeca/SAMBA_archive/day_volumes_S',\
          # 'Yyyyy/Rjjj.01','NZ')]

	# ('/Volumes/GeoPhysics_09/users-data/chambeca/Alpine_Fault_SAC/SAC_resampled',\
	# 'yyyymmdd','AF'),
        # ('/Volumes/GeoPhysics_09/users-data/chambeca/Alpine_Fault_SAC/SAC_resampled',\
         # 'yyyymmdd','NZ')]
                                    # Full path for the waveform database
                                    # Files must be in daylong format
                                    # To allow for data from multiple directories
                                    # for different networks this must be a list
                                    # of tuples of the form
                                    # [(path, directory format, netcode)]

plot=False                          # boolean, True for plotting the daily
                                    # cccsum values for each template.
cores=40                            # Value for number of parallel jobs to run
                                    # must be int
debug=1                             # Debug level, 0-5 with 0 being less verbose
