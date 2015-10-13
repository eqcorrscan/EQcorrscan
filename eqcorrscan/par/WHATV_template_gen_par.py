#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for template generation
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for template generation

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
import glob

sfilebase='/Volumes/GeoPhysics_09/users-data/chambeca/my_programs/Building/EQ_correl_codes/templates/WHATV'
                # Location of nordic s-files, string

sfiles=glob.glob(sfilebase+'/????/??/*L.S*')# List of s-files
# sfiles=[sfile.split('/')[-1] for sfile in sfiles]
#sfiles=[]

# tfiles=glob.glob('templates/pan_templates/brightness_group_*')  # List of pre-existing template files
tfiles=[]
samp_rate=100.0                      # Desired sampling rate, will be carried
                                    # through entire process, float in Hz

lowcut=5.0                          # Lowcut for bandpass filter in Hz, float

highcut=15.0                         # Highcut for bandpass filter in Hz, float

filter_order=3                      # Number of corners for the filter

length=1.0                          # Length of template window after pick

swin='all'                          # use s-picks or not, if True the
                                    # templates will include only channels with
                                    # picks, if False, all channels will be
                                    # used for stations with p-picks.
saveloc='/Volumes/GeoPhysics_09/users-data/chambeca/my_programs/Building/EQ_correl_codes/templates/EQcorrscan'
               # Location to save the templates to in
                                    # miniseed format
debug=3                             # Debugging level from 0-5 with 5 being most
                                    # verbose
