#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for template generation
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for template generation
"""
import glob

sfilebase='.'                       # Location of nordic s-files, string

sfiles=[]                           # List of s-files

tfiles=glob.glob('stack_templates/brightness_group_*')  # List of pre-existing template files

samp_rate=20.0                      # Desired sampling rate, will be carried
                                    # through entire process, float in Hz

lowcut=2.0                          # Lowcut for bandpass filter in Hz, float

highcut=8.0                         # Highcut for bandpass filter in Hz, float

filter_order=3                      # Number of corners for the filter

length=3.0                          # Length of template window after pick

swin='all'                          # Boolean, use s-picks or not, if True the
                                    # templates will include only channels with
                                    # picks, if False, all channels will be
                                    # used for stations with p-picks.
saveloc='./templates/brightness_alltime'               # Location to save the templates to in
                                    # miniseed format
debug=3                             # Debugging level from 0-5 with 5 being most
                                    # verbose
