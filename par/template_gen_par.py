#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for template generation
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for template generation
"""

sfilebase='.'                       # Location of nordic s-files, string

sfiles=['63877_1_2_matlab_renamed',
        '17208_1_2_3_matlab_renamed',
        '60070_1_2_matlab_renamed',
        '37575_1_2_matlab_renamed',
        '55115_1_1_matlab_renamed',
        '60905_1_matlab_renamed',
        '61044_1_2_3_matlab_renamed',
        '54905_1_2_3_4_5_matlab_renamed',
        '61220_1_2_3_matlab_renamed',
        '30442_1_2_matlab_renamed',
        '55432_1_2_3_4_matlab_renamed',
        '55200_1_2_3_4_5_matlab_renamed',
        '61100_1_matlab_renamed',
        '59966_1_2_matlab_renamed'] # Names of s-files to use, list of strings,
                                    # or list of already generated templates

samp_rate=100.0                     # Desired sampling rate, will be carried
                                    # through entire process, float in Hz

lowcut=2.0                          # Lowcut for bandpass filter in Hz, float

highcut=8.0                         # Highcut for bandpass filter in Hz, float

filter_order=3                      # Number of corners for the filter

length=6.0                          # Length of template window after pick

swin='all'                          # Boolean, use s-picks or not, if True the
                                    # templates will include only channels with
                                    # picks, if False, all channels will be
                                    # used for stations with p-picks.
saveloc='./templates'               # Location to save the templates to in
                                    # miniseed format
debug=4                             # Debugging level from 0-5 with 5 being most
                                    # verbose
