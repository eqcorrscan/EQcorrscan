#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for match filter python code
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for bright_lights python code
"""
from obspy import UTCDateTime
stations=['EORO','WHYM','COSA','FOZ',
          'GOVA','LABE','MTFO','RPZ']   # List of stations to use for template
                                        # generation

nllpath='./grid'                       # Path to nonlinloc .csv grid files of
                                        # travel times
volume=[(-43.8,-43.55),
        (169.5,170.05),
        (15,40)]
# volume=[(-43.8,-43.75),(169.5,169.52),(18,20)]
                                        # List of tuples in the form:
                                        # [(minlat,maxlat),(minlong,laxlong),
                                        # (mindepth,maxdepth)]
                                        # Where coordinates are decimal degrees
                                        # in lat and long and km in depth
resolution=(0.02,2)                     # Horizontal and vertical resolution
                                        # for resampled grid in decimal
                                        # degrees and km respectively.
                                        # with depth increasing down
startdate=UTCDateTime('2009-07-14')     # obspy UTCDateTime object giving
                                        # date to start looking for templates
                                        # beyond
enddate=UTCDateTime('2009-07-15')       # As above, but the end date
threshold=200                           # Threshold value, if threshtype is
                                        # set to MAD then this will by the
                                        # MAD to the power of this value
thresh_type='MAD'                       # Threshold type, can be either MAD
                                        # or abs (case-sensitive)
phase='S'                               # Can be either P or S, will use this
                                        # when reading in the travel-time
                                        # grids.
nodesimthresh=0.0625                    # Minimum cumulative difference in
                                        # network moveout, should be about the
                                        # period of twice the maximum frequency
                                        # of the signal you want to detect
coherance=0.2                           # Coherance threshold to remove
                                        # incoherant peaks in the network
                                        # response.
