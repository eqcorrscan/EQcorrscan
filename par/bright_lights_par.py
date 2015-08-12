#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for match filter python code
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for bright_lights python code
"""
from obspy import UTCDateTime
import numpy as np
import matplotlib.path as mplPath
stations=['EORO','WHYM','COSA','FOZ',
          'GOVA','LABE','MTFO','RPZ',
          'COVA','FRAN','POCR','SOLU',
          'WHAT','POCR2','WHAT2']   # List of stations to use for template
                                        # generation

nllpath='./grid/3D'                       # Path to nonlinloc .csv grid files of
                                        # travel times
corners=mplPath.Path(np.array([[169.7, -44.0519],\
                           [170.3, -43.7],\
                           [170.3, -43.3],\
                           [169.4, -43.75]]))
                                        # Numpy array to be converted into a
                                        # matplotlib path - should descirbe
                                        # the horizontal polygon to search within
                                        # Where coordinates are decimal degrees
                                        # in lat and long
mindepth=15
maxdepth=45                             # Depth cuts in km, cannot use these polygonally
resolution=(0.02,2)                     # Horizontal and vertical resolution
                                        # for resampled grid in decimal
                                        # degrees and km respectively.
                                        # with depth increasing down
dates=[UTCDateTime('2009-03-26')+i \
 for i in xrange(0, int(UTCDateTime('2015-03-09') - UTCDateTime('2009-03-26')),\
                 86400)] # Example of a generator expression for all-time
# dates=[UTCDateTime('2009-05-12'), UTCDateTime('2009-07-14'),\
       # UTCDateTime('2009-07-15'), UTCDateTime('2009-11-15'),\
       # UTCDateTime('2010-07-05'), UTCDateTime('2010-07-14'),\
       # UTCDateTime('2010-08-20'), UTCDateTime('2010-08-31'),\
       # UTCDateTime('2010-10-05'), UTCDateTime('2011-08-03'),\
       # UTCDateTime('2011-09-02'), UTCDateTime('2011-09-04'),\
       # UTCDateTime('2013-03-28')] #tremor days
                                        # List of dates to run through, can be
                                        # made in any pythonic way, but must be
                                        # a list of obspy.UTCDateTime objects
threshold=10                            # Threshold value, if threshtype is
                                        # set to MAD then this will by the
                                        # MAD to the power of this value
thresh_type='RMS'                       # Threshold type, can be either MAD
                                        # or abs (case-sensitive)
phase='P'                               # Can be either P or S, will use this
                                        # when reading in the travel-time
                                        # grids.
ps_ratio=1.68                           # P to S ratio if only one grid type
                                        # provided, P(vel) is S(vel)*ps_ratio
nodesimthresh=0.0625                    # Minimum cumulative difference in
                                        # network moveout, should be about the
                                        # period of twice the maximum frequency
                                        # of the signal you want to detect
coherance=0.10                          # Coherance threshold to remove
                                        # incoherant peaks in the network
                                        # response.
plotsave=True                           # Save plots, set to True to not show
                                        # any plots
clip_level=8                            # Level to clip energy amplitudes to as
                                        # a multiplier of the mean energy amplitude
cores=40                                # Number of cores to use for brightness search
