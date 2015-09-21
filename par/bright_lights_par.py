#!/usr/bin/python

#------------------------------------------------------------------------------
#   Purpose:    Outline definitions for match filter python code
#   Author:     Calum John Chamberlain
#------------------------------------------------------------------------------

"""
Outline definitions for bright_lights python code

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
import numpy as np
import matplotlib.path as mplPath
stations=['EORO','WHYM','COSA','FOZ','LBZ','JCZ','WVZ',
          'GOVA','LABE','MTFO','RPZ',
          'COVA','FRAN','POCR','SOLU',
          'WHAT','POCR2','WHAT2']   # List of stations to use for template
                                        # generation

nllpath='./grid/Caro_larger_grid/3D'    # Path to nonlinloc .csv grid files of
                                        # travel times
# corners=mplPath.Path(np.array([[170., -44.],\
                          # [170.1, -44.],\
                          # [170.1, -43.8],\
                          # [170, -43.8]]))
# corners=mplPath.Path(np.array([[169.7, -44.0519],\
                          # [170.3, -43.7],\
                          # [170.3, -43.3],\
                          # [169.4, -43.75]]))
corners=mplPath.Path(np.array([[169.55, -44.2],\
                               [170.8, -43.6],\
                               [170.5, -43.16],\
                               [169.2, -43.76]]))

                                        # Numpy array to be converted into a
                                        # matplotlib path - should descirbe
                                        # the horizontal polygon to search within
                                        # Where coordinates are decimal degrees
                                        # in lat and long
mindepth=15
maxdepth=45
# maxdepth=20                             # Depth cuts in km, cannot use these polygonally
resolution=(0.02,2)                     # Horizontal and vertical resolution
                                        # for resampled grid in decimal
                                        # degrees and km respectively.
                                        # with depth increasing down
dates=[UTCDateTime('2013-03-28')+i \
 for i in xrange(0, int(UTCDateTime('2013-05-09') - UTCDateTime('2013-03-28')),\
                 86400)] # Example of a generator expression for all-time
trem_dates=[UTCDateTime('2009-05-12')]#, UTCDateTime('2009-07-14'),\
       # UTCDateTime('2009-07-15'), UTCDateTime('2009-11-15'),\
       # UTCDateTime('2010-07-05'), UTCDateTime('2010-07-14'),\
       # UTCDateTime('2010-08-20'), UTCDateTime('2010-08-31'),\
       # UTCDateTime('2010-10-05'), UTCDateTime('2011-08-03'),\
       # UTCDateTime('2011-09-02'), UTCDateTime('2011-09-04'),\
       # UTCDateTime('2013-03-28')] #tremor days
# dates=trem_dates+dates
dates=trem_dates
                                        # List of dates to run through, can be
                                        # made in any pythonic way, but must be
                                        # a list of obspy.UTCDateTime objects
threshold=1.5                           # Threshold value, if threshtype is
                                        # set to MAD then this will by the
                                        # MAD to the power of this value
thresh_type='RMS'                       # Threshold type, can be either MAD
                                        # or abs (case-sensitive)
phase='P'                               # Can be either P or S, will use this
                                        # when reading in the travel-time
                                        # grids.
ps_ratio=1.68                           # P to S ratio if only one grid type
                                        # provided, P(vel) is S(vel)*ps_ratio
nodesimthresh=0.5                       # Minimum cumulative difference in
# nodesimthresh=0.0625                    # Minimum cumulative difference in
                                        # network moveout, should be about the
                                        # period of twice the maximum frequency
                                        # of the signal you want to detect
coherance=0.12                          # Coherance threshold to remove
                                        # incoherant peaks in the network
                                        # response.
plotsave=True                           # Save plots, set to True to not show
                                        # any plots
clip_level=1.5                          # Level to clip energy amplitudes to as
                                        # a multiplier of the mean energy amplitude
cores=20                                 # Number of cores to use for brightness search
