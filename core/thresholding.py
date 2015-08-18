#!/usr/bin/python
"""
Part of the EQcorrscan package written by Calum Chamberlain of Victoria
University of Wellington in 2015.  The main goal of this package is to generate
waveform templates from continuous data and search for repeats of this waveform
in continuous data using a match-filter cross-correlation routine.

The prupose of this module specifically is to determine thresholds for the use
in various detection routines within the core modules.

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
import numpy as np
from obspy import read

def coherance_test(stream, stations, nodes, lags, wlen):
    """
    Function to determine how the coherance of a day of contious multi-channel,
    multi-station seismic data varies with varying grid nodes.  It would be
    wise to apply this to a day of seismic data thought to be quiet from within
    the grid given, and to a day thought to be loud with energy sources within
    the grid to see how this function varies.

    :type stream: obspy.Stream
    :param stream: Day of multi-channel, multi-station seismic data
    :type stations: list of str
    :param station: List fo the stations in the same order as the delays
    :type nodes: list of tuple
    :param nodes: A list of the node points as (lat, lon, depth)
    :type lags: np.array
    :param lags: Array of arrays where lags[i][:] refers to station[i] for all
                nodes and lags[i][j] refers to the lag for station[i] at node[j]
                should be in seconds.
    :type wlen: int
    :param wlen: Length of window to determine coherance for - should be the
                same as your determined template length - in seconds

    :return: Coherance, array of arrays where coherance[i][:] refers to the
            running daily coherance at node[i]
    """
    from core.bright_lights import coherance
    from copy import deepcopy
    import sys
    # Convert wlen to samples
    wlen = int(wlen * stream[0].stats.sampling_rate)
    #Set up coherance array
    dailycoherance=np.array([np.array([np.nan]*(len(stream[0].data)-wlen))]\
            *len(nodes), dtype=np.float32) # Use half precision as we don't need
                                        # true doubles
    # Keep a sacred copy of the stream, as we will be applying the lags
    # directly to a copy of this
    copyof_stream=deepcopy(stream)
    for i in xrange(len(nodes)):
        sys.stdout.write("\r"+str(float(i)/len(nodes)*100)+"% \r")
        sys.stdout.flush()
        for j in xrange(len(stations)):
            #Apply delays
            st = copyof_stream.select(station = stations[j])
            for tr in st:
                tr.data = tr.data[int(lags[j][i]*tr.stats.sampling_rate):]
                pad = np.zeros(int(lags[j][i]*tr.stats.sampling_rate))
                tr.data = np.concatenate((pad, tr.data))
        for j in xrange(len(dailycoherance[i])):
            # Compute the coherance for each window
            window = deepcopy(copyof_stream)
            for tr in window:
                tr.data = tr.data[j:j+wlen]
            sys.stdout.write("\r"+str((float(i+1)*float(j+1))/(len(nodes)*len(dailycoherance[i]))*100)+"% \r")
            sys.stdout.flush()
            dailycoherance[i][j] = coherance(window)
    return dailycoherance

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        raise IOError("Needs arguments, e.g. --coherance 2009/07/14")
    else:
        if sys.argv[1] == "--coherance":
            date = sys.argv[2]
            import glob
            from obspy import read
            stream = read('test_data/*'+date.split('/')[0]+'-'+\
                    date.split('/')[1]+'-'+date.split('/')[2]+'-processed.ms')
            if not len(stream) == 0:
                import os
                path = os.path.dirname(os.path.abspath(__file__))
                sys.path.insert(0,path[0:len(path)-5])
                from par import template_gen_par as templatedef
                from par import bright_lights_par as brightdef
                from core import bright_lights
                from utils import EQcorrscan_plotting as plotting
                # Use the brightness function to search for possible templates
                # First read in the travel times
                print 'Reading in the original grids'
                stations, allnodes, alllags = \
                        bright_lights._read_tt(brightdef.nllpath,brightdef.stations,\
                                    brightdef.phase)
                # Resample the grid to allow us to run it quickly!
                print 'Cutting the grid'
                stations, nodes, lags = bright_lights._resample_grid(stations, allnodes,
                                                     alllags,
                                                     brightdef.volume,
                                                     brightdef.resolution)
                # Remove lags that have a similar network moveout, e.g. the sum of the
                # differences in moveouts is small.
                print "Removing simlar lags"
                stations, nodes, lags = bright_lights._rm_similarlags(stations, nodes, lags,
                                                      brightdef.nodesimthresh)
                print "Plotting new grid"
                plotting.threeD_gridplot(nodes)
                dailycoherance = coherance_test(stream, stations, nodes, lags, \
                                templatedef.length)
            else:
                raise IOError("No traces read in for this day, are they processed?")
        else:
            raise IOError("I only know --coherance at the moment")
