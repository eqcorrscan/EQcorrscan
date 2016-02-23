r"""
The prupose of this module specifically is to determine thresholds for the \
use in various detection routines within the core modules.

** Unfinished **

"""
import numpy as np


def coherence_test(stream, stations, nodes, lags, wlen):
    """
    Function to determine how the coherence of a day of contious \
    multi-channel, multi-station seismic data varies with varying grid nodes.

    It would be
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
    :param lags: Array of arrays where lags[i][:] refers to station[i] for \
        all nodes and lags[i][j] refers to the lag for station[i] at node[j] \
        should be in seconds.
    :type wlen: int
    :param wlen: Length of window to determine coherence for - should be the \
        same as your determined template length - in seconds.

    :return: coherence, array of arrays where coherence[i][:] refers to the \
        running daily coherence at node[i].
    """
    from core.bright_lights import coherence
    from copy import deepcopy
    import sys
    # Convert wlen to samples
    wlen = int(wlen * stream[0].stats.sampling_rate)
    # Set up coherence array
    dailycoherence = np.array([np.array([np.nan]*(len(stream[0].data)-wlen))]
                              * len(nodes), dtype=np.float32)
    # Use half precision as we don't need true doubles
    # Keep a sacred copy of the stream, as we will be applying the lags
    # directly to a copy of this
    copyof_stream = deepcopy(stream)
    for i in xrange(len(nodes)):
        sys.stdout.write("\r"+str(float(i)/len(nodes)*100)+"% \r")
        sys.stdout.flush()
        for j in xrange(len(stations)):
            # Apply delays
            st = copyof_stream.select(station=stations[j])
            for tr in st:
                tr.data = tr.data[int(lags[j][i]*tr.stats.sampling_rate):]
                pad = np.zeros(int(lags[j][i]*tr.stats.sampling_rate))
                tr.data = np.concatenate((pad, tr.data))
        for j in xrange(len(dailycoherence[i])):
            # Compute the coherence for each window
            window = deepcopy(copyof_stream)
            for tr in window:
                tr.data = tr.data[j:j+wlen]
            sys.stdout.write("\r"+str((float(i + 1) * float(j + 1)) /
                                      (len(nodes) * len(dailycoherence[i])) *
                                      100)+"% \r")
            sys.stdout.flush()
            dailycoherence[i][j] = coherence(window)
    return dailycoherence

# if __name__ == '__main__':
#     import sys
#     if len(sys.argv) == 1:
#         raise IOError("Needs arguments, e.g. --coherence 2009/07/14")
#     else:
#         if sys.argv[1] == "--coherence":
#             date = sys.argv[2]
#             stream = read('test_data/*'+date.split('/')[0]+'-' +
#                           date.split('/')[1]+'-'+date.split('/')[2] +
#                           '-processed.ms')
#             if not len(stream) == 0:
#                 import os
#                 path = os.path.dirname(os.path.abspath(__file__))
#                 sys.path.insert(0, path[0:len(path)-5])
#                 from core import bright_lights
#                 from utils import EQcorrscan_plotting as plotting
#                 # Use the brightness function to search for possible
#                 # templates
#                 # First read in the travel times
#                 print 'Reading in the original grids'
#                 stations, allnodes, alllags = \
#                     bright_lights._read_tt(nllpath, stations,
#                                            brightdef.phase)
#                 # Resample the grid to allow us to run it quickly!
#                 print 'Cutting the grid'
#                 stations, nodes, lags = \
#                     bright_lights._resample_grid(stations, allnodes,
#                                                  alllags,
#                                                  brightdef.volume,
#                                                  brightdef.resolution)
#                 # Remove lags that have a similar network moveout, e.g.
#                 # the sum of the differences in moveouts is small.
#                 print "Removing simlar lags"
#                 stations, nodes, lags = \
#                     bright_lights._rm_similarlags(stations, nodes, lags,
#                                                   brightdef.nodesimthresh)
#                 print "Plotting new grid"
#                 plotting.threeD_gridplot(nodes)
#                 dailycoherence = coherence_test(stream, stations, nodes,
#                                                 lags,
#                                                 templatedef.length)
#             else:
#                 raise IOError("No traces read in for this day, are they " +
#                               "processed?")
#         else:
#             raise IOError("I only know --coherence at the moment")
