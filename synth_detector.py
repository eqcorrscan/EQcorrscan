"""
Script to utilise the match_filter routines of EQcorrscan and the synth_seis
routines in the utils of EQcorrscan.  This script will read in a grid of
travel-times and generate a series of synthetic templates for this grid.
It will then use this selection of synthetic templates in a match-filter
routine to detect similar, real-earthquakes in coninuous seismic data.

These detected events can then be stacked and used as further templates...
ad infinitium...
"""

import sys, glob, datetime as dt, numpy as np
instance=0
Split=False
startdate=False
parallel=True
if len(sys.argv) == 2:
    flag=str(sys.argv[1])
    if flag == '--debug':
        Test=True
        Prep=False
    elif flag == '--debug-prep':
        Test=False
        Prep=True
    else:
        raise IOError("I don't recognise your arguments I know --debug and "+\
        "--debug-prep, and --instance, --splits, --startdate, --enddate")
elif len(sys.argv) == 5:
    args=sys.argv[1:len(sys.argv)]
    if args[0] == '--instance' or args[2]=='--instance':
        # Arguments to allow the code to be run in multiple instances
        Split=True
        Test=False
        Prep=False
        for i in xrange(len(args)):
            if args[i] == '--instance':
                instance=int(args[i+1])
                print 'I will run this for instance '+str(instance)
            elif args[i] == '--splits':
                splits=int(args[i+1])
                print 'I will divide the days into '+str(splits)+' chunks'
    elif args[0] == '--startdate' or args[2] == '--startdate':
        Split=False
        Test=False
        Prep=False
        for i in xrange(len(args)):
            if args[i]== '--startdate':
                startdate=dt.datetime.strptime(str(args[i+1]),'%Y/%m/%d')
                print 'Start date is: '+dt.datetime.strftime(startdate, '%Y/%m/%d')
            elif args[i] == '--enddate':
                enddate=dt.datetime.strptime(str(args[i+1]), '%Y/%m/%d')
                print 'End date is: '+dt.datetime.strftime(enddate, '%Y/%m/%d')
    else:
        raise IOError("I don't recognise your arguments,  I know --debug and --debug-prep, and --instance, --splits, --startdate, --enddate")
elif len(sys.argv) == 7:
    args=sys.argv[1:len(sys.argv)]
    Split=False
    Test=False
    Prep=False
    for i in xrange(len(args)):
        if args[i]== '--startdate':
            startdate=dt.datetime.strptime(str(args[i+1]),'%Y/%m/%d')
            print 'Start date is: '+dt.datetime.strftime(startdate, '%Y/%m/%d')
        elif args[i] == '--enddate':
            enddate=dt.datetime.strptime(str(args[i+1]), '%Y/%m/%d')
            print 'End date is: '+dt.datetime.strftime(enddate, '%Y/%m/%d')
        elif args[i] == '--instance':
            instance=int(args[i+1])
            print 'This instance is given the flag '+str(instance)
elif not len(sys.argv) == 1:
    raise ValueError("I only take one argument, no arguments, or two flags with arguments")
else:
    Test=False
    Prep=False
    Split=False

from par import template_gen_par as templatedef
from par import match_filter_par as matchdef
from par import bright_lights_par as brightdef
from utils import seismo_logs
if brightdef.plotsave:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.ioff()
#from par import lagcalc as lagdef
from obspy import UTCDateTime, Stream, read as obsread
# First generate the templates
from core import bright_lights, match_filter
from utils import pre_processing, synth_seis
from utils import EQcorrscan_plotting as plotting
from obspy.signal.filter import bandpass
from joblib import Parallel, delayed
import warnings

templates=[]
delays=[]
stations=[]
print 'Template generation parameters are:'
print '     sfilebase: '+templatedef.sfilebase
print '     samp_rate: '+str(templatedef.samp_rate)+' Hz'
print '     lowcut: '+str(templatedef.lowcut)+' Hz'
print '     highcut: '+str(templatedef.highcut)+' Hz'
print '     length: '+str(templatedef.length)+' s'
print '     swin: '+templatedef.swin+'\n'


# Use the brightness function to search for possible templates
# First read in the travel times
print 'Reading in the original grids'
stations, allnodes, alltravel_times = \
            bright_lights._read_tt(brightdef.nllpath,brightdef.stations,\
                                    brightdef.phase, phaseout='S', \
                                    ps_ratio=brightdef.ps_ratio, lags_switch=False)
print 'I have read in '+str(len(allnodes))+' nodes'

# We now have a grid of travel-times which can then be used to generate synthetic\
# seismograms using the utils.synth_seis functions.

# We should trim the grid to the area we want to work in
print 'Cutting the grid'
stations, nodes, travel_times = bright_lights._resample_grid(stations, allnodes,
                                                     alltravel_times,
                                                     brightdef.mindepth,
                                                     brightdef.maxdepth,
                                                     brightdef.corners,
                                                     brightdef.resolution)
del allnodes, alltravel_times
# Check that we still have a grid!
if len(nodes) == 0:
    raise IOError("You have no nodes left")

# Call the template generation function
synth_templates=synth_seis.template_grid(stations, nodes, travel_times, 'S', \
        PS_ratio=brightdef.ps_ratio, samp_rate=templatedef.samp_rate)

# Write out the synthetics!
i=0
for synth in synth_templates:
    # We need the data to be in int32
    for tr in synth:
        tr.filter('bandpass', freqmin=templatedef.lowcut,\
                    freqmax=templatedef.lowcut)
        tr.data=(tr.data*1000).astype(np.int32)
    synth.write('templates/synthetics/'+str(nodes[i][0])+'_'+str(nodes[i][1])+\
                '_'+str(nodes[i][2])+'_template.ms', format='MSEED',\
                encoding='STEIM2', reclen=512)
    i+=1
