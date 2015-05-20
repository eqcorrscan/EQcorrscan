#!/usr/bin/python

# simple script to check if COSA amplitude picks are useful

import glob
from utils import readSfile
max_period=0.025
sfiles=glob.glob('../COSAL/*/*/*L.S*')
for sfile in sfiles:
    # print 'Reading '+sfile
    picks=readSfile.readpicks(sfile)
    for pick in picks:
        if pick.peri != 'nan' and pick.peri < max_period:
            print sfile+' '+pick.station+' may be dodgy'

