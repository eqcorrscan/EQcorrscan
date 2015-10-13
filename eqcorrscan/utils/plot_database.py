#!/usr/bin/python
"""
Script to plot stations and events from a seisan database.

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

def volume_plot(stationpath, database, limits):
    """
    Function to read in station information from a file and earthquake info
    from sfiles.

    :type stationpath: str
    :type databse: str
    """
    from utils import Sfile_util
    import glob
    sfiles = glob.glob(database+'/*/*/*')
    eqlocs=[]
    for sfile in sfiles:
        try:
            eqlocs+=[(Sfile_util.readheader(sfile).latitude,\
                    Sfile_util.readheader(sfile).longitude,\
                    Sfile_util.readheader(sfile).depth)]
        except:
            continue
    stalocs=[]
    f = open(stationpath, 'r')
    for line in f:
        stalocs+=[(float(line.split(',')[1]),\
                float(line.split(',')[0]), float(line.split(',')[4])/1000)]
    f.close()
    from utils import EQcorrscan_plotting
    EQcorrscan_plotting.threeD_seismplot(stalocs, eqlocs, limits)
    return

if __name__ == "__main__":
    import sys
    if not len(sys.argv) == 3:
        raise IOError("Needs two arguments, station .csv file and database path")
    stationpath = str(sys.argv[1])
    database = str(sys.argv[2])
    import glob
    if not glob.glob(stationpath):
        raise IOError(stationpath+" does not exist")
    if not glob.glob(database):
        raise IOError(database+" does not exist")
    volume_plot(stationpath, database)
