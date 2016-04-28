"""
Functions for common operations on whole obspy.Inventory classes

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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import warnings
import seaborn as sns
import matplotlib.pyplot as plt


def inv2nlloc(inventory, out_file='tmp.txt'):
    r"""
    Function to write an obspy.Inventory class to the correct format for\
    NLLoc Grid2Time.

    Note: Not very general. Requires coords in LATLON currently and output\
    needs to be manually put into nlloc.in file. Could be developed more.
    """
    #Loop over each station and write the NLLOC accepted format
    import csv
    with open(out_file, 'wb') as f:
        csvwriter = csv.writer(f, delimiter=' ', escapechar=' ',
                               quoting=csv.QUOTE_NONE)
        for net in inventory:
            for sta in net:
                name = str(sta.code)
                lat = str(sta.latitude)
                lon = str(sta.longitude)
                #Elevation in km
                elev = sta.elevation / 1000
                # Account for borehole depths
                depth = sta.channels[0].depth
                elev = str(elev - depth)
                #Not entirely sure why adding the ' ' to this produced 3 spaces...
                csvwriter.writerow(['GTSRCE', name, 'LATLON', lat, lon, '0',
                                    ' ', elev])


def calc_sta_spacing(inventory, method='average'):
    """
    Calculate the station spacing for an obspy.Inventory class. Will return\
    either a float or a list of floats depending upon which method is\
    specified.

    Note: function will ignore separate networks and compute distances for all\
    stations in the inventory.
    """
    from eqcorrscan.utils.mag_calc import dist_calc

    # Merge all station into one list for ease of iteration
    sta_list = []
    for net in inventory:
        for sta in net:
            sta_list.append(sta)
    # Now compute average station spacing
    dist_list = []
    for i, master_sta in enumerate(sta_list):
        master_tup = (master_sta.latitude, master_sta.longitude,
                      master_sta.elevation // 1000 -
                      master_sta.channels[0].depth // 1000)
        # Loop over remaining station pairs
        for slave_sta in sta_list[i+1:]:
            slave_tup = (slave_sta.latitude, slave_sta.longitude,
                         slave_sta.elevation // 1000 -
                         slave_sta.channels[0].depth // 1000)
            dist_list.append(dist_calc(master_tup, slave_tup))
    if method == 'average':
        return np.mean(dist_list)
    elif method == 'all':
        return dist_list
