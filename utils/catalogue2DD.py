#!/usr/bin/python
"""
Module written by Calum Chamberlain as part of the EQcorrscan package.

This module contains functions to convert a seisan catalogue to files ready for
relocation in hypoDD - it will generate both a catalogue (dt.ct) file, event
file (event.dat), station information file (station.dat), and a correlation
oiutput file correlated every event in the catalogue with every other event to
optimize the picks (dt.cc).

The correlation routine relies on obspy's xcorrPickCorrection function from the
obspy.signal.cross_correlation module.  This function optimizes picks to better
than sample accuracy by interpolating the correlation function and finding the
maximum of this rather than the true maximum correlation value.  The output
from this function is stored in the dt.cc file.

Information for the station.dat file is read from SEISAN's STATION0.HYP file

Earthquake picks and locations are taken from the catalogued s-files - these
must be pre-located before entering this routine as origin times and hypocentre
locations are needed for event.dat files.
"""
from utils import Sfile_util
def readSTATION0(path):
    """
    Function to read the STATION0.HYP file on the path given.

    :type path: String
    :param path: Path to the STATION0.HYP file

    :returns: List of tuples of station, lat, long, elevation
    """
    stalist=[]
    f=open(path,'r')
    return stalist

def write_event(sfile_list):
    """
    Function to write out an event.dat file of the events

    :type sfile_list: List
    :param sfile_list: List of s-files to sort and put into the database

    :returns: List of tuples of event ID (int) and Sfile name
    """
    event_list=[]
    sfile_list.sort()
    i=0
    f=open('event.dat','w')
    for sfile in sfile_list:
        i+=1
        event_list.append(i, sfile)
        ev_info=Sfile_util.readheader(sfile)
        f.write(ev_info.time+'\n')
    return event_list
