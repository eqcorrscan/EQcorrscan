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
import Sfile_util
import numpy as np

def _cc_round(num, dp):
    """
    Convenience function to take a float and round it to dp padding with zeros
    to return a string

    :type num: float
    :type dp: int

    :returns: string
    """
    num=round(num, dp)
    num='{0:.{1}f}'.format(num, dp)
    return num

def _av_weight(W1, W2):
    """
    Function to convert from two seisan weights (0-4) to one hypoDD weight(0-1)

    :type W1: str
    :type W2: str

    :returns: str
    """
    if W1==' ':
        W1=1
    else:
        W1=1-int(W1)/4.0
    if W2==' ':
        W2=1
    else:
        W2=1-int(W2)/4.0
    W=(W1+W2)/2
    return _cc_round(W,4)

def _separation(loc1, loc2):
    """
    Function to calculate the distance between two points in the earth

    :type loc1: tuple (float, float, float)
    :param loc1: First point location as lat, long, depth in deg, deg, km
    :type loc2: tuple (float, float, float)
    :param loc2: First point location as lat, long, depth in deg, deg, km

    :returns: distance in km (float)
    """
    R=6371.009  # Radius of the Earth in km
    dlat=np.radians(abs(loc1[0]-loc2[0]))
    dlong=np.radians(abs(loc1[1]-loc2[1]))
    ddepth=abs(loc1[2]-loc2[2])
    mean_lat=np.radians((loc1[0]+loc2[0])/2)
    dist=R*np.sqrt(dlat**2+(np.cos(mean_lat)*dlong)**2)
    dist=np.sqrt(dist**2+ddepth**2)
    return dist

def readSTATION0(path, stations):
    """
    Function to read the STATION0.HYP file on the path given.  Outputs written
    in station.dat file.

    :type path: String
    :param path: Path to the STATION0.HYP file
    :type station: List
    :param station: Stations to look for

    :returns: List of tuples of station, lat, long, elevation
    """
    stalist=[]
    f=open(path+'/STATION0.HYP','r')
    for line in f:
        if line[2:6].strip() in stations:
            station=line[2:6].strip()
            print station
            lat=line[6:14] # Format is either ddmm.mmS/N or ddmm(.)mmmS/N
            if lat[-1]=='S':
                NS=-1
            else:
                NS=1
            if lat[4]=='.':
                lat=(int(lat[0:2])+float(lat[2:-1])/60)*NS
            else:
                lat=(int(lat[0:2])+float(lat[2:4]+'.'+lat[4:-1])/60)*NS
            lon=line[14:23]
            if lon[-1]=='S':
                EW=-1
            else:
                EW=1
            if lon[5]=='.':
                lon=(int(lon[0:3])+float(lon[3:-1])/60)*EW
            else:
                lon=(int(lon[0:3])+float(lon[3:5]+'.'+lon[5:-1])/60)*EW
            elev=float(line[23:-1].strip())
            stalist.append((station, lat, lon, elev))
    f.close()
    f=open('station.dat','w')
    for sta in stalist:
        f.write(sta[0]+'   '+_cc_round(sta[1],4)+'   '+_cc_round(sta[2],4)+\
                '   '+_cc_round(sta[3],4)+'\n')
    f.close()
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
        event_list.append((i, sfile))
        evinfo=Sfile_util.readheader(sfile)
        f.write(str(evinfo.time.year)+str(evinfo.time.month).zfill(2)+\
                str(evinfo.time.day).zfill(2)+'  '+\
                str(evinfo.time.hour).rjust(2)+str(evinfo.time.minute).zfill(2)+\
                str(evinfo.time.second).zfill(2)+\
                str(evinfo.time.microsecond).zfill(2)+'  '+\
                str(evinfo.latitude).ljust(8,'0')+'   '+\
                str(evinfo.longitude).ljust(8,'0')+'  '+\
                str(evinfo.depth).rjust(7).ljust(9,'0')+'   '+\
                str(evinfo.Mag_1)+'    0.00    0.00   '+\
                str(evinfo.t_RMS).ljust(4,'0')+\
                str(i).rjust(11)+'\n')
    f.close()
    return event_list

def write_catalogue(event_list, max_sep=10, min_link=8):
    """
    Function to write the dt.ct file needed by hypoDD - takes input event list
    from write_event as a list of tuples of event id and sfile.  It will read
    the pick information from the seisan formated s-file using the Sfile_util
    utilities.

    :type event_list: List of tuple
    :param event_list: List of tuples of event_id (int) and sfile (String)
    :type max_sep: float
    :param max_sep: Maximum seperation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired

    :returns: List stations
    """
    f=open('dt.ct','w')
    stations=[]
    for i in xrange(len(event_list)):
        master_sfile=event_list[i][1]
        master_event_id=event_list[i][0]
        master_picks=Sfile_util.readpicks(master_sfile)
        master_ori_time=Sfile_util.readheader(master_sfile).time
        master_location=(Sfile_util.readheader(master_sfile).latitude,\
                         Sfile_util.readheader(master_sfile).longitude,\
                         Sfile_util.readheader(master_sfile).depth)
        for j in xrange(i+1,len(event_list)):
            # Use this tactic to only output unique event pairings
            slave_sfile=event_list[j][1]
            slave_event_id=event_list[j][0]
            # Write out the header line
            event_text='#'+str(master_event_id).rjust(10)+\
                    str(slave_event_id).rjust(10)+'\n'
            slave_picks=Sfile_util.readpicks(slave_sfile)
            slave_ori_time=Sfile_util.readheader(slave_sfile).time
            slave_location=(Sfile_util.readheader(slave_sfile).latitude,\
                         Sfile_util.readheader(slave_sfile).longitude,\
                         Sfile_util.readheader(slave_sfile).depth)
            if _separation(master_location, slave_location) > max_sep:
                break
            links=0 # Count the number of linkages
            for pick in master_picks:
                if pick.phase not in ['P','S']:
                    continue # Only use P and S picks, not amplitude or 'other'
                # Find station, phase pairs
                slave_matches=[p for p in slave_picks if p.station==pick.station\
                               and p.phase==pick.phase]
                # Loop through the matches
                for slave_pick in slave_matches:
                    links+=1
                    event_text+=pick.station.rjust(4)+\
                            _cc_round(pick.time-master_ori_time,3).rjust(11)+\
                            _cc_round(slave_pick.time-slave_ori_time,3).rjust(8)+\
                            _av_weight(pick.weight, slave_pick.weight)+' '+\
                            pick.phase+'\n'
                    stations.append(pick.station)
            if links >= min_link:
                f.write(event_text)
    f.close()
    return list(set(stations))

def write_correlations(event_list, wavbase, extract_len, pre_pick, shift_len,\
                       lowcut=1.0, highcut=10.0, max_sep=10, min_link=8):
    """
    Function to write a dt.cc file for hypoDD input - takes an input list of
    events and computes pick refienements by correlation

    :type event_list: List of tuple
    :param event_list: List of tuples of event_id (int) and sfile (String)
    :type wavbase: string
    :param wavbase: Path to the seisan wave directory that the wavefiles in the
                    S-files are stored
    :type extract_len: float
    :param extract_len: Length in seconds to extract around the pick
    :type pre_pick: float
    :param pre_pick: Time before the pick to start the correclation window
    :type shift_len: float
    :param shift_len: Time to allow pick to vary
    :type lowcut: float
    :param lowcut: Lowcut in Hz - default=1.0
    :type highcut: float
    :param highcut: Highcut in Hz - deafult=10.0
    :type max_sep: float
    :param max_sep: Maximum seperation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired
    """
    from obspy.signal.cross_correlation import xcorrPickCorrection
    from obspy import read
    f=open('dt.cc','w')
    for i in xrange(len(event_list)):
        master_sfile=event_list[i][1]
        master_event_id=event_list[i][0]
        master_picks=Sfile_util.readpicks(master_sfile)
        master_ori_time=Sfile_util.readheader(master_sfile).time
        master_location=(Sfile_util.readheader(master_sfile).latitude,\
                         Sfile_util.readheader(master_sfile).longitude,\
                         Sfile_util.readheader(master_sfile).depth)
        master_wavefiles=Sfile_util.readwavename(master_sfile)
        masterstream=read(wavbase+'/*/*/'+master_wavefiles[0])
        if len(master_wavefiles)>1:
            for wavefile in master_wavefiles:
                masterstream+=read(wavbase+'/*/*/'+wavefile)
        for j in xrange(i+1,len(event_list)):
            # Use this tactic to only output unique event pairings
            slave_sfile=event_list[j][1]
            slave_event_id=event_list[j][0]
            slave_wavefiles=Sfile_util.readwavename(slave_sfile)
            slavestream=read(wavbase+'/*/*/'+slave_wavefiles[0])
            if len(slave_wavefiles)>1:
                for wavefile in slave_wavefiles:
                    slavestream+=read(wavbase+'/*/*/'+wavefile)
            # Write out the header line
            event_text='#'+str(master_event_id).rjust(10)+\
                    str(slave_event_id).rjust(10)+' 0.0   \n'
            slave_picks=Sfile_util.readpicks(slave_sfile)
            slave_ori_time=Sfile_util.readheader(slave_sfile).time
            slave_location=(Sfile_util.readheader(slave_sfile).latitude,\
                         Sfile_util.readheader(slave_sfile).longitude,\
                         Sfile_util.readheader(slave_sfile).depth)
            if _separation(master_location, slave_location) > max_sep:
                break
            links=0
            for pick in master_picks:
                if pick.phase not in ['P','S']:
                    continue # Only use P and S picks, not amplitude or 'other'
                # Find station, phase pairs
                slave_matches=[p for p in slave_picks if p.station==pick.station\
                               and p.phase==pick.phase]
                print pick.station+'.'+pick.channel
                mastertr=masterstream.select(station=pick.station, \
                                             channel='*'+pick.channel[-1])[0]
                # Loop through the matches
                for slave_pick in slave_matches:
                    slavetr=slavestream.select(station=slave_pick.station,\
                                               channel='*'+slave_pick.channel[-1])[0]
                    # Correct the picks
                    try:
                        correction, cc = xcorrPickCorrection(pick.time, mastertr,\
                                                             slave_pick.time,\
                                                             slavetr,pre_pick,\
                                                             extract_len-pre_pick, shift_len,\
                                                             filter="bandpass",\
                                                             filter_options={'freqmin':lowcut,
                                                                             'freqmax':highcut},plot=False)
                        # Get the differntial travel time using the corrected time.

                        dt=(pick.time-master_ori_time)-\
                                (slave_pick.time+correction-slave_ori_time)
                        event_text+=pick.station.rjust(4)+\
                                _cc_round(dt,3).rjust(11)+\
                                _cc_round(cc*cc,3).rjust(8)+\
                                ' '+pick.phase+'\n'
                        links+=1
                    except:
                        continue
            if links >= min_link:
                f.write(event_text)
    f.close()
