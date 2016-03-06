"""
This module contains functions to convert seisan catalogue to files to hypoDD
input files.

These functions will generate both a catalogue (dt.ct) file, event
file (event.dat), station information file (station.dat), and a correlation
oiutput file correlated every event in the catalogue with every other event to
optimize the picks (dt.cc).

The correlation routine relies on obspy's xcorr_pick_correction function from
the obspy.signal.cross_correlation module.  This function optimizes picks to
better than sample accuracy by interpolating the correlation function and
finding the maximum of this rather than the true maximum correlation value.
The output from this function is stored in the dt.cc file.

Information for the station.dat file is read from SEISAN's STATION0.HYP file

Earthquake picks and locations are taken from the catalogued s-files - these
must be pre-located before entering this routine as origin times and hypocentre
locations are needed for event.dat files.

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
from eqcorrscan.utils import sfile_util
import os


def _cc_round(num, dp):
    """
    Convenience function to take a float and round it to dp padding with zeros
    to return a string

    :type num: float
    :param num: Number to round
    :type dp: int
    :param dp: Number of decimal places to round to.

    :returns: str
    """
    num = round(num, dp)
    num = '{0:.{1}f}'.format(num, dp)
    return num


def _av_weight(W1, W2):
    """
    Function to convert from two seisan weights (0-4) to one hypoDD \
    weight(0-1).

    :type W1: str
    :param W1: Seisan input weight (0-4)
    :type W2: str
    :param W2: Seisan input weight (0-4)

    :returns: str
    """
    if W1 == ' ':
        W1 = 1
    elif W1 == '9':
        W1 = 0
    else:
        W1 = 1 - int(W1) / 4.0

    if W2 == ' ':
        W2 = 1
    elif W2 == '9':
        W2 = 0
    else:
        W2 = 1 - int(W2) / 4.0

    W = (W1 + W2) / 2
    return _cc_round(W, 4)


def readSTATION0(path, stations):
    """
    Function to read the STATION0.HYP file on the path given.  Outputs written
    in station.dat file.

    :type path: str
    :param path: Path to the STATION0.HYP file
    :type station: list
    :param station: Stations to look for

    :returns: List of tuples of station, lat, long, elevation
    """
    stalist = []
    f = open(path + '/STATION0.HYP', 'r')
    for line in f:
        if line[1:6].strip() in stations:
            station = line[1:6].strip()
            lat = line[6:14]  # Format is either ddmm.mmS/N or ddmm(.)mmmS/N
            if lat[-1] == 'S':
                NS = -1
            else:
                NS = 1
            if lat[4] == '.':
                lat = (int(lat[0:2]) + float(lat[2:-1]) / 60) * NS
            else:
                lat = (int(lat[0:2])+float(lat[2:4] + '.' + lat[4:-1]) /
                       60) * NS
            lon = line[14:23]
            if lon[-1] == 'S':
                EW = -1
            else:
                EW = 1
            if lon[5] == '.':
                lon = (int(lon[0:3]) + float(lon[3:-1]) / 60) * EW
            else:
                lon = (int(lon[0:3]) + float(lon[3:5]+'.'+lon[5:-1]) /
                       60) * EW
            elev = float(line[23:-1].strip())
            # Note, negative altitude can be indicated in 1st column
            if line[0] == '-':
                elev *= -1
            stalist.append((station, lat, lon, elev))
    f.close()
    f = open('station.dat', 'w')
    for sta in stalist:
        line = ''.join([sta[0].ljust(5), _cc_round(sta[1], 4).ljust(10),
                        _cc_round(sta[2], 4).ljust(10),
                        _cc_round(sta[3]/1000, 4).rjust(7), '\n'])
        f.write(line)
    f.close()
    return stalist


def sfiles_to_event(sfile_list):
    """
    Function to write out an event.dat file of the events

    :type sfile_list: list
    :param sfile_list: List of s-files to sort and put into the database

    :returns: List of tuples of event ID (int) and Sfile name
    """
    from obspy.core.event import Catalog
    event_list = []
    sort_list = [(sfile_util.readheader(sfile).origins[0].time, sfile)
                 for sfile in sfile_list]
    sort_list.sort(key=lambda tup: tup[0])
    sfile_list = [sfile[1] for sfile in sort_list]
    catalog = Catalog()
    for i, sfile in enumerate(sfile_list):
        event_list.append((i, sfile))
        catalog.append(sfile_util.readheader(sfile))
    # Hand off to sister function
    write_event(catalog)
    return event_list


def write_event(catalog):
    """
    Function to write obspy.core.Catalog to a hypoDD format event.dat file.

    :type catalog: osbpy.core.Catalog
    :param catalog: A catalog of obspy events
    """
    f = open('event.dat', 'w')
    for i, event in enumerate(catalog):
        evinfo = event.origins[0]
        Mag_1 = event.magnitudes[0].mag or ' '
        if event.origins[0].time_errors:
            t_RMS = event.origins[0].time_errors.Time_Residual_RMS or ' '
        else:
            t_RMS = ' '
        f.write(str(evinfo.time.year)+str(evinfo.time.month).zfill(2) +
                str(evinfo.time.day).zfill(2)+'  ' +
                str(evinfo.time.hour).rjust(2) +
                str(evinfo.time.minute).zfill(2) +
                str(evinfo.time.second).zfill(2) +
                str(evinfo.time.microsecond)[0:2].zfill(2)+'  ' +
                str(evinfo.latitude).ljust(8, '0')+'   ' +
                str(evinfo.longitude).ljust(8, '0')+'  ' +
                str(evinfo.depth / 1000).rjust(7).ljust(9, '0')+'   ' +
                str(Mag_1)+'    0.00    0.00   ' +
                str(t_RMS).ljust(4, '0') +
                str(i).rjust(11)+'\n')
    f.close()
    return


def write_catalog(event_list, max_sep=1, min_link=8):
    """
    Function to write the dt.ct file needed by hypoDD - takes input event list
    from write_event as a list of tuples of event id and sfile.  It will read
    the pick information from the seisan formated s-file using the sfile_util
    utilities.

    :type event_list: list of tuple
    :param event_list: List of tuples of event_id (int) and sfile (String)
    :type max_sep: float
    :param max_sep: Maximum seperation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired

    :returns: list stations

    .. note:: Currently we have not implemented a method for taking \
        unassociated event objects and wavefiles.  As such if you have events \
        with associated wavefiles you are advised to generate Sfiles for each \
        event using the sfile_util module prior to this step.
    """
    from eqcorrscan.utils.mag_calc import dist_calc
    f = open('dt.ct', 'w')
    f2 = open('dt.ct2', 'w')
    fphase = open('phase.dat', 'w')
    stations = []
    evcount = 0
    for i, master in enumerate(event_list):
        master_sfile = master[1]
        master_event_id = master[0]
        master_event = sfile_util.readpicks(master_sfile)
        master_ori_time = master_event.origins[0].time
        master_location = (master_event.origins[0].latitude,
                           master_event.origins[0].longitude,
                           master_event.origins[0].depth / 1000)
        if len(master_event.magnitudes) > 0:
            master_magnitude = master_event.magnitudes[0].mag or ' '
        else:
            master_magnitude = ' '
        header = '# '+master_ori_time.strftime('%Y  %m  %d  %H  %M  %S.%f') +\
            ' '+str(master_location[0]).ljust(8)+' ' +\
            str(master_location[1]).ljust(8)+' ' +\
            str(master_location[2]).ljust(4)+' ' +\
            str(master_magnitude).ljust(4)+' 0.0 0.0 0.0' +\
            str(master_event_id).rjust(4)
        fphase.write(header+'\n')
        for pick in master_event.picks:
            if pick.phase_hint[0].upper() in ['P', 'S']:
                weight = [arrival.time_weight
                          for arrival in master_event.origins[0].arrivals
                          if arrival.pick_id == pick.resource_id][0]
                # Convert seisan weight to hypoDD 0-1 weights
                if weight == 0:
                    weight = 1.0
                elif weight == 9:
                    weight = 0.0
                else:
                    weight = 1 - weight / 4.0
                fphase.write(pick.waveform_id.station_code+'  ' +
                             _cc_round(pick.time -
                                       master_ori_time, 3).rjust(6) +
                             '   '+str(weight).ljust(5)+pick.phase_hint+'\n')
        for j in range(i+1, len(event_list)):
            # Use this tactic to only output unique event pairings
            slave_sfile = event_list[j][1]
            slave_event_id = event_list[j][0]
            # Write out the header line
            event_text = '#'+str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10)+'\n'
            event_text2 = '#'+str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10)+'\n'
            slave_event = sfile_util.readpicks(slave_sfile)
            slave_ori_time = slave_event.origins[0].time
            slave_location = (slave_event.origins[0].latitude,
                              slave_event.origins[0].longitude,
                              slave_event.origins[0].depth / 1000)
            if dist_calc(master_location, slave_location) > max_sep:
                continue
            links = 0  # Count the number of linkages
            for pick in master_event.picks:
                if pick.phase_hint[0].upper() not in ['P', 'S']:
                    continue
                    # Only use P and S picks, not amplitude or 'other'
                # Added by Carolin
                slave_matches = [p for p in slave_event.picks
                                 if p.phase_hint == pick.phase_hint
                                 and p.waveform_id.station_code.upper() ==
                                 pick.waveform_id.station_code.upper()]
                # Loop through the matches
                for slave_pick in slave_matches:
                    links += 1
                    master_weight = [arrival.time_weight
                                     for arrival in master_event.
                                     origins[0].arrivals
                                     if arrival.pick_id == pick.resource_id][0]
                    slave_weight = [arrival.time_weight
                                    for arrival in slave_event.
                                    origins[0].arrivals
                                    if arrival.pick_id ==
                                    slave_pick.resource_id][0]
                    event_text += pick.waveform_id.station_code.ljust(5) +\
                        _cc_round(pick.time-master_ori_time, 3).rjust(11) +\
                        _cc_round(slave_pick.time-slave_ori_time, 3).rjust(8) +\
                        _av_weight(master_weight, slave_weight).rjust(7) +\
                        ' '+pick.phase_hint+'\n'
                    # Added by Carolin
                    event_text2 += pick.waveform_id.station_code.ljust(5) +\
                        _cc_round(pick.time-master_ori_time, 3).rjust(11) +\
                        _cc_round(slave_pick.time-slave_ori_time, 3).rjust(8) +\
                        _av_weight(master_weight, slave_weight).rjust(7) +\
                        ' '+pick.phase_hint+'\n'
                    stations.append(pick.waveform_id.station_code)
            if links >= min_link:
                f.write(event_text)
                f2.write(event_text2)
                evcount += 1
    print('You have '+str(evcount)+' links')
    # f.write('\n')
    f.close()
    f2.close()
    fphase.close()
    return list(set(stations))


def write_correlations(event_list, wavbase, extract_len, pre_pick, shift_len,
                       lowcut=1.0, highcut=10.0, max_sep=4, min_link=8,
                       coh_thresh=0.0, coherence_weight=True, plotvar=False):
    """
    Function to write a dt.cc file for hypoDD input - takes an input list of
    events and computes pick refienements by correlation.

    :type event_list: list of tuple
    :param event_list: List of tuples of event_id (int) and sfile (String)
    :type wavbase: str
    :param wavbase: Path to the seisan wave directory that the wavefiles in the
                    S-files are stored
    :type extract_len: float
    :param extract_len: Length in seconds to extract around the pick
    :type pre_pick: float
    :param pre_pick: Time before the pick to start the correlation window
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
    :type coherence_weight: bool
    :param coherence_weight: Use coherence to weight the dt.cc file, or the \
        raw cross-correlation value, defaults to false which uses the cross-\
        correlation value.
    :type plotvar: bool
    :param plotvar: To show the pick-correction plots, defualts to False.

    .. warning:: This is not a fast routine!

    .. warning:: In contrast to seisan's \
        corr routine, but in accordance with the hypoDD manual, this outputs \
        corrected differential time.

    .. note:: Currently we have not implemented a method for taking \
        unassociated event objects and wavefiles.  As such if you have events \
        with associated wavefiles you are advised to generate Sfiles for each \
        event using the sfile_util module prior to this step.
    """
    import obspy
    if int(obspy.__version__.split('.')[0]) > 0:
        from obspy.signal.cross_correlation import xcorr_pick_correction
    else:
        from obspy.signal.cross_correlation import xcorrPickCorrection \
            as xcorr_pick_correction
    import matplotlib.pyplot as plt
    from obspy import read
    from eqcorrscan.utils.mag_calc import dist_calc
    import glob
    import warnings

    corr_list = []
    f = open('dt.cc', 'w')
    f2 = open('dt.cc2', 'w')
    for i, master in enumerate(event_list):
        master_sfile = master[1]
        master_event_id = master[0]
        master_picks = sfile_util.readpicks(master_sfile).picks
        master_event = sfile_util.readheader(master_sfile)
        master_ori_time = master_event.origins[0].time
        master_location = (master_event.origins[0].latitude,
                           master_event.origins[0].longitude,
                           master_event.origins[0].depth)
        master_wavefiles = sfile_util.readwavename(master_sfile)
        masterpath = glob.glob(wavbase + os.sep + master_wavefiles[0])
        if masterpath:
            masterstream = read(masterpath[0])
        if len(master_wavefiles) > 1:
            for wavefile in master_wavefiles:
                try:
                    masterstream += read(os.join(wavbase, wavefile))
                except:
                    continue
                    raise IOError("Couldn't find wavefile")
        for j in range(i+1, len(event_list)):
            # Use this tactic to only output unique event pairings
            slave_sfile = event_list[j][1]
            slave_event_id = event_list[j][0]
            slave_wavefiles = sfile_util.readwavename(slave_sfile)
            try:
                # slavestream=read(wavbase+'/*/*/'+slave_wavefiles[0])
                slavestream = read(wavbase + os.sep + slave_wavefiles[0])
            except:
                # print(slavestream)
                raise IOError('No wavefile found: '+slave_wavefiles[0]+' ' +
                              slave_sfile)
            if len(slave_wavefiles) > 1:
                for wavefile in slave_wavefiles:
                    # slavestream+=read(wavbase+'/*/*/'+wavefile)
                    try:
                        slavestream += read(wavbase+'/'+wavefile)
                    except:
                        continue
            # Write out the header line
            event_text = '#'+str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10)+' 0.0   \n'
            event_text2 = '#'+str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10)+' 0.0   \n'
            slave_picks = sfile_util.readpicks(slave_sfile).picks
            slave_event = sfile_util.readheader(slave_sfile)
            slave_ori_time = slave_event.origins[0].time
            slave_location = (slave_event.origins[0].latitude,
                              slave_event.origins[0].longitude,
                              slave_event.origins[0].depth)
            if dist_calc(master_location, slave_location) > max_sep:
                continue
            links = 0
            phases = 0
            for pick in master_picks:
                if pick.phase_hint[0].upper() not in ['P', 'S']:
                    continue
                    # Only use P and S picks, not amplitude or 'other'
                # Find station, phase pairs
                # Added by Carolin
                slave_matches = [p for p in slave_picks
                                 if p.phase_hint == pick.phase_hint
                                 and p.waveform_id.station_code ==
                                 pick.waveform_id.station_code]

                if masterstream.select(station=pick.waveform_id.station_code,
                                       channel='*' +
                                       pick.waveform_id.channel_code[-1]):
                    mastertr = masterstream.\
                        select(station=pick.waveform_id.station_code,
                               channel='*' +
                               pick.waveform_id.channel_code[-1])[0]
                else:
                    print('No waveform data for ' +
                          pick.waveform_id.station_code + '.' +
                          pick.waveform_id.channel_code)
                    print(pick.waveform_id.station_code +
                          '.' + pick.waveform_id.channel_code +
                          ' ' + slave_sfile+' ' + master_sfile)
                    break
                # Loop through the matches
                for slave_pick in slave_matches:
                    if slavestream.select(station=slave_pick.waveform_id.
                                          station_code,
                                          channel='*'+slave_pick.waveform_id.
                                          channel_code[-1]):
                        slavetr = slavestream.\
                            select(station=slave_pick.waveform_id.station_code,
                                   channel='*'+slave_pick.waveform_id.
                                   channel_code[-1])[0]
                    else:
                        print('No slave data for ' +
                              slave_pick.waveform_id.station_code + '.' +
                              slave_pick.waveform_id.channel_code)
                        print(pick.waveform_id.station_code +
                              '.' + pick.waveform_id.channel_code +
                              ' ' + slave_sfile + ' ' + master_sfile)
                        break
                    # Correct the picks
                    try:
                        correction, cc =\
                            xcorr_pick_correction(pick.time, mastertr,
                                                  slave_pick.time,
                                                  slavetr, pre_pick,
                                                  extract_len - pre_pick,
                                                  shift_len, filter="bandpass",
                                                  filter_options={'freqmin':
                                                                  lowcut,
                                                                  'freqmax':
                                                                  highcut},
                                                  plot=plotvar)
                        # Get the differntial travel time using the
                        # corrected time.
                        # Check that the correction is within the allowed shift
                        # This can occur in the obspy routine when the
                        # correlation function is increasing at the end of the
                        # window.
                        if abs(correction) > shift_len:
                            warnings.warn('Shift correction too large, ' +
                                          'will not use')
                            continue
                        correction = (pick.time - master_ori_time) -\
                            (slave_pick.time + correction - slave_ori_time)
                        links += 1
                        if cc * cc >= coh_thresh:
                            if coherence_weight:
                                weight = cc * cc
                            else:
                                weight = cc
                            phases += 1
                            # added by Caro
                            event_text += pick.waveform_id.station_code.\
                                ljust(5) + _cc_round(correction, 3).\
                                rjust(11) + _cc_round(weight, 3).rjust(8) +\
                                ' '+pick.phase_hint+'\n'
                            event_text2 += pick.waveform_id.station_code\
                                .ljust(5).upper() +\
                                _cc_round(correction, 3).rjust(11) +\
                                _cc_round(weight, 3).rjust(8) +\
                                ' '+pick.phase_hint+'\n'

                            # links+=1
                        corr_list.append(cc*cc)
                    except:
                        # Should warn here
                        msg = "Couldn't compute correlation correction"
                        warnings.warn(msg)
                        continue
            if links >= min_link and phases > 0:
                f.write(event_text)
                f2.write(event_text2)
    if plotvar:
        plt.hist(corr_list, 150)
        plt.show()
    # f.write('\n')
    f.close()
    f2.close()
    return
