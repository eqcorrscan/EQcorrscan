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

.. todo:
    Change from forced seisan usage, to using obspy events.  This happens
    internally already, but the insides should be extracted and thin wrappers
    written for seisan.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import warnings
import matplotlib.pyplot as plt

from obspy.core.event import Catalog
from obspy import read

from eqcorrscan.utils import sfile_util
from eqcorrscan.utils.mag_calc import dist_calc


def _cc_round(num, dp):
    """
    Convenience function to take a float and round it to dp padding with zeros
    to return a string

    :type num: float
    :param num: Number to round
    :type dp: int
    :param dp: Number of decimal places to round to.

    :returns: str

    >>> print(_cc_round(0.25364, 2))
    0.25
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

    .. rubric:: Example
    >>> print(_av_weight(1, 4))
    0.3750
    >>> print(_av_weight(0, 0))
    1.0000
    >>> print(_av_weight(' ', ' '))
    1.0000
    >>> print(_av_weight(-9, 0))
    0.5000
    >>> print(_av_weight(1, -9))
    0.3750
    """
    import warnings

    if str(W1) in [' ', '']:
        W1 = 1
    elif str(W1) in ['-9', '9', '9.0', '-9.0']:
        W1 = 0
    else:
        W1 = 1 - (int(W1) / 4.0)
    if W1 < 0:
        warnings.warn('Negative weight found, setting to zero')
        W2 = 0

    if str(W2) in [' ', '']:
        W2 = 1
    elif str(W2) in ['-9', '9', '9.0', '-9.0']:
        W2 = 0
    else:
        W2 = 1 - (int(W2) / 4.0)
    if W2 < 0:
        warnings.warn('Negative weight found, setting to zero')
        W2 = 0

    W = (W1 + W2) / 2
    if W < 0:
        print('Weight 1: ' + str(W1))
        print('Weight 2: ' + str(W2))
        print('Final weight: ' + str(W))
        raise IOError('Negative average weight calculated, setting to zero')
    return _cc_round(W, 4)


def readSTATION0(path, stations):
    """
    Read a Seisan STATION0.HYP file on the path given.

    Outputs the information, and writes to station.dat file.

    :type path: str
    :param path: Path to the STATION0.HYP file
    :type stations: list
    :param stations: Stations to look for

    :returns: List of tuples of station, lat, long, elevation
    :rtype: list

    >>> readSTATION0('eqcorrscan/tests/test_data', ['WHFS', 'WHAT2'])
    [('WHFS', -43.261, 170.359, 60.0), ('WHAT2', -43.2793, \
170.36038333333335, 95.0)]
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
                lat = (int(lat[0:2]) + float(lat[2:4] + '.' + lat[4:-1]) /
                       60) * NS
            lon = line[14:23]
            if lon[-1] == 'S':
                EW = -1
            else:
                EW = 1
            if lon[5] == '.':
                lon = (int(lon[0:3]) + float(lon[3:-1]) / 60) * EW
            else:
                lon = (int(lon[0:3]) + float(lon[3:5] + '.' + lon[5:-1]) /
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
                        _cc_round(sta[3] / 1000, 4).rjust(7), '\n'])
        f.write(line)
    f.close()
    return stalist


def sfiles_to_event(sfile_list):
    """
    Write an event.dat file from a list of Seisan events

    :type sfile_list: list
    :param sfile_list: List of s-files to sort and put into the database

    :returns: List of tuples of event ID (int) and Sfile name
    """
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
    Write obspy.core.event.Catalog to a hypoDD format event.dat file.

    :type catalog: obspy.core.event.Catalog
    :param catalog: A catalog of obspy events.
    """
    f = open('event.dat', 'w')
    for i, event in enumerate(catalog):
        evinfo = event.origins[0]
        Mag_1 = event.magnitudes[0].mag or ' '
        if event.origins[0].time_errors.Time_Residual_RMS:
            t_RMS = event.origins[0].time_errors.Time_Residual_RMS
        else:
            print('No time residual in header')
            t_RMS = 0.0
        f.write(str(evinfo.time.year) + str(evinfo.time.month).zfill(2) +
                str(evinfo.time.day).zfill(2) + '  ' +
                str(evinfo.time.hour).rjust(2) +
                str(evinfo.time.minute).zfill(2) +
                str(evinfo.time.second).zfill(2) +
                str(evinfo.time.microsecond)[0:2].zfill(2) + '  ' +
                str(evinfo.latitude).ljust(8, str('0')) + '   ' +
                str(evinfo.longitude).ljust(8, str('0')) + '  ' +
                str(evinfo.depth / 1000).rjust(7).ljust(9, str('0')) + '   ' +
                str(Mag_1) + '    0.00    0.00   ' +
                str(t_RMS).ljust(4, str('0')) +
                str(i).rjust(11) + '\n')
    f.close()
    return


def write_catalog(event_list, max_sep=8, min_link=8):
    """
    Generate a dt.ct for hypoDD for a series of events.

    Takes input event list from
    :func:`eqcorrscan.utils.catalog_to_dd.write_event` as a list of tuples of
    event id and sfile.  It will read the pick information from the seisan
    formated s-file using the sfile_util utilities.

    :type event_list: list
    :param event_list: List of tuples of event_id (int) and sfile (String)
    :type max_sep: float
    :param max_sep: Maximum separation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired

    :returns: list of stations that have been used in this catalog

    .. note::
        We have not yet implemented a method for taking unassociated event
        objects and wavefiles.  As such if you have events with associated
        wavefiles you are advised to generate Sfiles for each event using
        the :mod:`eqcorrscan.utils.sfile_util` module prior to this step.
    """
    # Cope with possibly being passed a zip in python 3.x
    event_list = list(event_list)
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
        header = '# ' + \
            master_ori_time.strftime('%Y  %m  %d  %H  %M  %S.%f') +\
            ' ' + str(master_location[0]).ljust(8) + ' ' +\
            str(master_location[1]).ljust(8) + ' ' +\
            str(master_location[2]).ljust(4) + ' ' +\
            str(master_magnitude).ljust(4) + ' 0.0 0.0 0.0' +\
            str(master_event_id).rjust(4)
        fphase.write(header + '\n')
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
                fphase.write(pick.waveform_id.station_code + '  ' +
                             _cc_round(pick.time -
                                       master_ori_time, 3).rjust(6) +
                             '   ' + str(weight).ljust(5) +
                             pick.phase_hint + '\n')
        for j in range(i + 1, len(event_list)):
            # Use this tactic to only output unique event pairings
            slave_sfile = event_list[j][1]
            slave_event_id = event_list[j][0]
            # Write out the header line
            event_text = '#' + str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10) + '\n'
            event_text2 = '#' + str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10) + '\n'
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
                                 if p.phase_hint == pick.phase_hint and
                                 p.waveform_id.station_code.upper() ==
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
                    master_weight = str(int(master_weight))
                    slave_weight = str(int(slave_weight))
                    event_text += pick.waveform_id.station_code.ljust(5) +\
                        _cc_round(pick.time - master_ori_time, 3).rjust(11) +\
                        _cc_round(slave_pick.time -
                                  slave_ori_time, 3).rjust(8) +\
                        _av_weight(master_weight, slave_weight).rjust(7) +\
                        ' ' + pick.phase_hint + '\n'
                    # Added by Carolin
                    event_text2 += pick.waveform_id.station_code.ljust(5) +\
                        _cc_round(pick.time - master_ori_time, 3).rjust(11) +\
                        _cc_round(slave_pick.time -
                                  slave_ori_time, 3).rjust(8) +\
                        _av_weight(master_weight, slave_weight).rjust(7) +\
                        ' ' + pick.phase_hint + '\n'
                    stations.append(pick.waveform_id.station_code)
            if links >= min_link:
                f.write(event_text)
                f2.write(event_text2)
                evcount += 1
    print('You have ' + str(evcount) + ' links')
    # f.write('\n')
    f.close()
    f2.close()
    fphase.close()
    return list(set(stations))


def write_correlations(event_list, wavbase, extract_len, pre_pick, shift_len,
                       lowcut=1.0, highcut=10.0, max_sep=8, min_link=8,
                       cc_thresh=0.0, plotvar=False, debug=0):
    """
    Write a dt.cc file for hypoDD input for a given list of events.

    Takes an input list of events and computes pick refinements by correlation.
    Outputs two files, dt.cc and dt.cc2, each provides a different weight,
    dt.cc uses weights of the cross-correlation, and dt.cc2 provides weights
    as the square of the cross-correlation.

    :type event_list: list
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
    :param highcut: Highcut in Hz - default=10.0
    :type max_sep: float
    :param max_sep: Maximum separation between event pairs in km
    :type min_link: int
    :param min_link: Minimum links for an event to be paired
    :type cc_thresh: float
    :param cc_thresh: Threshold to include cross-correlation results.
    :type plotvar: bool
    :param plotvar: To show the pick-correction plots, defualts to False.
    :type debug: int
    :param debug: Variable debug levels from 0-5, higher=more output.

    .. warning:: This is not a fast routine!

    .. warning::
        In contrast to seisan's corr routine, but in accordance with the
        hypoDD manual, this outputs corrected differential time.

    .. note::
        Currently we have not implemented a method for taking
        unassociated event objects and wavefiles.  As such if you have events \
        with associated wavefiles you are advised to generate Sfiles for each \
        event using the sfile_util module prior to this step.

    .. note::
        There is no provision to taper waveforms within these functions, if you
        desire this functionality, you should apply the taper before calling
        this.  Note the :func:`obspy.Trace.taper` functions.
    """
    import obspy
    if int(obspy.__version__.split('.')[0]) > 0:
        from obspy.signal.cross_correlation import xcorr_pick_correction
    else:
        from obspy.signal.cross_correlation import xcorrPickCorrection \
            as xcorr_pick_correction
    warnings.filterwarnings(action="ignore",
                            message="Maximum of cross correlation " +
                                    "lower than 0.8: *")
    corr_list = []
    f = open('dt.cc', 'w')
    f2 = open('dt.cc2', 'w')
    k_events = len(list(event_list))
    for i, master in enumerate(event_list):
        master_sfile = master[1]
        if debug > 1:
            print('Computing correlations for master: %s' % master_sfile)
        master_event_id = master[0]
        master_picks = sfile_util.readpicks(master_sfile).picks
        master_event = sfile_util.readheader(master_sfile)
        master_ori_time = master_event.origins[0].time
        master_location = (master_event.origins[0].latitude,
                           master_event.origins[0].longitude,
                           master_event.origins[0].depth / 1000.0)
        master_wavefiles = sfile_util.readwavename(master_sfile)
        masterpath = glob.glob(wavbase + os.sep + master_wavefiles[0])
        if masterpath:
            masterstream = read(masterpath[0])
        if len(master_wavefiles) > 1:
            for wavefile in master_wavefiles:
                try:
                    masterstream += read(os.join(wavbase, wavefile))
                except:
                    raise IOError("Couldn't find wavefile")
                    continue
        for j in range(i + 1, k_events):
            # Use this tactic to only output unique event pairings
            slave_sfile = event_list[j][1]
            if debug > 2:
                print('Comparing to event: %s' % slave_sfile)
            slave_event_id = event_list[j][0]
            slave_wavefiles = sfile_util.readwavename(slave_sfile)
            try:
                slavestream = read(wavbase + os.sep + slave_wavefiles[0])
            except:
                raise IOError('No wavefile found: ' + slave_wavefiles[0] +
                              ' ' + slave_sfile)
            if len(slave_wavefiles) > 1:
                for wavefile in slave_wavefiles:
                    try:
                        slavestream += read(wavbase + os.sep + wavefile)
                    except IOError:
                        print('No waveform found: %s' %
                              (wavbase + os.sep + wavefile))
                        continue
            # Write out the header line
            event_text = '#' + str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10) + ' 0.0   \n'
            event_text2 = '#' + str(master_event_id).rjust(10) +\
                str(slave_event_id).rjust(10) + ' 0.0   \n'
            slave_picks = sfile_util.readpicks(slave_sfile).picks
            slave_event = sfile_util.readheader(slave_sfile)
            slave_ori_time = slave_event.origins[0].time
            slave_location = (slave_event.origins[0].latitude,
                              slave_event.origins[0].longitude,
                              slave_event.origins[0].depth / 1000.0)
            if dist_calc(master_location, slave_location) > max_sep:
                if debug > 0:
                    print('Seperation exceeds max_sep: %s' %
                          (dist_calc(master_location, slave_location)))
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
                                 if p.phase_hint == pick.phase_hint and
                                 p.waveform_id.station_code ==
                                 pick.waveform_id.station_code]

                if masterstream.select(station=pick.waveform_id.station_code,
                                       channel='*' +
                                       pick.waveform_id.channel_code[-1]):
                    mastertr = masterstream.\
                        select(station=pick.waveform_id.station_code,
                               channel='*' +
                               pick.waveform_id.channel_code[-1])[0]
                elif debug > 1:
                    print('No waveform data for ' +
                          pick.waveform_id.station_code + '.' +
                          pick.waveform_id.channel_code)
                    print(pick.waveform_id.station_code +
                          '.' + pick.waveform_id.channel_code +
                          ' ' + slave_sfile + ' ' + master_sfile)
                    break
                # Loop through the matches
                for slave_pick in slave_matches:
                    if slavestream.select(station=slave_pick.waveform_id.
                                          station_code,
                                          channel='*' + slave_pick.waveform_id.
                                          channel_code[-1]):
                        slavetr = slavestream.\
                            select(station=slave_pick.waveform_id.station_code,
                                   channel='*' + slave_pick.waveform_id.
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
                        # Get the differential travel time using the
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
                        if cc >= cc_thresh:
                            weight = cc
                            phases += 1
                            # added by Caro
                            event_text += pick.waveform_id.station_code.\
                                ljust(5) + _cc_round(correction, 3).\
                                rjust(11) + _cc_round(weight, 3).rjust(8) +\
                                ' ' + pick.phase_hint + '\n'
                            event_text2 += pick.waveform_id.station_code\
                                .ljust(5) + _cc_round(correction, 3).\
                                rjust(11) +\
                                _cc_round(weight * weight, 3).rjust(8) +\
                                ' ' + pick.phase_hint + '\n'
                            if debug > 3:
                                print(event_text)
                        else:
                            print('cc too low: %s' % cc)
                        corr_list.append(cc * cc)
                    except:
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


def read_phase(ph_file):
    """
    Read hypoDD phase files into Obspy catalog class.

    :type ph_file: str
    :param ph_file: Phase file to read event info from.

    :returns: Catalog of events from file.
    :rtype: :class:`obspy.core.event.Catalog`

    >>> from obspy.core.event.catalog import Catalog
    >>> catalog = read_phase('eqcorrscan/tests/test_data/tunnel.phase')
    >>> isinstance(catalog, Catalog)
    True
    """
    from obspy.core.event import Catalog
    ph_catalog = Catalog()
    f = open(ph_file, 'r')
    # Topline of each event is marked by # in position 0
    for line in f:
        if line[0] == '#':
            if 'event_text' not in locals():
                event_text = {'header': line.rstrip(),
                              'picks': []}
            else:
                ph_catalog.append(_phase_to_event(event_text))
                event_text = {'header': line.rstrip(),
                              'picks': []}
        else:
            event_text['picks'].append(line.rstrip())
    ph_catalog.append(_phase_to_event(event_text))
    return ph_catalog


def _phase_to_event(event_text):
    """
    Function to convert the text for one event in hypoDD phase format to \
    event object.

    :type event_text: dict
    :param event_text: dict of two elements, header and picks, header is a \
        str, picks is a list of str.

    :returns: obspy.core.event.Event
    """
    from obspy.core.event import Event, Origin, Magnitude
    from obspy.core.event import Pick, WaveformStreamID, Arrival
    from obspy import UTCDateTime
    ph_event = Event()
    # Extract info from header line
    # YR, MO, DY, HR, MN, SC, LAT, LON, DEP, MAG, EH, EZ, RMS, ID
    header = event_text['header'].split()
    ph_event.origins.append(Origin())
    ph_event.origins[0].time =\
        UTCDateTime(year=int(header[1]), month=int(header[2]),
                    day=int(header[3]), hour=int(header[4]),
                    minute=int(header[5]), second=int(header[6].split('.')[0]),
                    microsecond=int(float(('0.' + header[6].split('.')[1])) *
                                    1000000))
    ph_event.origins[0].latitude = float(header[7])
    ph_event.origins[0].longitude = float(header[8])
    ph_event.origins[0].depth = float(header[9]) * 1000
    ph_event.origins[0].time_errors['Time_Residual_RMS'] = float(header[13])
    ph_event.magnitudes.append(Magnitude())
    ph_event.magnitudes[0].mag = float(header[10])
    ph_event.magnitudes[0].magnitude_type = 'M'
    # Extract arrival info from picks!
    for i, pick_line in enumerate(event_text['picks']):
        pick = pick_line.split()
        _waveform_id = WaveformStreamID(station_code=pick[0])
        pick_time = ph_event.origins[0].time + float(pick[1])
        ph_event.picks.append(Pick(waveform_id=_waveform_id,
                                   phase_hint=pick[3],
                                   time=pick_time))
        ph_event.origins[0].arrivals.append(Arrival(phase=ph_event.picks[i],
                                                    pick_id=ph_event.picks[i].
                                                    resource_id))
        ph_event.origins[0].arrivals[i].time_weight = float(pick[2])
    return ph_event


if __name__ == '__main__':
    import doctest
    doctest.testmod()
