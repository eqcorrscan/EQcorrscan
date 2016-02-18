#!/usr/bin/python
"""
Part of the EQcorrscan module to read nordic format s-files and write them
EQcorrscan is a python module designed to run match filter routines for
seismology, within it are routines for integration to seisan and obspy.
With obspy integration (which is necessary) all main waveform formats can be
read in and output.

Note that these functions do not provide full functionality between quakeML
and seisan s-files.  Currently (as of version 0.0.9) these only convert pick
times and phase information, along with amplitude information for local
magnitudes between seisan and quakeML.  Location information including
hypocentre, origin time and magnitudes are also handled.

A series of wrappers and conversions is included between the legacy PICK
and EVENTINFO classes, however these will be depreciated along with these
classes for version 0.1.0.  Users should transition to using obspy.core.event
classes as these have more support and functionality.

We have not implimented any handling of focal mechanism solutions between
the two formats.

Code written by Calum John Chamberlain and Chet Hopp both of
Victoria University of Wellington, 2015 & 2016.

Copyright 2015, 2016 the authors.

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


from obspy import UTCDateTime
import numpy as np


class PICK:
    """
    Pick information for seisan implimentation, note all fields can be left\
    blank to obtain a default pick: picks have a print function which will\
    print them as they would be seen in an S-file.

    Attributes:
        :type station: str
        :param station: Station name, less than five charectars required as\
         standard
        :type channel: str
        :param channel: Two or three charactar channel name, stored as two\
            charactars in S-file
        :type impulsivity: str
        :param impulsivity: either 'C' or 'D' for compressive and dilatational
        :type phase: str
        :param phase: Any allowable phase name in two characters
        :type weight: int
        :param weight: 0-4 with 0=100%, 4=0%, use weight=9 for unknown timing
        :type polarity: str
        :type time: obspy.UTCDateTime()
        :param time: Pick time as an obspy.UTCDateTime object
        :type coda: int
        :param coda: Length of coda in seconds
        :type amplitude: float
        :param amplitude: Amplitude (zero-peak), type is given in phase
        :type peri: float
        :param peri: Period of amplitude
        :type azimuth: float
        :param azimuth: Direction of approach in degrees
        :type velocity: float
        :param velocity: Phase velocity (km/s)
        :type AIN: int
        :param AIN: Angle of incidence.
        :type SNR: float
        :param SNR: Signal to noise ratio
        :type azimuthres: int
        :param azimuthres: Residual azimuth
        :type timeres: float
        :param timeres: Time residual in seconds
        :type finalweight: int
        :param finalweight: Final weight used in location
        :type distance: float
        :param distance: Source-reciever distance in km
        :type CAZ: int
        :param CAZ: Azimuth at source.

    .. rubric:: Note: Depreciated legacy function, use the obspy.core.event \
    classes. This will be removed in future releases.
    """
    pickcount = 0

    def __init__(self, station=' ', channel=' ', impulsivity=' ', phase=' ',
                 weight=999, polarity=' ', time=UTCDateTime(0),
                 coda=999, amplitude=float('NaN'),
                 peri=float('NaN'), azimuth=float('NaN'),
                 velocity=float('NaN'), AIN=999, SNR=float('NaN'),
                 azimuthres=999, timeres=float('NaN'),
                 finalweight=999, distance=float('NaN'),
                 CAZ=999):
        self.station = station
        self.channel = channel
        self.impulsivity = impulsivity
        self.phase = phase
        self.weight = weight
        self.polarity = polarity
        self.time = time
        self.coda = coda
        self.amplitude = amplitude
        self.peri = peri
        self.azimuth = azimuth
        self.velocity = velocity
        self.AIN = AIN
        self.SNR = SNR
        self.azimuthres = azimuthres
        self.timeres = timeres
        self.finalweight = finalweight
        self.distance = distance
        self.CAZ = CAZ
        self.pickcount += 1

    def __repr__(self):
        return "PICK()"

    def __str__(self):
        if self.distance >= 100.0:
            self.distance = _int_conv(self.distance)
        elif 10.0 < self.distance < 100.0:
            self.distance = round(self.distance, 1)
            round_len = 1
        elif self.distance < 10.0:
            self.distance = round(self.distance, 2)
            round_len = 2
        else:
            round_len = False
        if self.peri < 10.0:
            peri_round = 2
        elif self.peri >= 10.0:
            peri_round = 1
        else:
            peri_round = False
        if not self.AIN == '':
            if not np.isnan(self.AIN):
                dummy = int(self.AIN)
            else:
                dummy = self.AIN
        else:
            dummy = self.SNR
        print(_str_conv(self.weight).rjust(1))
        print_str = ' ' + self.station.ljust(5) +\
            self.channel[0]+self.channel[len(self.channel)-1] +\
            ' ' + self.impulsivity +\
            self.phase.ljust(4) +\
            _str_conv(self.weight).rjust(1) + ' ' +\
            self.polarity.rjust(1) + ' ' +\
            str(self.time.hour).rjust(2) +\
            str(self.time.minute).rjust(2) +\
            str(self.time.second).rjust(3) + '.' +\
            str(float(self.time.microsecond) /
                (10 ** 4)).split('.')[0].zfill(2) +\
            _str_conv(int(self.coda)).rjust(5)[0:5] +\
            _str_conv(round(self.amplitude, 1)).rjust(7)[0:7] +\
            _str_conv(self.peri, rounded=peri_round).rjust(5) +\
            _str_conv(self.azimuth).rjust(6) +\
            _str_conv(self.velocity).rjust(5) +\
            _str_conv(dummy).rjust(4) +\
            _str_conv(int(self.azimuthres)).rjust(3) +\
            _str_conv(self.timeres, rounded=2).rjust(5) +\
            _str_conv(int(self.finalweight)).rjust(2) +\
            _str_conv(self.distance, rounded=round_len).rjust(5) +\
            _str_conv(int(self.CAZ)).rjust(4)+' '
        return print_str

    def write(self, filename):
        """
        Public function to write the pick to a file

        :type filename: str
        :param filename: Path to file to write to - will append to file
        """
        import os
        import warnings
        if os.path.isfile(filename):
            open_as = 'a'
        else:
            warnings.warn('File does not exist, no header')
            open_as = 'w'

        with open(filename, open_as) as f:
            pickstr = self.__str__()
            f.write(pickstr+'\n')
        return


class EVENTINFO:
    """
    Header information for seisan events, again all fields can be left blank \
    for a default empty header.  The print function for header will print \
    important information, but not as seen in an S-file.

    For more information on parameters see the seisan manual.

    Attributes:
        :type time: obspy.UTCDateTime
        :param time: Event origin time
        :type loc_mod_ind: str
        :param loc_mod_ind:
        :type dist_ind: str
        :param dist_ind: Distance flag, usually 'L' for local, 'R' for \
            regional and 'D' for distant
        :type ev_id: str
        :param ev_id: Often blank, 'E' denotes explosion and fixes depth to 0km
        :type latitude: float
        :param latitude: Hypocentre latitude in decimal degrees
        :type longitude: float
        :param lognitude: Hypocentre longitude in decimal degrees
        :type depth: float
        :param depth: hypocentre depth in km
        :type depth_ind: str
        :param depth_ind:
        :type loc_ind: str
        :param loc_ind:
        :type agency: str
        :param agency: Reporting agency, three letters
        :type nsta: int
        :param nsta: Number of stations recording
        :type t_RMS: float
        :param t_RMS: Root-mean-squared time residual
        :type Mag_1: float
        :param Mag_1: first magnitude
        :type Mag_1_type: str
        :param Mag_1_type: Type of magnitude for Mag_1 ('L', 'C', 'W')
        :type Mag_1_agency: str
        :param Mag_1_agency: Reporting agency for Mag_1
        :type Mag_2: float
        :param Mag_2: second magnitude
        :type Mag_2_type: str
        :param Mag_2_type: Type of magnitude for Mag_2 ('L', 'C', 'W')
        :type Mag_2_agency: str
        :param Mag_2_agency: Reporting agency for Mag_2
        :type Mag_3: float
        :param Mag_3: third magnitude
        :type Mag_3_type: str
        :param Mag_3_type: Type of magnitude for Mag_3 ('L', 'C', 'W')
        :type Mag_3_agency: str
        :param Mag_3_agency: Reporting agency for Mag_3

    .. rubric:: Note: Depreciated legacy function, use the obspy.core.event \
    classes. This will be removed in future releases.
    """
    def __init__(self, time=UTCDateTime(0), loc_mod_ind=' ', dist_ind=' ',
                 ev_id=' ', latitude=float('NaN'), longitude=float('NaN'),
                 depth=float('NaN'), depth_ind=' ', loc_ind=' ', agency=' ',
                 nsta=0, t_RMS=float('NaN'), Mag_1=float('NaN'),
                 Mag_1_type=' ', Mag_1_agency=' ', Mag_2=float('NaN'),
                 Mag_2_type=' ', Mag_2_agency=' ', Mag_3=float('NaN'),
                 Mag_3_type=' ', Mag_3_agency=' '):
        self.time = time
        self.loc_mod_ind = loc_mod_ind
        self.dist_ind = dist_ind
        self.ev_id = ev_id
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.depth_ind = depth_ind
        self.loc_ind = loc_ind
        self.agency = agency
        self.nsta = nsta
        self.t_RMS = t_RMS
        self.Mag_1 = Mag_1
        self.Mag_1_type = Mag_1_type
        self.Mag_1_agency = Mag_1_agency
        self.Mag_2 = Mag_2
        self.Mag_2_type = Mag_2_type
        self.Mag_2_agency = Mag_2_agency
        self.Mag_3 = Mag_3
        self.Mag_3_type = Mag_3_type
        self.Mag_3_agency = Mag_3_agency

    def __repr__(self):
        return "HEADER()"

    def __str__(self):
        print_str = str(self.time) + ' ' + str(self.latitude) + ',' +\
            str(self.longitude) + ' ' + str(self.depth) + ' ' +\
            self.Mag_1_type + ':' + str(self.Mag_1) + self.Mag_2_type + ':' +\
            str(self.Mag_2) + self.Mag_3_type + ':' + str(self.Mag_3) + '  '\
            + self.agency
        return print_str


def _int_conv(string):
    """
    Convenience tool to convert from string to integer, if empty string \
    return a 999 rather than an error.
    """
    try:
        intstring = int(string)
    except:
        intstring = 999
    return intstring


def _float_conv(string):
    """
    Convenience tool to convert from string to float, if empty string return \
    NaN rather than an error
    """
    try:
        floatstring = float(string)
    except:
        # cjh Changed this to deal with QuakeML validation issues with 'NaN'
        floatstring = float('999')
    return floatstring


def _str_conv(number, rounded=False):
    """
    Convenience tool to convert a number, either float or into into a string, \
    if the int is 999, or the float is NaN, returns empty string.
    """
    if (type(number) == float and np.isnan(number)) or number == 999:
        string = ' '
    elif type(number) == str:
        return number
    elif not rounded:
        if number < 100000:
            string = str(number)
        else:
            exponant = int('{0:.2E}'.format(number).split('E+')[-1]) - 1
            divisor = 10 ** exponant
            string = '{0:.1f}'.format(number / divisor) + 'e' + str(exponant)
    elif rounded == 2:
        string = '{0:.2f}'.format(number)
    elif rounded == 1:
        string = '{0:.1f}'.format(number)
    return string


def _evmagtonor(mag_type):
    """
    Convenience tool to switch from obspy event magnitude types to seisan \
    syntax
    """
    if mag_type == 'ML':
        mag = 'L'
    elif mag_type == 'mB':
        mag = 'b'
    elif mag_type == 'Ms':
        mag = 's'
    elif mag_type == 'MS':
        mag = 'S'
    elif mag_type == 'MW':
        mag = 'W'
    elif mag_type == 'MbLg':
        mag = 'G'
    elif mag_type == 'Mc':
        mag = 'C'
    else:
        raise IOError(mag_type + ' is not convertable')
    return mag


def _nortoevmag(mag_type):
    """
    Convenience tool to switch from nordic type magnitude notation to obspy \
    event magnitudes.
    """
    if mag_type == 'L':
        mag = 'ML'
    elif mag_type == 'b':
        mag = 'mB'
    elif mag_type == 's':
        mag = 'Ms'
    elif mag_type == 'S':
        mag = 'MS'
    elif mag_type == 'W':
        mag = 'MW'
    elif mag_type == 'G':
        mag = 'MbLg'
    elif mag_type == 'C':
        mag = 'Mc'
    else:
        raise IOError(mag_type + ' is not convertable')
    return mag


def readheader(sfile):
    """
    Function to read the header information from a seisan nordic format S-file.
    Returns an obspy.core.event.Catalog type: note this changed for version \
    0.1.0 from the inbuilt class types.

    :type sfile: str
    :param sfile: Path to the s-file

    :returns: :class: obspy.core.event.Event
    """
    import warnings
    from obspy.core.event import Event, Origin, Magnitude, Comment
    from obspy.core.event import EventDescription, CreationInfo
    f = open(sfile, 'r')
    # Base populate to allow for empty parts of file
    new_event = Event()
    topline = f.readline()
    if topline[79] == ' ' or topline[79] == '1':
        # Topline contains event information
        try:
            sfile_seconds = int(topline[16:18])
            if sfile_seconds == 60:
                sfile_seconds = 0
                add_seconds = 60
            else:
                add_seconds = 0
            new_event.origins.append(Origin())
            new_event.origins[0].time = UTCDateTime(int(topline[1:5]),
                                                    int(topline[6:8]),
                                                    int(topline[8:10]),
                                                    int(topline[11:13]),
                                                    int(topline[13:15]),
                                                    sfile_seconds,
                                                    int(topline[19:20])*100000)\
                + add_seconds
        except:
            warnings.warn("Couldn't read a date from sfile: "+sfile)
            new_event.origins.append(Origin(time=UTCDateTime(0)))
        # new_event.loc_mod_ind=topline[20]
        new_event.event_descriptions.append(EventDescription())
        new_event.event_descriptions[0].text = topline[21:23]
        # new_event.ev_id=topline[22]
        if not _float_conv(topline[23:30]) == 999:
            new_event.origins[0].latitude = _float_conv(topline[23:30])
            new_event.origins[0].longitude = _float_conv(topline[31:38])
            new_event.origins[0].depth = _float_conv(topline[39:43])
        # new_event.depth_ind = topline[44]
        # new_event.loc_ind = topline[45]
        new_event.creation_info = CreationInfo(agency_id=topline[45:48].
                                               strip())
        ksta = Comment(text='Number of stations=' +
                       topline[49:51].strip())
        new_event.origins[0].comments.append(ksta)
        # new_event.origins[0].nsta??? = _int_conv(topline[49:51])
        new_event.origins[0].time_errors['Time_Residual_RMS'] = \
            _float_conv(topline[52:55])
        # Read in magnitudes if they are there.
        if len(topline[59].strip()) > 0:
            new_event.magnitudes.append(Magnitude())
            new_event.magnitudes[0].mag = _float_conv(topline[56:59])
            new_event.magnitudes[0].magnitude_type = topline[59]
            new_event.magnitudes[0].creation_info = \
                CreationInfo(agency_id=topline[60:63].strip())
            new_event.magnitudes[0].origin_id = new_event.origins[0].\
                resource_id
        if len(topline[67].strip()) > 0:
            new_event.magnitudes.append(Magnitude())
            new_event.magnitudes[1].mag = _float_conv(topline[64:67])
            new_event.magnitudes[1].magnitude_type = topline[67]
            new_event.magnitudes[1].creation_info = \
                CreationInfo(agency_id=topline[68:71].strip())
            new_event.magnitudes[1].origin_id = new_event.origins[0].\
                resource_id
        if len(topline[75].strip()) > 0:
            new_event.magnitudes.append(Magnitude())
            new_event.magnitudes[2].mag = _float_conv(topline[72:75])
            new_event.magnitudes[2].magnitude_type = topline[75]
            new_event.magnitudes[2].creation_info = \
                CreationInfo(agency_id=topline[76:79].strip())
            new_event.magnitudes[2].origin_id = new_event.origins[0].\
                resource_id
    else:
        for line in f:
            if line[79] == '1':
                line = topline
                try:
                    new_event.origins.append(Origin())
                    new_event.origins[0].time = \
                        UTCDateTime(int(topline[1:5]),
                                    int(topline[6:8]),
                                    int(topline[8:10]),
                                    int(topline[11:13]),
                                    int(topline[13:15]),
                                    int(topline[16:18]),
                                    int(topline[19:20])*10)
                except:
                    new_event.origins.append(Origin(time=UTCDateTime(0)))
                # new_event.loc_mod_ind=topline[21]
                new_event.event_descriptions.append(EventDescription())
                new_event.event_descriptions[0].text = topline[21:23]
                # new_event.ev_id=topline[23]
                new_event.origins[0].latitude = _float_conv(topline[23:30])
                new_event.origins[0].longitude = _float_conv(topline[31:38])
                new_event.origins[0].depth = _float_conv(topline[39:43])
                # new_event.depth_ind = topline[44]
                # new_event.loc_ind = topline[45]
                new_event.creation_info = \
                    CreationInfo(agency_id=topline[45:48].strip())
                ksta = Comment(text='Number of stations=' +
                                    topline[49:51].strip())
                new_event.origins[0].comments.append(ksta)
                # new_event.origins[0].nsta??? = _int_conv(topline[49:51])
                new_event.origins[0].time_errors['Time_Residual_RMS'] = \
                    _float_conv(topline[52:55])
                # Read in magnitudes if they are there.
                if len(topline[59].strip()) > 0:
                    new_event.magnitudes.append(Magnitude())
                    new_event.magnitudes[0].mag = _float_conv(topline[56:59])
                    new_event.magnitudes[0].magnitude_type = topline[59]
                    new_event.magnitudes[0].creation_info = \
                        CreationInfo(agency_id=topline[60:63].strip())
                    new_event.magnitudes[0].origin_id = new_event.origins[0].\
                        resource_id
                if len(topline[67].strip()) > 0:
                    new_event.magnitudes.append(Magnitude())
                    new_event.magnitudes[1].mag = _float_conv(topline[64:67])
                    new_event.magnitudes[1].magnitude_type = topline[67]
                    new_event.magnitudes[1].creation_info = \
                        CreationInfo(agency_id=topline[68:71].strip())
                    new_event.magnitudes[1].origin_id = new_event.origins[0].\
                        resource_id
                if len(topline[75].strip()) > 0:
                    new_event.magnitudes.append(Magnitude())
                    new_event.magnitudes[2].mag = _float_conv(topline[72:75])
                    new_event.magnitudes[2].magnitude_type = topline[75]
                    new_event.magnitudes[2].creation_info = \
                        CreationInfo(agency_id=topline[76:79].strip())
                    new_event.magnitudes[2].origin_id = new_event.origins[0].\
                        resource_id
            if line[79] == '7':
                break
    f.close()
    # convert the nordic notation of magnitude to more general notation
    for _magnitude in new_event.magnitudes:
        _magnitude.magnitude_type = _nortoevmag(_magnitude.magnitude_type)
    return new_event


def readpicks(sfile):
    """
    Function to read pick information from the s-file and store this in an \
    obspy.event.Catalog type.  This was changed for version 0.1.0 from using \
    the inbuilt PICK class.

    :type sfile: String
    :param sfile: Path to sfile

    :return: obspy.core.event.Event

    .. warning:: Currently finalweight is unsupported, nor is velocity, \
    or angle of incidence.  This is because obspy.event stores slowness \
    in s/deg and takeoff angle, which would require computation from the \
    values stored in seisan.  Multiple weights are also not supported in \
    Obspy.event.
    """
    from obspy.core.event import Pick, WaveformStreamID, Arrival, Amplitude
    # Get wavefile name for use in resource_ids
    wav_names = readwavename(sfile)
    # First we need to read the header to get the timing info
    new_event = readheader(sfile)
    evtime = new_event.origins[0].time
    f = open(sfile, 'r')
    pickline = []
    # Set a default, ignored later unless overwritten
    SNR = 999
    if 'headerend' in locals():
        del headerend
    for lineno, line in enumerate(f):
        if 'headerend' in locals():
            if len(line.rstrip('\n').rstrip('\r')) in [80, 79] and \
               (line[79] == ' ' or line[79] == '4' or line[79] == '\n'):
                pickline += [line]
        elif line[79] == '7':
            header = line
            headerend = lineno
    amplitude_index = 0
    for pick_index, line in enumerate(pickline):
        if line[18:28].strip() == '':  # If line is empty miss it
            continue
        station = line[1:6].strip()
        channel = line[6:8].strip()
        network = 'NA'  # No network information provided in Sfile.
        weight = line[14]
        if weight == '_':
            phase = line[10:17]
            weight = 0
            polarity = ''
        else:
            phase = line[10:14].strip()
            polarity = line[16]
            if weight == ' ':
                weight = 0
        if polarity == '':
            polarity = "undecidable"
        elif polarity == 'C':
            polarity = "positive"
        elif polarity == 'D':
            polarity = 'negative'
        else:
            polarity = "undecidable"
        try:
            time = UTCDateTime(evtime.year, evtime.month, evtime.day,
                               int(line[18:20]), int(line[20:22]),
                               int(line[23:28].split('.')[0]),
                               int(line[23:28].split('.')[1])*10000)
        except (ValueError):
            time = UTCDateTime(evtime.year, evtime.month, evtime.day,
                               int(line[18:20]), int(line[20:22]), 0, 0)
            time += 60  # Add 60 seconds on to the time, this copes with s-file
            # preference to write seconds in 1-60 rather than 0-59 which
            # datetime objects accept
        coda = _int_conv(line[28:33])
        amplitude = _float_conv(line[33:40])
        peri = _float_conv(line[41:45])
        azimuth = _float_conv(line[46:51])
        velocity = _float_conv(line[52:56])
        if header[57:60] == 'AIN':
            AIN = _float_conv(line[57:60])
        elif header[57:60] == 'SNR':
            SNR = _float_conv(line[57:60])
        azimuthres = _int_conv(line[60:63])
        timeres = _float_conv(line[63:68])
        finalweight = _int_conv(line[68:70])
        distance = _float_conv(line[70:75])
        CAZ = _int_conv(line[76:79])
        # Create a new obspy.event.Pick class for this pick
        _waveform_id = WaveformStreamID(station_code=station,
                                        channel_code=channel,
                                        network_code=network)
        new_event.picks.append(Pick(waveform_id=_waveform_id,
                                    phase_hint=phase,
                                    polarity=polarity, time=time))
        if line[9] == 'I':
            new_event.picks[pick_index].onset = 'impulsive'
        elif line[9] == 'E':
            new_event.picks[pick_index].onset = 'emergent'
        if line[15] == 'A':
            new_event.picks[pick_index].evaluation_mode = 'automatic'
        else:
            new_event.picks[pick_index].evaluation_mode = 'manual'
        # Note these two are not always filled - velocity conversion not yet
        # implimented, needs to be converted from km/s to s/deg
        # if not velocity == 999.0:
            # new_event.picks[pick_index].horizontal_slowness = 1.0 / velocity
        if not azimuth == 999:
            new_event.picks[pick_index].backazimuth = azimuth
        del _waveform_id
        # Create new obspy.event.Amplitude class which references above Pick
        # only if there is an amplitude picked.
        if not amplitude == 999.0:
            new_event.amplitudes.append(Amplitude(generic_amplitude=amplitude,
                                                  period=peri,
                                                  pick_id=new_event.
                                                  picks[pick_index].
                                                  resource_id,
                                                  waveform_id=new_event.
                                                  picks[pick_index].
                                                  waveform_id))
            if new_event.picks[pick_index].phase_hint == 'IAML':
                # Amplitude for local magnitude
                new_event.amplitudes[amplitude_index].type = 'AML'
                # Set to be evaluating a point in the trace
                new_event.amplitudes[amplitude_index].category = 'point'
                # Default AML unit in seisan is nm (Page 139 of seisan
                # documentation, version 10.0)
                new_event.amplitudes[amplitude_index].generic_amplitude /=\
                    10**9
                new_event.amplitudes[amplitude_index].unit = 'm'
                new_event.amplitudes[amplitude_index].magnitude_hint = 'ML'
            else:
                # Generic amplitude type
                new_event.amplitudes[amplitude_index].type = 'A'
            if not SNR == 999.0:
                new_event.amplitudes[amplitude_index].snr = SNR
            amplitude_index += 1
        elif not coda == 999:
            # Create an amplitude instance for code duration also
            new_event.amplitudes.append(Amplitude(generic_amplitude=coda,
                                                  pick_id=new_event.
                                                  picks[pick_index].
                                                  resource_id,
                                                  waveform_id=new_event.
                                                  picks[pick_index].
                                                  waveform_id))
            # Amplitude for coda magnitude
            new_event.amplitudes[amplitude_index].type = 'END'
            # Set to be evaluating a point in the trace
            new_event.amplitudes[amplitude_index].category = 'duration'
            new_event.amplitudes[amplitude_index].unit = 's'
            new_event.amplitudes[amplitude_index].magnitude_hint = 'Mc'
            if SNR and not SNR == 999.0:
                new_event.amplitudes[amplitude_index].snr = SNR
            amplitude_index += 1
        # Create new obspy.event.Arrival class referencing above Pick
        new_event.origins[0].arrivals.append(Arrival(phase=new_event.
                                                     picks[pick_index].
                                                     phase_hint,
                                                     pick_id=new_event.
                                                     picks[pick_index].
                                                     resource_id))
        if weight != 999:
            new_event.origins[0].arrivals[pick_index].time_weight =\
                weight
        if azimuthres != 999:
            new_event.origins[0].arrivals[pick_index].backazimuth_residual =\
                azimuthres
        if timeres != 999:
            new_event.origins[0].arrivals[pick_index].time_residual =\
                timeres
        if distance != 999:
            new_event.origins[0].arrivals[pick_index].distance =\
                distance
        if CAZ != 999:
            new_event.origins[0].arrivals[pick_index].azimuth =\
                CAZ
    f.close()
    # Write event to catalog object for ease of .write() method
    return new_event


def readwavename(sfile):
    """
    Convenience function to extract the waveform filename from the s-file, \
    returns a list of waveform names found in the s-file as multiples can \
    be present.

    :type sfile: str
    :param sfile: Path to the sfile

    :returns: List of str
    """
    f = open(sfile)
    wavename = []
    for line in f:
        if len(line) == 81 and line[79] == '6':
            wavename.append(line[1:79].strip())
    f.close()
    return wavename


def blanksfile(wavefile, evtype, userID, outdir, overwrite=False,
               evtime=False):
    """
    Module to generate an empty s-file with a populated header for a given \
    waveform.

    :type wavefile: String
    :param wavefile: Wavefile to associate with this S-file, the timing of \
        the S-file will be taken from this file if evtime is not set.
    :type evtype: String
    :param evtype: L,R,D
    :type userID: String
    :param userID: 4-charectar SEISAN USER ID
    :type outdir: String
    :param outdir: Location to write S-file
    :type overwrite: Bool
    :param overwrite: Overwrite an existing S-file, default=False
    :type evtime: UTCDateTime
    :param evtime: If given this will set the timing of the S-file

    :returns: String, S-file name
    """

    from obspy import read as obsread
    import sys
    import os
    import datetime

    if not evtime:
        try:
            st = obsread(wavefile)
            evtime = st[0].stats.starttime
        except:
            print('Wavefile: '+wavefile +
                  ' is invalid, try again with real data.')
            sys.exit()
    # Check that user ID is the correct length
    if len(userID) != 4:
        print 'User ID must be 4 characters long'
        sys.exit()
    # Check that outdir exists
    if not os.path.isdir(outdir):
        print 'Out path does not exist, I will not create this: '+outdir
        sys.exit()
    # Check that evtype is one of L,R,D
    if evtype not in ['L', 'R', 'D']:
        print 'Event type must be either L, R or D'
        sys.exit()

    # Generate s-file name in the format dd-hhmm-ss[L,R,D].Syyyymm
    sfile = outdir + '/' + str(evtime.day).zfill(2) + '-' +\
        str(evtime.hour).zfill(2) +\
        str(evtime.minute).zfill(2) + '-' +\
        str(evtime.second).zfill(2) + evtype + '.S' +\
        str(evtime.year) +\
        str(evtime.month).zfill(2)
    # Check is sfile exists
    if os.path.isfile(sfile) and not overwrite:
        print 'Desired sfile: ' + sfile + ' exists, will not overwrite'
        for i in range(1, 10):
            sfile = outdir + '/' + str(evtime.day).zfill(2) + '-' +\
                str(evtime.hour).zfill(2) +\
                str(evtime.minute).zfill(2) + '-' +\
                str(evtime.second+i).zfill(2) + evtype + '.S' +\
                str(evtime.year) +\
                str(evtime.month).zfill(2)
            if not os.path.isfile(sfile):
                break
        else:
            print 'Tried generated files up to 20s in advance and found they'
            print 'all exist, you need to clean your stuff up!'
            sys.exit()
        # sys.exit()
    f = open(sfile, 'w')
    # Write line 1 of s-file
    f.write(' ' + str(evtime.year) + ' ' +
            str(evtime.month).rjust(2) +
            str(evtime.day).rjust(2) + ' ' +
            str(evtime.hour).rjust(2) +
            str(evtime.minute).rjust(2) + ' ' +
            str(float(evtime.second)).rjust(4) + ' ' +
            evtype + '1'.rjust(58) + '\n')
    # Write line 2 of s-file
    f.write(' ACTION:ARG ' + str(datetime.datetime.now().year)[2:4] + '-' +
            str(datetime.datetime.now().month).zfill(2) + '-' +
            str(datetime.datetime.now().day).zfill(2) + ' ' +
            str(datetime.datetime.now().hour).zfill(2) + ':' +
            str(datetime.datetime.now().minute).zfill(2) + ' OP:' +
            userID.ljust(4) + ' STATUS:'+'ID:'.rjust(18) +
            str(evtime.year) +
            str(evtime.month).zfill(2) +
            str(evtime.day).zfill(2) +
            str(evtime.hour).zfill(2) +
            str(evtime.minute).zfill(2) +
            str(evtime.second).zfill(2) +
            'I'.rjust(6) + '\n')
    # Write line 3 of s-file
    f.write(' ' + wavefile + '6'.rjust(79-len(wavefile)) + '\n')
    # Write final line of s-file
    f.write(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU' +
            ' VELO AIN AR TRES W  DIS CAZ7\n')
    f.close()
    print 'Written s-file: ' + sfile
    return sfile


def eventtoSfile(event, userID, evtype, outdir, wavefiles, explosion=False,
                 overwrite=False):
    """
    Function to take an obspy.event and write the relevant information to a \
    nordic formatted s-file

    :type event: obspy.event.core.Catalog
    :param event: A single obspy event
    :type userID: str
    :param userID: Up to 4 character user ID
    :type evtype: str
    :param evtype: Single character string to describe the event, either L, R \
        or D.
    :type outdir: str
    :param outdir: Path to directory to write to
    :type wavefiles: list of str
    :param wavefiles: Waveforms to associate the sfile with
    :type explosion: bool
    :param explosion: Note if the event is an explosion, will be marked by an \
        E.
    :type overwrite: bool
    :param overwrite: force to overwrite old files, defaults to False

    :returns: str: name of sfile written
    """
    import datetime
    import os
    from obspy.core.event import Catalog, Event
    # First we need to work out what to call the s-file and open it
    # Check that user ID is the correct length
    if len(userID) != 4:
        raise IOError('User ID must be 4 characters long')
    # Check that outdir exists
    if not os.path.isdir(outdir):
        raise IOError('Out path does not exist, I will not create this: ' +
                      outdir)
    # Check that evtype is one of L,R,D
    if evtype not in ['L', 'R', 'D']:
        raise IOError('Event type must be either L, R or D')
    if explosion:
        evtype += 'E'
    # Check that there is one event
    if isinstance(event, Catalog) and len(event) == 1:
        event = event[0]
    elif isinstance(event, Event):
        event = event
    else:
        raise IOError('Needs a single event')
    if isinstance(wavefiles, str):
        wavefiles = [wavefiles]
    elif isinstance(wavefiles, list):
        wavefiles = wavefiles
    else:
        raise IOError(wavefiles + ' is neither string or list')
    # Determine name from origin time
    sfilename = event.origins[0].time.datetime.strftime('%d-%H%M-%S') +\
        evtype + '.S' + event.origins[0].time.datetime.strftime('%Y%m')
    evtime = event.origins[0].time
    # Check that the file doesn't exist
    if not overwrite and os.path.isfile(outdir+'/'+sfilename):
        raise IOError(outdir+'/'+sfilename +
                      ' already exists, will not overwrite')
    sfile = open(outdir + '/' + sfilename, 'w')
    # Write the header info.
    lat = '{0:.3f}'.format(event.origins[0].get('latitude')) or ''
    lon = '{0:.3f}'.format(event.origins[0].get('longitude')) or ''
    depth = '{0:.1f}'.format(event.origins[0].get('depth')) or ''
    agency = event.creation_info.get('agency_id') or ''
    timerms = '{0:.1f}'.format(event.origins[0].time_errors.Time_Residual_RMS)\
        or ''
    mag_1 = '{0:.1f}'.format(event.magnitudes[0].mag) or ''
    mag_1_type = _evmagtonor(event.magnitudes[0].magnitude_type) or ''
    mag_1_agency = event.magnitudes[0].creation_info.agency_id or ''
    mag_2 = '{0:.1f}'.format(event.magnitudes[1].mag) or ''
    mag_3 = '{0:.1f}'.format(event.magnitudes[2].mag) or ''
    mag_2_type = _evmagtonor(event.magnitudes[1].magnitude_type) or ''
    mag_2_agency = event.magnitudes[1].creation_info.agency_id or ''
    mag_3_type = _evmagtonor(event.magnitudes[2].magnitude_type) or ''
    mag_3_agency = event.magnitudes[2].creation_info.agency_id or ''
    # Work out how many stations were used
    if len(event.picks) > 0:
        stations = [pick.waveform_id.station_code for pick in event.picks]
        ksta = str(len(list(set(stations))))
    else:
        ksta = ''
    sfile.write(' ' + str(evtime.year) + ' ' +
                str(evtime.month).rjust(2) +
                str(evtime.day).rjust(2) + ' ' +
                str(evtime.hour).rjust(2) +
                str(evtime.minute).rjust(2) + ' ' +
                str(float(evtime.second)).rjust(4) + ' ' +
                evtype.ljust(2) + lat.rjust(7) + ' ' + lon.rjust(7) +
                depth.rjust(5) + agency.rjust(5) + ksta.rjust(3) +
                timerms.rjust(4) +
                mag_1.rjust(4) + mag_1_type.rjust(1) +
                mag_1_agency[0:3].rjust(3) +
                mag_2.rjust(4) + mag_2_type.rjust(1) +
                mag_2_agency[0:3].rjust(3) +
                mag_3.rjust(4) + mag_3_type.rjust(1) +
                mag_3_agency[0:3].rjust(3) + '1' + '\n')
    # Write line 2 of s-file
    sfile.write(' ACTION:ARG ' + str(datetime.datetime.now().year)[2:4] + '-' +
                str(datetime.datetime.now().month).zfill(2) + '-' +
                str(datetime.datetime.now().day).zfill(2) + ' ' +
                str(datetime.datetime.now().hour).zfill(2) + ':' +
                str(datetime.datetime.now().minute).zfill(2) + ' OP:' +
                userID.ljust(4) + ' STATUS:'+'ID:'.rjust(18) +
                str(evtime.year) +
                str(evtime.month).zfill(2) +
                str(evtime.day).zfill(2) +
                str(evtime.hour).zfill(2) +
                str(evtime.minute).zfill(2) +
                str(evtime.second).zfill(2) +
                'I'.rjust(6) + '\n')
    # Write line 3 of s-file
    for wavefile in wavefiles:
        sfile.write(' ' + wavefile + '6'.rjust(79-len(wavefile)) + '\n')
    # Write final line of s-file
    sfile.write(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU' +
                ' VELO AIN AR TRES W  DIS CAZ7\n')
    sfile.close()
    # Now call the populateSfile function
    if len(event.picks) > 0:
        populateSfile(sfilename, event)
    return sfilename


def populateSfile(sfile, event):
    """
    Module to populate a blank nordic format S-file with pick information, \
    arguments required are the filename of the blank s-file and the picks \
    where picks is a dictionary of picks including station, channel, \
    impulsivity, phase, weight, polarity, time, coda, amplitude, peri, \
    azimuth, velocity, SNR, azimuth residual, Time-residual, final weight, \
    epicentral distance & azimuth from event to station.

    This is a full pick line information from the seisan manual, P. 341

    :type sfile: str
    :param sfile: Path to S-file to populate, must have a header already
    :type event: :class: obspy.event.core.Catalog
    :param picks: A single event to be written to a single S-file.
    """
    from obspy.core.event import Catalog, Event
    # first check that the event is only one event
    if isinstance(event, Catalog) and len(event) == 1:
        event = event[0]
    elif isinstance(event, Event):
        event = event
    else:
        raise AttributeError('More than one event in the catalog, use a ' +
                             ' different method')
    f = open(sfile, 'r')
    # Find type 7 line, under which picks should be - if there are already
    # picks there we should preserve them
    body = ''
    header = ''
    if 'headerend' in locals():
        del headerend
    for lineno, line in enumerate(f):
        identifier = line[79]
        if 'headerend' in locals():
            body += line
        else:
            header += line
        if identifier == '7':
            headerend = lineno
    f.close()
    #
    # Now generate lines for the new picks
    newpicks = '\n'.join(nordpick(event))
    # Write all new and old info back in
    f = open(sfile, 'w')
    f.write(header)
    f.write(body)
    f.write(newpicks + '\n')
    f.write('\n'.rjust(81))
    # f.write('\n')
    f.close()
    return


def eventtopick(event):
    """
    Wrapper function to convert from obspy.core.event to legacy PICK and \
    EVENT classes.

    :type event: obspy.core.event.Event
    :param event: A single obspy event

    :returns: List of PICK(), and a single EVENTINFO()

    .. note:: This is a wrapper to simplify transition from PICK and \
    EVENT classes to obspy.core.event classes.  This will not be maintained \
    beyond v 0.1.0.

    .. versionadded:: 0.1.0
    """
    # Check that the event is a single event
    from obspy.core.event import Catalog, Event
    # first check that the event is only one event
    if isinstance(event, Catalog) and len(event) == 1:
        event = event[0]
    elif isinstance(event, Event):
        event = event
    else:
        raise AttributeError('More than one event in the catalog, use a ' +
                             ' different method')
    stations = [pick.waveform_id.station_code for pick in event.picks]
    nsta = len(list(set(stations)))
    # Generate the EVENTINFO object
    event_descriptions = event.event_descriptions[0].text or '   '
    if len(event_descriptions) == 2:
        event_descriptions = event_descriptions.rjust(3)
    evinfo = EVENTINFO(time=event.origins[0].time,
                       loc_mod_ind=event_descriptions[0],
                       dist_ind=event_descriptions[1],
                       ev_id=event_descriptions[2],
                       latitude=event.origins[0].latitude,
                       longitude=event.origins[0].longitude,
                       depth=event.origins[0].depth,
                       depth_ind=' ',
                       loc_ind=' ',
                       agency=event.creation_info.get('agency_id') or ' ',
                       nsta=nsta,
                       t_RMS=event.origins[0].time_errors.Time_Residual_RMS or
                       float('NaN'),
                       Mag_1=event.magnitudes[0].mag or float('NaN'),
                       Mag_1_type=_evmagtonor(event.magnitudes[0].
                                              magnitude_type) or ' ',
                       Mag_1_agency=event.magnitudes[0].creation_info.agency_id
                       or ' ',
                       Mag_2=event.magnitudes[1].mag or float('NaN'),
                       Mag_2_type=_evmagtonor(event.magnitudes[1].
                                              magnitude_type) or ' ',
                       Mag_2_agency=event.magnitudes[1].creation_info.agency_id
                       or ' ',
                       Mag_3=event.magnitudes[2].mag or float('NaN'),
                       Mag_3_type=_evmagtonor(event.magnitudes[2].
                                              magnitude_type) or ' ',
                       Mag_3_agency=event.magnitudes[2].creation_info.agency_id
                       or ' ')
    # Can make use of nordpick, which will remain in place for many versions?
    pick_strings = nordpick(event)
    # Then convert from pick-strings to PICK class
    picks = []
    evtime = event.origins[0].time
    for line in pick_strings:
        # Copied from old readpicks function
        station = line[1:6].strip()
        channel = line[6:8].strip()
        impulsivity = line[9]
        weight = line[14]
        if weight == '_':
            phase = line[10:17]
            weight = ''
            polarity = ''
        else:
            phase = line[10:14].strip()
            polarity = line[16]
        try:
            time = UTCDateTime(evtime.year, evtime.month, evtime.day,
                               int(line[18:20]), int(line[20:22]),
                               int(line[23:28].split('.')[0]),
                               int(line[23:28].split('.')[1]) * 10000)
            # Includes possible bug with seconds not being alligned.
            # Shouldn't happen in this case, but best to include for possible
            # copy/paste later!
        except (ValueError):
            time = UTCDateTime(evtime.year, evtime.month, evtime.day,
                               int(line[18:20]), int(line[20:22]), 0, 0)
            time += 60  # Add 60 seconds on to the time, this copes with s-file
        weight = _int_conv(weight)
        coda = _int_conv(line[28:33])
        amplitude = _float_conv(line[33:40])
        peri = _float_conv(line[41:45])
        azimuth = _float_conv(line[46:51])
        velocity = _float_conv(line[52:56])
        azimuthres = _int_conv(line[60:63])
        timeres = _float_conv(line[63:68])
        finalweight = _int_conv(line[68:70])
        distance = _float_conv(line[70:75])
        CAZ = _int_conv(line[76:79])
        AIN = 999
        SNR = float('NaN')
        picks += [PICK(station, channel, impulsivity, phase, weight, polarity,
                       time, coda, amplitude, peri, azimuth, velocity, AIN,
                       SNR, azimuthres, timeres, finalweight, distance, CAZ)]
    return picks, evinfo


def picktoevent(evinfo, picks):
    """
    Wrapper function to convert from EVENTINFO and PICK classes to \
    obspy.core.event.Event.

    :type evinfo: EVENTINFO
    :param evinfo: Event header info for a single event
    :type picks: List of PICK
    :param picks: List of picks associated with the event

    :returns: obspy.core.event.Event

    .. note:: This is a legacy support function, users should avoid \
    this as it will be removed for version 0.1.1.  Written to aid transition \
    from in-built classes to obspy.core.event classes.

    .. versionadded:: 0.1.0
    """
    from obspy.core.event import Event, Origin, Magnitude, Comment
    from obspy.core.event import EventDescription, CreationInfo
    from obspy.core.event import Pick, WaveformStreamID, Arrival, Amplitude
    # Cope with possible single pick case
    if not isinstance(picks, list):
        picks = [picks]
    # Convert the relevant evinfo fields to an Event instance
    event = Event()
    event.origins.append(Origin())
    event.origins[0].time = evinfo.time
    event.event_descriptions.append(EventDescription())
    event.event_descriptions[0].text = ''.join([evinfo.loc_mod_ind,
                                                evinfo.dist_ind,
                                                evinfo.ev_id]).strip()
    event.origins[0].latitude = evinfo.latitude
    event.origins[0].longitude = evinfo.longitude
    event.origins[0].depth = evinfo.depth
    event.creation_info = CreationInfo(agency_id=evinfo.agency)
    event.origins[0].comments.append(Comment(text='Number of stations=' +
                                             str(evinfo.nsta)))
    event.origins[0].time_errors['Time_Residual_RMS'] = \
        evinfo.t_RMS
    event.magnitudes.append(Magnitude())
    event.magnitudes[0].mag = evinfo.Mag_1
    event.magnitudes[0].magnitude_type = _nortoevmag(evinfo.Mag_1_type)
    event.magnitudes[0].creation_info = CreationInfo(agency_id=evinfo.
                                                     Mag_1_agency)
    event.magnitudes[0].origin_id = event.origins[0].resource_id
    event.magnitudes.append(Magnitude())
    event.magnitudes[1].mag = evinfo.Mag_2
    event.magnitudes[1].magnitude_type = _nortoevmag(evinfo.Mag_2_type)
    event.magnitudes[1].creation_info = CreationInfo(agency_id=evinfo.
                                                     Mag_2_agency)
    event.magnitudes[1].origin_id = event.origins[0].resource_id
    event.magnitudes.append(Magnitude())
    event.magnitudes[2].mag = evinfo.Mag_3
    event.magnitudes[2].magnitude_type = _nortoevmag(evinfo.Mag_3_type)
    event.magnitudes[2].creation_info = CreationInfo(agency_id=evinfo.
                                                     Mag_3_agency)
    event.magnitudes[2].origin_id = event.origins[0].resource_id
    # We now have all the header info converted that we can hold in EVENT class
    # Move on to the picks.
    amplitude_index = 0
    for pick_index, pick in enumerate(picks):
        _waveform_id = WaveformStreamID(station_code=pick.station,
                                        channel_code=pick.channel,
                                        network_code='NA')
        if pick.polarity == '':
            polarity = "undecidable"
        elif pick.polarity == 'C':
            polarity = "positive"
        elif pick.polarity == 'D':
            polarity = 'negative'
        else:
            polarity = "undecidable"
        event.picks.append(Pick(waveform_id=_waveform_id,
                                phase_hint=pick.phase,
                                polarity=polarity,
                                time=pick.time))
        if pick.impulsivity == 'I':
            event.picks[pick_index].onset = 'impulsive'
        elif pick.impulsivity == 'E':
            event.picks[pick_index].onset = 'emergent'
        if not np.isnan(pick.azimuth):
            event.picks[pick_index].backazimuth = pick.azimuth
        del _waveform_id
        if not np.isnan(pick.amplitude):
            event.amplitudes.append(Amplitude(generic_amplitude=pick.amplitude,
                                              period=pick.peri,
                                              pick_id=event.picks[pick_index].
                                              resource_id,
                                              waveform_id=event.
                                              picks[pick_index].waveform_id))
            if event.picks[pick_index].phase_hint == 'IAML':
                event.amplitudes[amplitude_index].type = 'AML'
                # Set to be evaluating a point in the trace
                event.amplitudes[amplitude_index].category = 'point'
                # Default AML unit in seisan is nm (Page 139 of seisan
                # documentation, version 10.0)
                event.amplitudes[amplitude_index].generic_amplitude /=\
                    10**9
                event.amplitudes[amplitude_index].unit = 'm'
                event.amplitudes[amplitude_index].magnitude_hint = 'ML'
            else:
                # Generic amplitude type
                event.amplitudes[amplitude_index].type = 'A'
            if not np.isnan(pick.SNR):
                event.amplitudes[amplitude_index].snr = pick.SNR
            amplitude_index += 1
        elif not pick.coda == 999:
            # Create an amplitude instance for code duration also
            event.amplitudes.append(Amplitude(generic_amplitude=pick.coda,
                                              pick_id=event.
                                              picks[pick_index].resource_id,
                                              waveform_id=event.
                                              picks[pick_index].waveform_id))
            # Amplitude for coda magnitude
            event.amplitudes[amplitude_index].type = 'END'
            # Set to be evaluating a point in the trace
            event.amplitudes[amplitude_index].category = 'duration'
            event.amplitudes[amplitude_index].unit = 's'
            event.amplitudes[amplitude_index].magnitude_hint = 'Mc'
            if pick.SNR and not np.isnan(pick.SNR):
                event.amplitudes[amplitude_index].snr = pick.SNR
            amplitude_index += 1
        # Generate Arrival objects for the pick
        event.origins[0].arrivals.append(Arrival(phase=event.picks[pick_index].
                                                 phase_hint,
                                                 pick_id=event.
                                                 picks[pick_index].
                                                 resource_id))
        if pick.weight != 999:
            event.origins[0].arrivals[pick_index].time_weight =\
                pick.weight
        if pick.azimuthres != 999:
            event.origins[0].arrivals[pick_index].backazimuth_residual =\
                pick.azimuthres
        if not np.isnan(pick.timeres):
            event.origins[0].arrivals[pick_index].time_residual =\
                pick.timeres
        if not np.isnan(pick.distance):
            event.origins[0].arrivals[pick_index].distance =\
                pick.distance
        if pick.CAZ != 999:
            event.origins[0].arrivals[pick_index].azimuth =\
                pick.CAZ
    return event


def nordpick(event):
    """
    Function to print information from an obspy.event class to nordic format.

    :type event: :class: obspy.core.event.Event
    :param event: A single obspy event.

    :returns: List of String

    .. note:: Currently finalweight is unsupported, nor is velocity, or \
    angle of incidence.  This is because obspy.event stores slowness in \
    s/deg and takeoff angle, which would require computation from the \
    values stored in seisan.  Multiple weights are also not supported in \
    Obspy.event.

    .. versionadded:: 0.1.0
    """

    pick_strings = []
    for pick in event.picks:
        # Convert string to short sting
        if pick.onset == 'impulsive':
            impulsivity = 'I'
        elif pick.onset == 'emergent':
            impulsivity = 'E'
        else:
            impulsivity = ' '

        # Convert string to short string
        if pick.polarity == 'positive':
            polarity = 'C'
        elif pick.polarity == 'negative':
            polarity = 'D'
        else:
            polarity = ' '
        # Extract velocity: Note that horizontal slowness in quakeML is stored
        # as s/deg
        if pick.horizontal_slowness:
            # velocity = 1.0 / pick.horizontal_slowness
            velocity = ' '  # Currently this conversion is unsupported.
        else:
            velocity = ' '
        # Extract azimuth
        if pick.backazimuth:
            azimuth = pick.backazimuth
        else:
            azimuth = ' '
        # Extract the correct arrival info for this pick - assuming only one
        # arrival per pick...
        arrival = [arrival for arrival in event.origins[0].arrivals
                   if arrival.pick_id == pick.resource_id]
        if len(arrival) > 0:
            arrival = arrival[0]
            # Extract weight - should be stored as 0-4, or 9 for seisan.
            if arrival.time_weight:
                weight = int(arrival.time_weight)
            else:
                weight = '0'
            # Extract azimuth residual
            if arrival.backazimuth_residual:
                azimuthres = int(arrival.backazimuth_residual)
            else:
                azimuthres = ' '
            # Extract time residual
            if arrival.time_residual:
                timeres = arrival.time_residual
            else:
                timeres = ' '
            # Extract distance
            if arrival.distance:
                distance = arrival.distance
                if distance >= 100.0:
                    distance = _int_conv(distance)
                elif 10.0 < distance < 100.0:
                    distance = round(distance, 1)
                    round_len = 1
                elif distance < 10.0:
                    distance = round(distance, 2)
                    round_len = 2
                else:
                    round_len = False
            else:
                distance = ' '
                round_len = False
            # Extract CAZ
            if arrival.azimuth:
                CAZ = int(arrival.azimuth)
            else:
                CAZ = ' '
        else:
            CAZ = ' '
            round_len = False
            distance = ' '
            timeres = ' '
            azimuthres = ' '
            azimuth = ' '
            weight = ' '
        # Extract amplitude: note there can be multiple amplitudes, but they
        # should be associated with different picks.
        amplitude = [amplitude for amplitude in event.amplitudes
                     if amplitude.pick_id == pick.resource_id]
        if len(amplitude) > 0:
            amplitude = amplitude[0]
            # Determine type of amplitude
            if amplitude.type != 'END':
                # Extract period
                if amplitude.period:
                    peri = amplitude.period
                    if peri < 10.0:
                        peri_round = 2
                    elif peri >= 10.0:
                        peri_round = 1
                    else:
                        peri_round = False
                else:
                    peri = ' '
                # Extract amplitude and convert units
                if amplitude.generic_amplitude:
                    amp = amplitude.generic_amplitude
                    if amplitude.unit in ['m', 'm/s', 'm/(s*s)', 'm*s']:
                        amp *= 10**9
                    # Otherwise we will assume that the amplitude is in counts
                else:
                    amp = np.nan
                coda = ' '
            else:
                coda = int(amplitude.generic_amplitude)
                peri = ' '
                amp = np.nan
        else:
            peri = ' '
            peri_round = False
            amp = np.nan
            coda = ' '
        # If the weight is 0 then we don't need to print it
        if weight == 0 or weight == '0':
            weight = 999  # this will return an empty string using _str_conv
        # Generate a print string and attach it to the list
        pick_strings.append(' ' + pick.waveform_id.station_code.ljust(5) +
                            pick.waveform_id.channel_code[0] +
                            pick.waveform_id.channel_code[len(pick.waveform_id.
                                                              channel_code)
                                                          - 1] +
                            ' ' + impulsivity + pick.phase_hint.ljust(4) +
                            _str_conv(int(weight)).rjust(1) + ' ' +
                            polarity.rjust(1) + ' ' +
                            str(pick.time.hour).rjust(2) +
                            str(pick.time.minute).rjust(2) +
                            str(pick.time.second).rjust(3) + '.' +
                            str(float(pick.time.microsecond) /
                            (10 ** 4)).split('.')[0].zfill(2) +
                            _str_conv(coda).rjust(5)[0:5] +
                            _str_conv(round(amp, 1)).rjust(7)[0:7] +
                            _str_conv(peri, rounded=peri_round).rjust(5) +
                            _str_conv(azimuth).rjust(6) +
                            _str_conv(velocity).rjust(5) +
                            _str_conv(' ').rjust(4) +
                            _str_conv(azimuthres).rjust(3) +
                            _str_conv(timeres, rounded=2).rjust(5) +
                            _str_conv(' ').rjust(2) +
                            _str_conv(distance, rounded=round_len).rjust(5) +
                            _str_conv(CAZ).rjust(4)+' ')
        # Note that currently finalweight is unsupported, nor is velocity, or
        # angle of incidence.  This is because obspy.event stores slowness in
        # s/deg and takeoff angle, which would require computation from the
        # values stored in seisan.  Multiple weights are also not supported in
        # Obspy.event
    return pick_strings


if __name__ == "__main__":
    import doctest
    doctest.testmod()

# if __name__ == '__main__':
#     # Read arguments
#     import sys
#     if len(sys.argv) != 6:
#         print('Requires 5 arguments: wavefile, evtype, userID, outdir,'
#               ' overwrite')
#         sys.exit()
#     else:
#         wavefile = str(sys.argv[1])
#         evtype = str(sys.argv[2])
#         userID = str(sys.argv[3])
#         outdir = str(sys.argv[4])
#         overwrite = str(sys.argv[5])
#     sfile = blanksfile(wavefile, evtype, userID, outdir, overwrite)
#     print sfile
