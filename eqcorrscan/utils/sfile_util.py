"""
Part of the EQcorrscan module to read nordic format s-files and write them.
Maps between Nordic format and obspy Event objects.

Note that these functions do not provide full functionality between quakeML
and seisan s-files.  Currently (as of version 0.1.1) these only convert pick
times and phase information, along with amplitude information for local
magnitudes between seisan and quakeML.  Location information including
hypocentre, origin time and magnitudes are also handled.

A series of wrappers and conversions is included between the legacy PICK
and EVENTINFO classes, however these will be depreciated along with these
classes for version 0.1.0.  Users should transition to using obspy.core.event
classes as these have more support and functionality.

We have not implemented any handling of focal mechanism solutions between
the two formats.

.. note:: Pick time-residuals are handled in event.origins[0].arrivals, with \
    the arrival.pick_id linking the arrival (which contain calculated \
    information) with the pick.resource_id (where the pick contains only \
    physical measured information).

:Example:

>>> from eqcorrscan.utils import sfile_util
>>> event = sfile_util.readpicks('eqcorrscan/tests/test_data/REA/TEST_/' +
...                              '01-0411-15L.S201309')
>>> pick = event.picks[0]
>>> time_rms = [arrival.time_residual for arrival in event.origins[0].arrivals
...             if arrival.pick_id == pick.resource_id]

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
from six import string_types
from obspy import UTCDateTime
import numpy as np
import warnings


def _int_conv(string):
    """
    Convenience tool to convert from string to integer.
    If empty string return a 999 rather than an error.

    >>> _int_conv('12')
    12
    >>> _int_conv('')
    999
    """
    try:
        intstring = int(string)
    except:
        intstring = 999
    return intstring


def _float_conv(string):
    """
    Convenience tool to convert from string to float.
    If empty string return NaN rather than an error.

    >>> _float_conv('12')
    12.0
    >>> _float_conv('')
    999.0
    >>> _float_conv('12.324')
    12.324
    """
    try:
        floatstring = float(string)
    except:
        # cjh Changed this to deal with QuakeML validation issues with 'NaN'
        floatstring = float('999')
    return floatstring


def _str_conv(number, rounded=False):
    """
    Convenience tool to convert a number, either float or int into a string.
    If the int is 999, or the float is NaN, returns empty string.

    >>> _str_conv(12.3)
    '12.3'
    >>> _str_conv(12.34546, rounded=1)
    '12.3'
    """
    if (isinstance(number, float) and np.isnan(number)) or number == 999:
        string = ' '
    elif isinstance(number, string_types):
        return str(number)
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
    return str(string)


def _evmagtonor(mag_type):
    """
    Switch from obspy event magnitude types to seisan syntax.

    >>> _evmagtonor('mB')
    'b'
    """
    if mag_type in ['ML', 'MLv']:
        # MLv is local magnitude on vertical component
        mag = 'L'
    elif mag_type == 'mB':
        mag = 'b'
    elif mag_type == 'Ms':
        mag = 's'
    elif mag_type == 'MS':
        mag = 'S'
    elif mag_type in ['MW', 'Mw']:
        mag = 'W'
    elif mag_type == 'MbLg':
        mag = 'G'
    elif mag_type in ['Mc', 'MC']:
        mag = 'C'
    elif mag_type == 'M':
        mag = 'W'  # Convert generic magnitude to moment magnitude
        msg = ('Converting generic magnitude to being stored as moment mag')
        warnings.warn(msg)
    else:
        warnings.warn(mag_type + ' is not convertable')
        return ''
    return str(mag)


def _nortoevmag(mag_type):
    """
    Switch from nordic type magnitude notation to obspy event magnitudes.

    >>> _nortoevmag('b')
    'mB'
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
        warnings.warn(mag_type + ' is not convertable')
        return ''
    return str(mag)


def readheader(sfile):
    """
    Read header information from a seisan nordic format S-file.
    Returns an obspy.core.event.Catalog type: note this changed for version \
    0.1.0 from the inbuilt class types.

    :type sfile: str
    :param sfile: Path to the s-file

    :returns: :class: obspy.core.event.Event

    >>> event = readheader('eqcorrscan/tests/test_data/REA/TEST_/' +
    ...                    '01-0411-15L.S201309')
    >>> print(event.origins[0].time)
    2013-09-01T04:11:15.700000Z
    """
    import warnings
    from obspy.core.event import Event, Origin, Magnitude, Comment
    from obspy.core.event import EventDescription, CreationInfo
    f = open(sfile, 'r')
    # Base populate to allow for empty parts of file
    new_event = Event()
    topline = f.readline()
    if not len(topline.rstrip()) == 80:
        raise IOError('s-file has a corrupt header, not 80 char long')
    f.seek(0)
    for line in f:
        if line[79] in [' ', '1']:
            topline = line
            break
        if line[79] == '7':
            raise IOError('No header found, corrupt s-file?')
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
                                                int(topline[19:20]) *
                                                100000)\
            + add_seconds
    except:
        warnings.warn("Couldn't read a date from sfile: " + sfile)
        new_event.origins.append(Origin(time=UTCDateTime(0)))
    # new_event.loc_mod_ind=topline[20]
    new_event.event_descriptions.append(EventDescription())
    new_event.event_descriptions[0].text = topline[21:23]
    # new_event.ev_id=topline[22]
    if not _float_conv(topline[23:30]) == 999:
        new_event.origins[0].latitude = _float_conv(topline[23:30])
        new_event.origins[0].longitude = _float_conv(topline[31:38])
        new_event.origins[0].depth = _float_conv(topline[39:43]) * 1000
    else:
        # The origin 'requires' a lat & long
        new_event.origins[0].latitude = float('NaN')
        new_event.origins[0].longitude = float('NaN')
        new_event.origins[0].depth = float('NaN')
    # new_event.depth_ind = topline[44]
    # new_event.loc_ind = topline[45]
    new_event.creation_info = CreationInfo(agency_id=topline[45:48].
                                           strip())
    ksta = Comment(text='Number of stations=' +
                   topline[49:51].strip())
    new_event.origins[0].comments.append(ksta)
    # new_event.origins[0].nsta??? = _int_conv(topline[49:51])
    if not _float_conv(topline[51:55]) == 999:
        new_event.origins[0].time_errors['Time_Residual_RMS'] = \
            _float_conv(topline[51:55])
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
    f.close()
    # convert the nordic notation of magnitude to more general notation
    for _magnitude in new_event.magnitudes:
        _magnitude.magnitude_type = _nortoevmag(_magnitude.magnitude_type)
    # Set the useful things like preferred magnitude and preferred origin
    new_event.preferred_origin_id = str(new_event.origins[0].resource_id)
    if len(new_event.magnitudes) > 1:
        try:
            # Select moment first, then local, then
            mag_filter = ['MW', 'Mw', 'ML', 'Ml', 'MB', 'Mb',
                          'MS', 'Ms', 'Mc', 'MC']
            _magnitudes = [(m.magnitude_type, m.resource_id)
                           for m in new_event.magnitudes]
            preferred_magnitude = sorted(_magnitudes,
                                         key=lambda x: mag_filter.index(x[0]))
            new_event.preferred_magnitude_id = str(preferred_magnitude[0][1])
        except ValueError:
            # If there is a magnitude not specified in filter
            new_event.preferred_magnitude_id =\
                str(new_event.magnitudes[0].resource_id)
    elif len(new_event.magnitudes) == 1:
        new_event.preferred_magnitude_id =\
            str(new_event.magnitudes[0].resource_id)
    return new_event


def read_event(sfile):
    """
    Read all information from a Nordic formatted s-file.
    Maps to readpicks and readheader to read origin and pick information,
    outputs an obspy.core.event.Event.

    :type sfile: str
    :param sfile: Nordic formatted file to open and read from.

    :returns: event
    :rtype: obspy.core.event.Event.
    """
    event = readpicks(sfile)
    return event


def read_select(select_file):
    """
    Read a catalog of events from a Nordic formatted select file.
    Generates a series of temporary files for each event in the select file.

    :type select_file: str
    :param select_file: Nordic formatted select.out file to open

    :return: catalog of events
    :rtype: obspy.core.event.Catalog
    """
    from obspy.core.event import Catalog
    from tempfile import NamedTemporaryFile
    import os

    catalog = Catalog()
    event_str = []
    f = open(select_file, 'r')
    for line in f:
        if len(line.rstrip()) > 0:
            event_str.append(line)
        elif len(event_str) > 0:
            # Write to a temporary file then read from it
            tmp_sfile = NamedTemporaryFile(mode='w', delete=False)
            for event_line in event_str:
                tmp_sfile.write(event_line)
            tmp_sfile.close()
            catalog += read_event(tmp_sfile.name)
            os.remove(tmp_sfile.name)
            event_str = []
    return catalog


def readpicks(sfile):
    """
    Read all pick information from the s-file to an obspy.event.Catalog type.

    .. note:: This was changed for version 0.1.0 from using the inbuilt \
    PICK class.

    :type sfile: str
    :param sfile: Path to sfile

    :return: obspy.core.event.Event

    .. warning:: Currently finalweight is unsupported, nor is velocity, \
    or angle of incidence.  This is because obspy.event stores slowness \
    in s/deg and takeoff angle, which would require computation from the \
    values stored in seisan.  Multiple weights are also not supported in \
    Obspy.event.

    .. rubric:: Example

    >>> event = readpicks('eqcorrscan/tests/test_data/REA/TEST_/' +
    ...                   '01-0411-15L.S201309')
    >>> print(event.origins[0].time)
    2013-09-01T04:11:15.700000Z
    >>> print(event.picks[0].time)
    2013-09-01T04:11:17.240000Z
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
                               int(line[23:28].split('.')[1]) * 10000)
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
    Extract the waveform filename from the s-file.
    Returns a list of waveform names found in the s-file as multiples can \
    be present.

    :type sfile: str
    :param sfile: Path to the sfile

    :returns: List of strings of wave paths
    :rtype: list

    >>> readwavename('eqcorrscan/tests/test_data/REA/TEST_/' +
    ...              '01-0411-15L.S201309')
    ['2013-09-01-0410-35.DFDPC_024_00']
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
    Generate an empty s-file with a populated header for a given waveform.

    :type wavefile: str
    :param wavefile: Wavefile to associate with this S-file, the timing of \
        the S-file will be taken from this file if evtime is not set.
    :type evtype: str
    :param evtype: Event type letter code, e.g. L, R, D
    :type userID: str
    :param userID: 4-character SEISAN USER ID
    :type outdir: str
    :param outdir: Location to write S-file
    :type overwrite: bool
    :param overwrite: Overwrite an existing S-file, default=False
    :type evtime: obspy.core.utcdatetime.UTCDateTime
    :param evtime: If given this will set the timing of the S-file

    :returns: str, S-file name

    >>> from eqcorrscan.utils.sfile_util import readwavename
    >>> import os
    >>> wavefile = os.path.join('eqcorrscan', 'tests', 'test_data', 'WAV',
    ...                         'TEST_', '2013-09-01-0410-35.DFDPC_024_00')
    >>> sfile = blanksfile(wavefile, 'L', 'TEST',
    ...                    '.', overwrite=True)
    Written s-file: ./01-0410-35L.S201309
    >>> readwavename(sfile)
    ['2013-09-01-0410-35.DFDPC_024_00']
    """

    from obspy import read as obsread
    import os
    import datetime

    if not evtime:
        try:
            st = obsread(wavefile)
            evtime = st[0].stats.starttime
        except:
            raise IOError('Wavefile: ' + wavefile +
                          ' is invalid, try again with real data.')
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

    # Generate s-file name in the format dd-hhmm-ss[L,R,D].Syyyymm
    sfile = outdir + '/' + str(evtime.day).zfill(2) + '-' +\
        str(evtime.hour).zfill(2) +\
        str(evtime.minute).zfill(2) + '-' +\
        str(evtime.second).zfill(2) + evtype + '.S' +\
        str(evtime.year) +\
        str(evtime.month).zfill(2)
    # Check is sfile exists
    if os.path.isfile(sfile) and not overwrite:
        print('Desired sfile: ' + sfile + ' exists, will not overwrite')
        for i in range(1, 10):
            sfile = outdir + '/' + str(evtime.day).zfill(2) + '-' +\
                str(evtime.hour).zfill(2) +\
                str(evtime.minute).zfill(2) + '-' +\
                str(evtime.second + i).zfill(2) + evtype + '.S' +\
                str(evtime.year) +\
                str(evtime.month).zfill(2)
            if not os.path.isfile(sfile):
                break
        else:
            msg = 'Tried generated files up to 20s in advance and found ' +\
                'all exist, you need to clean your stuff up!'
            raise IOError(msg)
        # sys.exit()
    f = open(sfile, 'w')
    # Write line 1 of s-file
    f.write(str(' ' + str(evtime.year) + ' ' +
                str(evtime.month).rjust(2) +
                str(evtime.day).rjust(2) + ' ' +
                str(evtime.hour).rjust(2) +
                str(evtime.minute).rjust(2) + ' ' +
                str(float(evtime.second)).rjust(4) + ' ' +
                evtype + '1'.rjust(58) + '\n'))
    # Write line 2 of s-file
    f.write(str(' ACTION:ARG ' + str(datetime.datetime.now().year)[2:4] + '-' +
                str(datetime.datetime.now().month).zfill(2) + '-' +
                str(datetime.datetime.now().day).zfill(2) + ' ' +
                str(datetime.datetime.now().hour).zfill(2) + ':' +
                str(datetime.datetime.now().minute).zfill(2) + ' OP:' +
                userID.ljust(4) + ' STATUS:' + 'ID:'.rjust(18) +
                str(evtime.year) +
                str(evtime.month).zfill(2) +
                str(evtime.day).zfill(2) +
                str(evtime.hour).zfill(2) +
                str(evtime.minute).zfill(2) +
                str(evtime.second).zfill(2) +
                'I'.rjust(6) + '\n'))
    # Write line 3 of s-file
    write_wavfile = wavefile.split(os.sep)[-1]
    f.write(str(' ' + write_wavfile + '6'.rjust(79 - len(write_wavfile)) +
                '\n'))
    # Write final line of s-file
    f.write(str(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU' +
                ' VELO AIN AR TRES W  DIS CAZ7\n'))
    f.close()
    print('Written s-file: ' + sfile)
    return sfile


def eventtosfile(event, userID, evtype, outdir, wavefiles, explosion=False,
                 overwrite=False):
    """
    Write an obspy.event to a nordic formatted s-file.

    :type event: obspy.core.event.Event
    :param event: A single obspy event
    :type userID: str
    :param userID: Up to 4 character user ID
    :type evtype: str
    :param evtype: Single character string to describe the event, either L, R \
        or D.
    :type outdir: str
    :param outdir: Path to directory to write to
    :type wavefiles: list
    :param wavefiles: Waveforms to associate the sfile with
    :type explosion: bool
    :param explosion: Note if the event is an explosion, will be marked by an \
        E.
    :type overwrite: bool
    :param overwrite: force to overwrite old files, defaults to False

    :returns: str: name of sfile written

    .. note:: Seisan can find waveforms either by their relative or absolute \
        path, or by looking for the file recursively in directories within \
        the WAV directory in your seisan install.  Because all lines need to \
        be less than 79 ccharacterslong (fortran hangover) in the s-files, \
        you will need to determine whether the full-path is okay or not.

    >>> import obspy
    >>> # Note that this example shows how to download from GeoNet which
    >>> # doesn't have full fdsn capability.
    >>> if int(obspy.__version__.split('.')[0]) >= 1:
    ...    from obspy.clients.fdsn import Client
    ...    from obspy import read_events
    ... else:
    ...    from obspy.fdsn import Client
    ...    from obspy import readEvents as read_events
    >>> client = Client('GEONET')
    >>> data_stream = client._download('http://quakeml.geonet.org.nz/' +
    ...                                'quakeml/1.2/2016p008122')
    >>> close = data_stream.seek(0, 0)
    >>> catalog = read_events(data_stream, format="quakeml")
    >>> data_stream.close()
    >>> eventtosfile(catalog[0], 'TEST', 'R', '.', ['DUMMY'], overwrite=True)
    '04-0007-55R.S201601'
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
    if isinstance(wavefiles, string_types):
        wavefiles = [str(wavefiles)]
    elif isinstance(wavefiles, list):
        wavefiles = wavefiles
    else:
        print(type(wavefiles))
        raise IOError(wavefiles + ' is neither string or list')
    # Determine name from origin time
    try:
        evtime = event.origins[0].time
    except IndexError:
        msg = 'list index out of range: Need at least one origin with at ' +\
            'least an origin time'
        raise IndexError(msg)
    if not evtime:
        msg = 'event has an origin, but time is not populated.  ' +\
            'This is required!'
        raise ValueError(msg)
    # Attempt to cope with possible pre-existing files
    range_list = []
    for i in range(30):  # Look +/- 30 seconds around origin time
        range_list.append(i)
        range_list.append(-1 * i)
    range_list = range_list[1:]
    for add_secs in range_list:
        sfilename = (evtime + add_secs).datetime.strftime('%d-%H%M-%S') +\
            evtype[0] + '.S' + (evtime + add_secs).datetime.strftime('%Y%m')
        if not os.path.isfile(outdir + os.sep + sfilename):
            sfile = open(outdir + os.sep + sfilename, 'w')
            break
        elif overwrite:
            sfile = open(outdir + os.sep + sfilename, 'w')
            break
    else:
        raise IOError(outdir + os.sep + sfilename +
                      ' already exists, will not overwrite')
    # Write the header info.
    if event.origins[0].latitude:
        if event.origins[0].latitude not in [float('NaN'), 999]:
            lat = '{0:.3f}'.format(event.origins[0].latitude)
        else:
            lat = ''
    else:
        lat = ''
    if event.origins[0].longitude:
        if event.origins[0].longitude not in [float('NaN'), 999]:
            lon = '{0:.3f}'.format(event.origins[0].longitude)
        else:
            lon = ''
    else:
        lon = ''
    if event.origins[0].depth:
        if event.origins[0].depth not in [float('NaN'), 999]:
            depth = '{0:.1f}'.format(event.origins[0].depth / 1000)
        else:
            depth = ''
    else:
        depth = ''
    if event.creation_info:
        try:
            agency = event.creation_info.get('agency_id')
            # If there is creation_info this may not raise an error annoyingly
            if agency is None:
                agency = ''
        except AttributeError:
            agency = ''
    else:
        agency = ''
    if len(agency) > 3:
        agency = agency[0:3]
    # Cope with differences in event uncertainty naming
    if event.origins[0].time_errors:
        try:
            timerms = '{0:.1f}'.format(event.origins[0].
                                       time_errors.Time_Residual_RMS)
        except AttributeError:
            timerms = '0.0'
    else:
        timerms = '0.0'
    try:
        mag_1 = '{0:.1f}'.format(event.magnitudes[0].mag) or ''
        mag_1_type = _evmagtonor(event.magnitudes[0].magnitude_type) or ''
        if event.magnitudes[0].creation_info:
            mag_1_agency = event.magnitudes[0].creation_info.agency_id or ''
        else:
            mag_1_agency = ''
    except IndexError:
        mag_1 = ''
        mag_1_type = ''
        mag_1_agency = ''
    try:
        mag_2 = '{0:.1f}'.format(event.magnitudes[1].mag) or ''
        mag_2_type = _evmagtonor(event.magnitudes[1].magnitude_type) or ''
        if event.magnitudes[1].creation_info:
            mag_2_agency = event.magnitudes[1].creation_info.agency_id or ''
        else:
            mag_2_agency = ''
    except IndexError:
        mag_2 = ''
        mag_2_type = ''
        mag_2_agency = ''
    try:
        mag_3 = '{0:.1f}'.format(event.magnitudes[2].mag) or ''
        mag_3_type = _evmagtonor(event.magnitudes[2].magnitude_type) or ''
        if event.magnitudes[2].creation_info:
            mag_3_agency = event.magnitudes[2].creation_info.agency_id or ''
        else:
            mag_3_agency = ''
    except IndexError:
        mag_3 = ''
        mag_3_type = ''
        mag_3_agency = ''
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
                userID.ljust(4) + ' STATUS:' + 'ID:'.rjust(18) +
                str(evtime.year) +
                str(evtime.month).zfill(2) +
                str(evtime.day).zfill(2) +
                str(evtime.hour).zfill(2) +
                str(evtime.minute).zfill(2) +
                str(evtime.second).zfill(2) +
                'I'.rjust(6) + '\n')
    # Write line 3 of s-file
    for wavefile in wavefiles:
        sfile.write(' ' + wavefile + '6'.rjust(79 - len(wavefile)) + '\n')
    # Write final line of s-file
    sfile.write(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU' +
                ' VELO AIN AR TRES W  DIS CAZ7\n')
    sfile.close()
    # Now call the populatesfile function
    if len(event.picks) > 0:
        populatesfile(outdir + '/' + sfilename, event)
    return str(sfilename)


def populatesfile(sfile, event):
    """
    Populate a blank nordic format S-file with pick information.
    arguments required are the filename of the blank s-file and the picks \
    where picks is a dictionary of picks including station, channel, \
    impulsivity, phase, weight, polarity, time, coda, amplitude, peri, \
    azimuth, velocity, SNR, azimuth residual, Time-residual, final weight, \
    epicentral distance & azimuth from event to station.

    This is a full pick line information from the seisan manual, P. 341

    :type sfile: str
    :param sfile: Path to S-file to populate, must have a header already
    :type event: obspy.core.event.Catalog
    :param picks: A single event to be written to a single S-file.

    >>> from eqcorrscan.utils.sfile_util import blanksfile, readpicks
    >>> sfile = blanksfile('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...                    '2013-09-01-0410-35.DFDPC_024_00', 'L', 'TEST',
    ...                    '.', overwrite=True)
    Written s-file: ./01-0410-35L.S201309
    >>> # Poor example, but we need an event, so we will use one we know is
    >>> # associated with the event...
    >>> event = readpicks('eqcorrscan/tests/test_data/REA/TEST_/' +
    ...                   '01-0411-15L.S201309')
    >>> populatesfile(sfile, event)
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


def nordpick(event):
    """
    Format information from an obspy.event class to nordic string format.

    :type event: obspy.core.event.Event
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
        if not pick.waveform_id:
            msg = 'No waveform id for pick, skipping'
            warnings.warn(msg)
            continue
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
                    distance = str(_int_conv(distance))
                elif 10.0 < distance < 100.0:
                    distance = _str_conv(round(distance, 1), 1)
                elif distance < 10.0:
                    distance = _str_conv(round(distance, 2), 2)
                else:
                    distance = _str_conv(distance, False)
            else:
                distance = ' '
            # Extract CAZ
            if arrival.azimuth:
                CAZ = int(arrival.azimuth)
            else:
                CAZ = ' '
        else:
            CAZ = ' '
            distance = ' '
            timeres = ' '
            azimuthres = ' '
            azimuth = ' '
            weight = 0
        if not pick.phase_hint:
            # Cope with some authorities not providing phase hints :(
            phase_hint = ' '
        else:
            phase_hint = pick.phase_hint
        # Extract amplitude: note there can be multiple amplitudes, but they
        # should be associated with different picks.
        amplitude = [amplitude for amplitude in event.amplitudes
                     if amplitude.pick_id == pick.resource_id]
        if len(amplitude) > 0:
            if len(amplitude) > 1:
                msg = 'Nordic files need one pick for each amplitude, ' + \
                      'using the first amplitude only'
                warnings.warn(msg)
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
                    peri_round = False
                # Extract amplitude and convert units
                if amplitude.generic_amplitude:
                    amp = amplitude.generic_amplitude
                    if amplitude.unit in ['m', 'm/s', 'm/(s*s)', 'm*s']:
                        amp *= 10**9
                    # Otherwise we will assume that the amplitude is in counts
                else:
                    amp = np.nan
                coda = ' '
                if amplitude.magnitude_hint == 'Ml':
                    phase_hint = 'IAML'
                    impulsivity = ' '
            else:
                coda = int(amplitude.generic_amplitude)
                peri = ' '
                peri_round = False
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
        channel_code = pick.waveform_id.channel_code or '   '
        pick_strings.append(' ' + pick.waveform_id.station_code.ljust(5) +
                            channel_code[0] + channel_code[-1] +
                            ' ' + impulsivity + phase_hint.ljust(4) +
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
                            _str_conv(timeres, rounded=2).rjust(5)[0:5] +
                            _str_conv(' ').rjust(2) +
                            distance.rjust(5) +
                            _str_conv(CAZ).rjust(4) + ' ')
        # Note that currently finalweight is unsupported, nor is velocity, or
        # angle of incidence.  This is because obspy.event stores slowness in
        # s/deg and takeoff angle, which would require computation from the
        # values stored in seisan.  Multiple weights are also not supported in
        # Obspy.event
    return pick_strings


def stationtoseisan(station):
    """
    Convert obspy inventory to string formatted for Seisan STATION0.HYP file.

    :type station: obspy.core.inventory.station.Station
    :param station: Inventory containing a single station.

    :returns: str

    .. note:: Only works to the low-precision level at the moment (see seisan \
        manual for explanation).
    """

    if station.latitude < 0:
        lat_str = 'S'
    else:
        lat_str = 'N'
    if station.longitude < 0:  # Stored in =/- 180, not 0-360
        lon_str = 'W'
    else:
        lon_str = 'E'
    if len(station.code) > 4:
        sta_str = station.code[0:4]
    else:
        sta_str = station.code.ljust(4)
    if len(station.channels) > 0:
        depth = station.channels[0].depth
    else:
        msg = 'No depth found in station.channels, have you set the level ' +\
              'of stationXML download to channel if using obspy.get_stations?'
        raise IOError(msg)
    elev = str(int(round(station.elevation - depth))).rjust(4)
    # lat and long are written in STATION0.HYP in deg,decimal mins
    lat = abs(station.latitude)
    lat_degree = int(lat)
    lat_decimal_minute = (lat - lat_degree) * 60
    lon = abs(station.longitude)
    lon_degree = int(lon)
    lon_decimal_minute = (lon - lon_degree) * 60
    lat = ''.join([str(int(abs(lat_degree))),
                   '{0:.2f}'.format(lat_decimal_minute).rjust(5)])
    lon = ''.join([str(int(abs(lon_degree))),
                   '{0:.2f}'.format(lon_decimal_minute).rjust(5)])
    station_str = ''.join(['  ', sta_str, lat, lat_str, lon, lon_str, elev])
    return station_str


if __name__ == "__main__":
    import doctest
    doctest.testmod()
