"""
Functions to read log-files for seismic data to determine whether there are \
timing issues present.  Designed to be used with the EQcorrscan package and \
to flag data that has more than a threshold timing issue.

.. note:: Currently only written to read RefTek rt130 log-files, and will not \
    read all parameters - only for use when checking logs during \
    cross-correlation.  For full log-file exploration, Passcal tools: logpeek \
    is useful.

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

import datetime as dt
import re
import os
import io
import warnings
import glob
import sys


def rt_time_log(logfile, startdate):
    """
    Open and read reftek raw log-file.

    Function to open and read a log-file as written by a RefTek RT130
    datalogger. The information within is then scanned for timing errors
    above the threshold.

    :type logfile: str
    :param logfile: The logfile to look in
    :type startdate: datetime.date
    :param startdate: The start of the file as a date - files contain timing \
        and the julian day, but not the year.

    :returns: List of tuple of (:class:`datetime.datetime`, float) as time \
        stamps and phase error.
    """
    if os.name == 'nt':
        f = io.open(logfile, 'rb')
    else:
        f = io.open(logfile, 'rb')
    phase_err = []
    lock = []
    # Extract all the phase errors
    for line_binary in f:
        try:
            line = line_binary.decode("utf8", "ignore")
        except UnicodeDecodeError:
            warnings.warn('Cannot decode line, skipping')
            continue
        if re.search("INTERNAL CLOCK PHASE ERROR", line):
            match = re.search("INTERNAL CLOCK PHASE ERROR", line)
            d_start = match.start() - 13
            phase_err.append((dt.datetime.strptime(str(startdate.year) +
                                                   ':' +
                                                   line[d_start:d_start + 12],
                                                   '%Y:%j:%H:%M:%S'),
                             float(line.rstrip().split()[-2]) *
                             0.000001))
        elif re.search("EXTERNAL CLOCK POWER IS TURNED OFF", line):
            match = re.search("EXTERNAL CLOCK POWER IS TURNED OFF", line)
            d_start = match.start() - 13
            lock.append((dt.datetime.strptime(str(startdate.year) +
                                              ':' + line[d_start:d_start + 12],
                                              '%Y:%j:%H:%M:%S'),
                        999))
    if len(phase_err) == 0 and len(lock) > 0:
        phase_err = lock
    f.close()
    return phase_err


def rt_location_log(logfile):
    """
    Extract location information from a RefTek raw log-file.

    Function to read a specific RefTek RT130 log-file and find all location
    information.

    :type logfile: str
    :param logfile: The logfile to look in

    :returns: list of tuples of lat, lon, elevation in decimal degrees and km.
    :rtype: list
    """
    if os.name == 'nt':
        f = open(logfile, 'rb')
    else:
        f = open(logfile, 'rb')
    locations = []
    for line_binary in f:
        try:
            line = line_binary.decode("utf8", "ignore")
        except UnicodeDecodeError:
            warnings.warn('Cannot decode line, skipping')
            print(line_binary)
            continue
        match = re.search("GPS: POSITION:", line)
        if match:
            # Line is of form:
            # jjj:hh:mm:ss GPS: POSITION: xDD:MM:SS.SS xDDD:MM:SS.SS xMMMMMMM
            loc = line[match.end() + 1:].rstrip().split(' ')
            lat_sign = loc[0][0]
            lat = loc[0][1:].split(':')
            lat = int(lat[0]) + (int(lat[1]) / 60.0) + (float(lat[2]) / 3600.0)
            if lat_sign == 'S':
                lat *= -1
            lon_sign = loc[1][0]
            lon = loc[1][1:].split(':')
            lon = int(lon[0]) + (int(lon[1]) / 60.0) + (float(lon[2]) / 3600.0)
            if lon_sign == 'W':
                lon *= -1
            elev_sign = loc[2][0]
            elev_unit = loc[2][-1]
            if not elev_unit == 'M':
                raise NotImplementedError('Elevation is not in M: unit=' +
                                          elev_unit)
            elev = int(loc[2][1:-1])
            if elev_sign == '-':
                elev *= -1
            # Convert to km
            elev /= 1000
            locations.append((lat, lon, elev))
    f.close()
    return locations


def flag_time_err(phase_err, time_thresh=0.02):
    """
    Find large time errors in list.

    Scan through a list of tuples of time stamps and phase errors
    and return a list of time stamps with timing errors above a threshold.

    .. note::
        This becomes important for networks cross-correlations, where
        if timing information is uncertain at one site, the relative arrival
        time (lag) will be incorrect, which will degrade the cross-correlation
        sum.

    :type phase_err: list
    :param phase_err: List of Tuple of float, datetime.datetime
    :type time_thresh: float
    :param time_thresh: Threshold to declare a timing error for

    :returns: List of :class:`datetime.datetime` when timing is questionable.
    """
    time_err = []
    for stamp in phase_err:
        if abs(stamp[1]) > time_thresh:
            time_err.append(stamp[0])
    return time_err


def check_all_logs(directory, time_thresh):
    """
    Check all the log-files in a directory tree for timing errors.

    :type directory: str
    :param directory: Directory to search within
    :type time_thresh: float
    :param time_thresh: Time threshold in seconds

    :returns: List of :class:`datetime.datetime` for which error timing is
        above threshold, e.g. times when data are questionable.
    :rtype: list
    """
    log_files = glob.glob(directory + '/*/0/000000000_00000000')
    print('I have ' + str(len(log_files)) + ' log files to scan')
    total_phase_errs = []
    for i, log_file in enumerate(log_files):
        startdate = dt.datetime.strptime(log_file.split('/')[-4][0:7],
                                         '%Y%j').date()
        total_phase_errs += rt_time_log(log_file, startdate)
        sys.stdout.write("\r" + str(float(i) / len(log_files) * 100) +
                         "% \r")
        sys.stdout.flush()
    time_errs = flag_time_err(total_phase_errs, time_thresh)
    time_errs.sort()
    return time_errs, total_phase_errs


if __name__ == "__main__":
    import doctest
    doctest.testmod()
