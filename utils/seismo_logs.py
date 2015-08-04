#!/usr/bin/python
"""
Functions to read log-files for seismic data to acertain whether there are
timing issues present.  Designed to be used with the EQcorrscan package and
to flag data that has more than a threshold timing issue.

Currently only written to read RefTek rt130 log-files.

Written by Calum Chamberlain, VUW 2015

All rights reserved, this is distributed under the LGPL licence.
"""

def Read_RT_log(logfile, startdate):
    """
    Function to open and read a log-file as written by a RefTek RT130 datalogger.\
    The information within is then scanned for timing errors above the threshold.

    :type logfile: String
    :param logfile: The logfile to look in
    :type startdate: :class: datetime.date
    :param startdate: The start of the file as a date - files contain timing\
            and the julian day, but not the year.
    :type time_thresh: float
    :param time_thresh: Threshold to raise a flag for the data in seconds

    :returns: List of tuple of :class: datetime.datetime, float as time stamps\
            and phase error.
    """
    import datetime as dt
    f=open(logfile,'r')
    phase_err=[]
    # Extract all the phase errors
    for line in f:
        if line[13:39]=="INTERNAL CLOCK PHASE ERROR":
            phase_err.append((dt.datetime.strptime(str(startdate.year)+\
                                                            ':'+line[0:12],\
                                                            '%Y:%j:%H:%M:%S'),\
                             float(line[43:45].strip())*\
                             0.000001))
    return phase_err

def Flag_time_err(phase_err, time_thresh=0.02):
    """
    Fucntion to scan through a list of tuples of time stamps and phase errors
    and return a list of time stamps with timing errors above a threshold

    :type phase_err: List of Tuple of float, datetime.datetime
    :type time_thresh: float
    :param time_thresh: Threshold to declare a timing error for

    :returns: List of datetime.datetime
    """
    time_err=[]
    for stamp in phase_err:
        if abs(stamp[1])>time_thresh:
            time_err.append(stamp[0])
    return time_err

def check_all_logs(directory, time_thresh):
    """
    Function to check all the log-files in a directory tree for timing errors.

    :type directory: String
    :param directory: Directory to search within
    :type time_thresh: float
    :param time_thresh: Time threshold in seconds

    :returns: List of :class: datetime.datetime for which timing is above\
            threshold
    """
    import glob, sys
    import datetime as dt
    log_files=glob.glob(directory+'/*/*/0/000000000_00000000')
    print 'I have '+str(len(log_files))+' log files to scan'
    total_phase_errs=[]
    i=0
    for log_file in log_files:
        startdate=dt.datetime.strptime(log_file.split('/')[-4][0:7], '%Y%j').date()
        total_phase_errs+=Read_RT_log(log_file, startdate)
        sys.stdout.write("\r"+str(float(i)/len(log_files)*100)+"% \r")
        sys.stdout.flush()
        i+=1
    time_errs=Flag_time_err(total_phase_errs, time_thresh)
    time_errs.sort()
    return time_errs, total_phase_errs
