"""
Helper functions for reading data from archives for the EQcorrscan package.

.. note:: These functions are tools to aid simplification of general scripts, \
    they do not cover all use cases, however if you have a use case you want \
    to see here, then let the authors know, or implement it yourself and \
    contribute it back to the project.

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


def read_data(archive, arc_type, day, stachans, length=86400):
    """
    Function to read the appropriate data from an archive for a day.

    :type archive: str
    :param archive:
        The archive source - if arc_type is seishub, this should be a url,
        if the arc_type is FDSN then this can be either a url or a known obspy
        client.  If arc_type is day_vols, then this is the path to the top
        directory.
    :type arc_type: str
    :param arc_type: The type of archive, can be: seishub, FDSN, day_volves
    :type day: datetime.date
    :param day: Date to retrieve data for
    :type stachans: list
    :param stachans: List of tuples of Stations and channels to try and get,
        will not fail if stations are not available, but will warn.
    :type length: float
    :param length: Data length to extract in seconds, defaults to 1 day.

    :returns: Stream of data
    :rtype: obspy.core.stream.Stream

    .. note:: A note on arc_types, if arc_type is day_vols, then this will \
        look for directories labelled in the IRIS DMC conventions of \
        Yyyyy/Rjjj.01/... where yyyy is the year and jjj is the julian day. \
        Data within these files directories should be stored as day-long, \
        single-channel files.  This is not implemented in the fasted way \
        possible to allow for a more general situation.  If you require more \
        speed you will need to re-write this.

    .. rubric:: Example

    >>> from eqcorrscan.utils.archive_read import read_data
    >>> from obspy import UTCDateTime
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> stachans = [('FOZ', 'HHZ'), ('JCZ', 'HHZ')]
    >>> st = read_data('GEONET', 'FDSN', t1, stachans)
    >>> print(st)
    2 Trace(s) in Stream:
    NZ.FOZ.10.HHZ | 2012-03-25T23:59:57.018393Z - 2012-03-27T00:00:00.688393Z \
| 100.0 Hz, 8640368 samples
    NZ.JCZ.10.HHZ | 2012-03-25T23:59:57.348391Z - 2012-03-27T00:00:02.958391Z \
| 100.0 Hz, 8640562 samples


    .. rubric:: Example, missing data

    >>> from eqcorrscan.utils.archive_read import read_data
    >>> from obspy import UTCDateTime
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> stachans = [('FOZ', 'HHZ'), ('GCSZ', 'HHZ')]
    >>> st = read_data('GEONET', 'FDSN', t1, stachans)
    >>> print(st)
    1 Trace(s) in Stream:
    NZ.FOZ.10.HHZ | 2012-03-25T23:59:57.018393Z - 2012-03-27T00:00:00.688393Z \
| 100.0 Hz, 8640368 samples


    .. rubric:: Example, local day-volumes

    >>> from eqcorrscan.utils.archive_read import read_data
    >>> from obspy import UTCDateTime
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> stachans = [('WHYM', 'SHZ'), ('EORO', 'SHZ')]
    >>> st = read_data('eqcorrscan/tests/test_data/day_vols', 'day_vols',
    ...                t1, stachans)
    >>> print(st)
    2 Trace(s) in Stream:
    AF.WHYM..SHZ | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.000000Z \
| 1.0 Hz, 86400 samples
    AF.EORO..SHZ | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.000000Z \
| 1.0 Hz, 86400 samples
    """
    import obspy
    import os
    from obspy.clients.fdsn.header import FDSNException
    if arc_type.lower() == 'seishub':
        if int(obspy.__version__.split('.')[0]) >= 1:
            from obspy.clients.seishub import Client
        else:
            from obspy.seishub import Client
    else:
        if int(obspy.__version__.split('.')[0]) >= 1:
            from obspy.clients.fdsn import Client
        else:
            from obspy.fdsn import Client
    from obspy import read, UTCDateTime
    import warnings

    st = []
    available_stations = _check_available_data(archive, arc_type, day)
    for station in stachans:
        if len(station[1]) == 2:
            # Cope with two char channel naming in seisan
            station_map = (station[0], station[1][0] + '*' + station[1][1])
            available_stations_map = [(sta[0], sta[1][0] + '*' + sta[1][-1])
                                      for sta in available_stations]
        else:
            station_map = station
            available_stations_map = available_stations
        if station_map not in available_stations_map:
            msg = ' '.join([station[0], station_map[1], 'is not available for',
                            day.strftime('%Y/%m/%d')])
            warnings.warn(msg)
            continue
        if arc_type.lower() in ['seishub', 'fdsn']:
            client = Client(archive)
            try:
                st += client.get_waveforms(network='*', station=station_map[0],
                                           location='*',
                                           channel=station_map[1],
                                           starttime=UTCDateTime(day),
                                           endtime=UTCDateTime(day) + length)
            except FDSNException:
                warnings.warn('No data on server despite station being ' +
                              'available...')
                continue
        elif arc_type.lower() == 'day_vols':
            wavfiles = _get_station_file(os.path.join(archive,
                                                      day.strftime('Y%Y' +
                                                                   os.sep +
                                                                   'R%j.01')),
                                         station_map[0], station_map[1])
            for wavfile in wavfiles:
                st += read(wavfile, starttime=day, endtime=day + length)
    st = obspy.Stream(st)
    return st


def _get_station_file(path_name, station, channel, debug=0):
    """
    Helper function to find the correct file.

    :type path_name: str
    :param path_name: Path to files to check.
    :type station: str
    :type channel: str

    :returns: list of filenames, str
    """
    import glob
    import os
    from multiprocessing import Pool, cpu_count
    pool = Pool(processes=cpu_count())
    wavfiles = glob.glob(path_name + os.sep + '*')

    results = [pool.apply_async(_parallel_checking_loop,
                                args=(wavfile, station, channel, debug))
               for wavfile in wavfiles]
    pool.close()
    out_files = [p.get() for p in results]
    pool.join()
    out_files = list(set(out_files))
    return out_files


def _parallel_checking_loop(wavfile, station, channel, debug=0):
    """
    Inner loop for parallel checks.

    :type wavfile: str
    :param wavfile: Wavefile path name to look in.
    :type station: str
    :param station: Channel name to check for
    :type channel: str
    :param channel: Channel name to check for
    :type debug: int
    :param debug: Debug level, if > 1, will output what it it working on.
    """
    from obspy import read
    if debug > 1:
        print('Checking ' + wavfile)
    st = read(wavfile, headonly=True)
    for tr in st:
        if tr.stats.station == station and tr.stats.channel == channel:
            return wavfile


def _check_available_data(archive, arc_type, day):
    """
    Function to check what stations are available in the archive for a given \
    day.

    :type archive: str
    :param archive: The archive source
    :type arc_type: str
    :param arc_type: The type of archive, can be:
    :type day: datetime.date
    :param day: Date to retrieve data for

    :returns: list of tuples of (station, channel) as available.

    .. note:: Currently the seishub options are untested.

    """
    from obspy import read, UTCDateTime
    import glob
    import os

    available_stations = []
    if arc_type.lower() == 'day_vols':
        wavefiles = glob.glob(os.path.join(archive, day.strftime('Y%Y'),
                                           day.strftime('R%j.01'), '*'))
        for wavefile in wavefiles:
            header = read(wavefile, headonly=True)
            available_stations.append((header[0].stats.station,
                                       header[0].stats.channel))
    elif arc_type.lower() == 'seishub':
        from obspy.clients.seishub import Client
        client = Client(archive)
        st = client.get_previews(starttime=UTCDateTime(day),
                                 endtime=UTCDateTime(day) + 86400)
        for tr in st:
            available_stations.append((tr.stats.station, tr.stats.channel))
    elif arc_type.lower() == 'fdsn':
        from obspy.clients.fdsn import Client
        client = Client(archive)
        inventory = client.get_stations(starttime=UTCDateTime(day),
                                        endtime=UTCDateTime(day) + 86400,
                                        level='channel')
        for network in inventory:
            for station in network:
                for channel in station:
                    available_stations.append((station.code,
                                               channel.code))
    return available_stations


if __name__ == '__main__':
    import doctest
    doctest.testmod()
