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
import os
import logging
import glob
import fnmatch

from obspy import read, UTCDateTime, Stream
from obspy.clients.fdsn.header import FDSNException
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.filesystem.sds import Client as SDSClient


Logger = logging.getLogger(__name__)


def read_data(archive, arc_type, day, stachans, length=86400):
    """
    Function to read the appropriate data from an archive for a day.

    :type archive: str
    :param archive:
        The archive source -  if arc_type is FDSN then this can be either a url
        or a known obspy client. If arc_type is day_vols, then this is the
        path to the top directory.
    :type arc_type: str
    :param arc_type: The type of archive, can be: FDSN, day_volumes
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

    >>> from obspy import UTCDateTime
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> stachans = [('JCNB', 'SP1')]
    >>> st = read_data('NCEDC', 'FDSN', t1, stachans)
    >>> print(st)
    1 Trace(s) in Stream:
    BP.JCNB.40.SP1 | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.\
950000Z | 20.0 Hz, 1728000 samples

    .. rubric:: Example, missing data

    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> stachans = [('JCNB', 'SP1'), ('GCSZ', 'HHZ')]
    >>> st = read_data('NCEDC', 'FDSN', t1, stachans)
    >>> print(st)
    1 Trace(s) in Stream:
    BP.JCNB.40.SP1 | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.\
950000Z | 20.0 Hz, 1728000 samples


    .. rubric:: Example, local day-volumes

    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> t1 = UTCDateTime(2012, 3, 26)
    >>> stachans = [('WHYM', 'SHZ'), ('EORO', 'SHZ')]
    >>> st = read_data(TEST_PATH + '/day_vols', 'day_vols',
    ...                t1, stachans)
    >>> print(st)
    2 Trace(s) in Stream:
    AF.WHYM..SHZ | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.000000Z \
| 1.0 Hz, 86400 samples
    AF.EORO..SHZ | 2012-03-26T00:00:00.000000Z - 2012-03-26T23:59:59.000000Z \
| 1.0 Hz, 86400 samples
    """
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
        # Allow matching of station-channel-tuples with wildcards
        stachan_in_station_map = False
        for av_station in available_stations_map:
            if fnmatch.fnmatch(av_station[0], station[0])\
                    and fnmatch.fnmatch(av_station[1], station[1]):
                stachan_in_station_map = True
                break
        if not stachan_in_station_map:
            msg = ' '.join([station[0], station_map[1], 'is not available for',
                            day.strftime('%Y/%m/%d')])
            Logger.warning(msg)
            continue
        elif arc_type.upper() == "FDSN":
            client = FDSNClient(archive)
            try:
                st += client.get_waveforms(
                    network='*', station=station_map[0], location='*',
                    channel=station_map[1], starttime=UTCDateTime(day),
                    endtime=UTCDateTime(day) + length)
            except FDSNException:
                Logger.warning('No data on server despite station being ' +
                               'available...')
                continue
        elif arc_type.upper() == "SDS":
            client = SDSClient(archive)
            st += client.get_waveforms(
                    network='*', station=station_map[0], location='*',
                    channel=station_map[1], starttime=UTCDateTime(day),
                    endtime=UTCDateTime(day) + length)
        elif arc_type.lower() == 'day_vols':
            wavfiles = _get_station_file(os.path.join(
                archive, day.strftime('Y%Y' + os.sep + 'R%j.01')),
                station_map[0], station_map[1])
            for wavfile in wavfiles:
                st += read(wavfile, starttime=day, endtime=day + length)
    st = Stream(st)
    return st


def _get_station_file(path_name, station, channel):
    """
    Helper function to find the correct file.

    :type path_name: str
    :param path_name: Path to files to check.
    :type station: str
    :type channel: str

    :returns: list of filenames, str
    """
    wavfiles = glob.glob(path_name + os.sep + '*')

    out_files = [_check_data(wavfile, station, channel)
                 for wavfile in wavfiles]
    out_files = list(set(out_files))
    return out_files


def _check_data(wavfile, station, channel):
    """
    Inner loop for parallel checks.

    :type wavfile: str
    :param wavfile: Wavefile path name to look in.
    :type station: str
    :param station: Channel name to check for
    :type channel: str
    :param channel: Channel name to check for
    """
    Logger.debug('Checking ' + wavfile)
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


    """
    available_stations = []
    if arc_type.lower() == 'day_vols':
        wavefiles = glob.glob(os.path.join(archive, day.strftime('Y%Y'),
                                           day.strftime('R%j.01'), '*'))
        for wavefile in wavefiles:
            header = read(wavefile, headonly=True)
            available_stations.append((header[0].stats.station,
                                       header[0].stats.channel))
    elif arc_type.lower() == 'fdsn':
        client = FDSNClient(archive)
        inventory = client.get_stations(starttime=UTCDateTime(day),
                                        endtime=UTCDateTime(day) + 86400,
                                        level='channel')
        for network in inventory:
            for station in network:
                for channel in station:
                    available_stations.append((station.code,
                                               channel.code))
    elif arc_type.upper() == "SDS":
        client = SDSClient(archive)
        nslc = client.get_all_nslc(sds_type=None, datetime=UTCDateTime(day))
        for item in nslc:
            available_stations.append((item[1], item[3]))
    return available_stations


if __name__ == '__main__':
    import doctest
    doctest.testmod()
