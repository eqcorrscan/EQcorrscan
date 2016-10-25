"""
Functions to aid downloading of GeoNet events into obspy catalog objects.

:copyright:
    Calum Chamberlain
:licence:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""


def get_geonet_ids(minlat, maxlat, minlon, maxlon, mindepth=None,
                   maxdepth=None, minmag=None, maxmag=None,
                   startdate=None, enddate=None):
    """Generate quakesearch URL query and extract publicID from returned csv.

    :type minlat: float
    :param minlat: Southern edge of bounding box.
    :type maxlat: float
    :param maxlat: Northern edge of bounding box.
    :type minlon: float
    :param minlon: Western edge of bounding box.
    :type maxlon: float
    :param maxlon: Eastern edge of bounding box.
    :type mindepth: float
    :param mindepth: Minimum depth for events (depth is positive down).
    :type maxdepth: float
    :param maxdepth: Maximum depth for events (depth is positive down).
    :type minmag: float
    :param minmag: Minimum magnitude to extract.
    :type maxmag: float
    :param maxmag: Maximum magnitude to extract.
    :type startdate: obspy.core.UTCDateTime
    :param startdate: Start to search.
    :type enddate: obspy.core.UTCDateTime
    :param enddate: End of search.

    :returns: list of str of event ids
    """
    import csv
    import sys
    import io
    if sys.version_info.major == 2:
        from urllib2 import urlopen
    else:
        from urllib.request import urlopen

    base_url = "http://quakesearch.geonet.org.nz/csv?"
    bbox_url = "bbox=" + ','.join([str(minlon), str(minlat),
                                   str(maxlon), str(maxlat)])
    url = base_url + bbox_url
    if mindepth:
        url += "&mindepth=" + str(mindepth)
    if maxdepth:
        url += "&maxdepth=" + str(maxdepth)
    if minmag:
        url += "&minmag=" + str(minmag)
    if maxmag:
        url += "&maxmag=" + str(maxmag)
    if startdate:
        startdate_url = "&startdate=" + startdate.strftime('%Y-%m-%dT%H:%M:%S')
        url += startdate_url
    if enddate:
        enddate_url = "&enddate=" + enddate.strftime('%Y-%m-%dT%H:%M:%S')
        url += enddate_url
    print("Downloading info from:")
    print(url)
    response = urlopen(url)
    if sys.version_info.major == 3:
        quake_search = csv.reader(io.TextIOWrapper(response))
    else:
        quake_search = csv.reader(response)

    header = next(quake_search)
    # Usually publicID is the first column, error if not true
    if not header[0] == 'publicid':
        raise IOError('Unexpected format, first column is not publicid')
    event_ids = [row[0] for row in quake_search]
    return event_ids


def _get_geonet_pubids(publicids, parallel=True):
    """
    Get GeoNet events while they haven't included get_events in fdsn.

    :type publicids: list
    :param publicids: List of public id numbers for events wanted.

    :returns: Catalog of events
    :rtype: obspy.core.event.Catalog
    """
    import obspy
    if int(obspy.__version__.split('.')[0]) > 0:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client
    from obspy.core.event import Catalog
    from multiprocessing import Pool, cpu_count

    client = Client('GEONET')
    catalog = Catalog()
    # Multi-process this bad-boy
    if not parallel:
        for publicid in publicids:
            catalog += _inner_get_event(publicid=publicid, client=client)
    else:
        pool = Pool(processes=cpu_count())
        results = [pool.apply_async(_inner_get_event, args=(publicid, client))
                   for publicid in publicids]
        pool.close()
        cat_list = [p.get() for p in results]
        pool.join()
        for ev in cat_list:
            catalog += ev
    return catalog


def _inner_get_event(publicid, client):
    """
    Inner loop for parallel processing

    :type publicid: str
    :param publicid: GeoNet public ID
    :return: catalog
    """
    import warnings
    from obspy.clients.fdsn.header import FDSNException
    from obspy import read_events
    try:
        data_stream = client._download('http://quakeml.geonet.org.nz/' +
                                       'quakeml/1.2/' + publicid)
        data_stream.seek(0, 0)
        catalog = read_events(data_stream, format="quakeml")
        data_stream.close()
    except FDSNException:
        warnings.warn('Unable to download event: ' + publicid)
    return catalog


def get_geonet_events(minlat, maxlat, minlon, maxlon, mindepth=None,
                      maxdepth=None, minmag=None, maxmag=None,
                      startdate=None, enddate=None):
    """Generate quakesearch URL query and extract publicID from returned csv.

    :type minlat: float
    :param minlat: Southern edge of bounding box.
    :type maxlat: float
    :param maxlat: Northern edge of bounding box.
    :type minlon: float
    :param minlon: Western edge of bounding box.
    :type maxlon: float
    :param maxlon: Eastern edge of bounding box.
    :type mindepth: float
    :param mindepth: Minimum depth for events (depth is positive down).
    :type maxdepth: float
    :param maxdepth: Maximum depth for events (depth is positive down).
    :type minmag: float
    :param minmag: Minimum magnitude to extract.
    :type maxmag: float
    :param maxmag: Maximum magnitude to extract.
    :type startdate: obspy.core.UTCDateTime
    :param startdate: Start to search.
    :type enddate: obspy.core.UTCDateTime
    :param enddate: End of search.

    :returns: catalog of events
    :rtype: obspy.core.event.Catalog
    """
    pubids = get_geonet_ids(minlat=minlat, maxlat=maxlat, minlon=minlon,
                            maxlon=maxlon, mindepth=mindepth,
                            maxdepth=maxdepth, minmag=minmag,
                            maxmag=maxmag, startdate=startdate,
                            enddate=enddate)
    catalog = _get_geonet_pubids(pubids)
    return catalog
