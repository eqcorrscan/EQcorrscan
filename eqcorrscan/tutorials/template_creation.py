"""
Simple tutorial detailing how to generate a series of templates from catalog\
data available online.
"""


def mktemplates(network_code='GEONET',
                publicIDs=['2016p008122', '2016p008353', '2016p008155',
                           '2016p008194'], plot=True):
    """Functional wrapper to make templates"""

    from collections import Counter
    from eqcorrscan.core import template_gen

    # This import section copes with namespace changes between obspy versions
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
        from obspy import read_events
    else:
        from obspy.fdsn import Client
        from obspy import readEvents as read_events
    from obspy.core.event import Catalog

    # We want to download some QuakeML files from the New Zealand GeoNet
    # network, GeoNet currently doesn't support FDSN event queries, so we
    # have to work around to download quakeml from their quakeml.geonet site.

    client = Client(network_code)
    # We want to download a few events from an earthquake sequence, these are
    # identified by publiID numbers, given as arguments

    catalog = Catalog()
    for publicID in publicIDs:
        if network_code == 'GEONET':
            data_stream = client._download('http://quakeml.geonet.org.nz/' +
                                           'quakeml/1.2/' + publicID)
            data_stream.seek(0, 0)
            catalog += read_events(data_stream, format="quakeml")
            data_stream.close()
        else:
            catalog += client.get_events(eventid=publicID,
                                         includearrivals=True)

    # Lets plot the catalog to see what we have
    if plot:
        catalog.plot(projection='local', resolution='h')

    # We don't need all the picks, lets take the information from the
    # five most used stations - note that this is done to reduce computational
    # costs.
    all_picks = []
    for event in catalog:
        all_picks += [(pick.waveform_id.station_code,
                       pick.waveform_id.channel_code) for pick in event.picks]
    # Python 3.x and python 2.7 do not sort in the same way, this is a cludge
    # to work around that...
    counted = Counter(all_picks).most_common()
    # Going to take an initial set that all have atleast 1 less pick than the
    # highest pick-count...
    all_picks = []
    for i in range(counted[0][1]):
        highest = [item[0] for item in counted if item[1] >= counted[0][1] - i]
        # Sort them by alphabetical order in station
        highest = sorted(highest, key=lambda tup: tup[0])
        all_picks += highest
        if len(all_picks) > 5:
            all_picks = all_picks[0:5]
            break

    for event in catalog:
        if len(event.picks) == 0:
            raise IOError('No picks found')
        event.picks = [pick for pick in event.picks
                       if (pick.waveform_id.station_code,
                           pick.waveform_id.channel_code) in all_picks]

    # Now we can generate the templates
    templates = template_gen.from_client(catalog=catalog,
                                         client_id=network_code,
                                         lowcut=2.0, highcut=9.0,
                                         samp_rate=20.0, filt_order=4,
                                         length=3.0, prepick=0.15,
                                         swin='all', debug=0, plot=plot)

    # We now have a series of templates! Using Obspys Stream.write() method we
    # can save these to disk for later use.  We will do that now for use in the
    # following tutorials.
    for i, template in enumerate(templates):
        template.write('tutorial_template_' + str(i) + '.ms', format='MSEED')
        # Note that this will warn you about data types.  As we don't care
        # at the moment, whatever obspy chooses is fine.
    return


if __name__ == '__main__':
    """Wrapper for template creation"""
    import sys
    import warnings
    if not len(sys.argv) > 1:
        warnings.warn('Needs a network ID followed by a list of event IDs, ' +
                      'will run the test case instead')
        mktemplates()
    else:
        net_code = sys.argv[1]
        idlist = list(sys.argv)[2:]
        print(idlist)
        mktemplates(net_code, idlist)
