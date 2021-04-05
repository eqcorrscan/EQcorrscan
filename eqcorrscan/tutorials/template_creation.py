"""
Simple tutorial detailing how to generate a series of templates from catalog\
data available online.
"""

import logging

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog

from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.core import template_gen

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")


def mktemplates(
        network_code='GEONET', plot=True, public_ids=None):
    """Functional wrapper to make templates"""
    public_ids = public_ids or [
        '2016p008122', '2016p008353', '2016p008155', '2016p008194']
    client = Client(network_code)
    # We want to download a few events from an earthquake sequence, these are
    # identified by publiID numbers, given as arguments

    catalog = Catalog()
    for public_id in public_ids:
        try:
            catalog += client.get_events(
                eventid=public_id, includearrivals=True)
        except TypeError:
            # Cope with some FDSN services not implementing includearrivals
            catalog += client.get_events(eventid=public_id)

    # Lets plot the catalog to see what we have
    if plot:
        catalog.plot(projection='local', resolution='h')

    # We don't need all the picks, lets take the information from the
    # five most used stations - note that this is done to reduce computational
    # costs.
    catalog = filter_picks(catalog, top_n_picks=5)
    # We only want the P picks in this example, but you can use others or all
    #  picks if you want.
    for event in catalog:
        for pick in event.picks:
            if pick.phase_hint == 'S':
                event.picks.remove(pick)

    # Now we can generate the templates
    templates = template_gen.template_gen(
        method='from_client', catalog=catalog, client_id=network_code,
        lowcut=2.0, highcut=9.0, samp_rate=20.0, filt_order=4, length=3.0,
        prepick=0.15, swin='all', process_len=3600, plot=plot)

    # We now have a series of templates! Using Obspy's Stream.write() method we
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
        mktemplates(network_code=net_code, public_ids=idlist)
