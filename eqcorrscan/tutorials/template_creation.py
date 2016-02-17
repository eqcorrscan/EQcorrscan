"""
Simple tutorial detailing how to generate a series of templates from catalog\
data available online.
"""

from eqcorrscan.utils import Sfile_util
from eqcorrscan.core import template_gen
from obspy.fdsn import Client
from obspy import readEvents, UTCDateTime
from obspy.core.event import Catalog

# We want to download some QuakeML files from the New Zealand GeoNet network,
# Obspy doesn't format the request in the way that GeoNet have written it, so
# we have to work around a bit.

client = Client("GEONET")
# We want to download a few events from an earthquake sequence, these are
# identified by publiID numbers
publicIDs = ['2016p008122', '2016p008353', '2016p008155', '2016p008194']
catalog = Catalog()
for publicID in publicIDs:
    data_stream = client._download('http://quakeml.geonet.org.nz/quakeml/' +
                                   '1.2/' + publicID)
    data_stream.seek(0, 0)
    catalog += readEvents(data_stream, format="quakeml")
    data_stream.close()

# Now we should have a catalog of a four events with pick and event information
# We need the seismic data to go with this.  We are going to focus on one day
# of data for this tutorial, all events are from the 4th of January.  We will
# download that day of data for five sites nearby the earthquakes.

t1 = UTCDateTime("2016-01-04T00:00:00.000")
t2 = t1 + 86400

bulk_info = [('NZ', 'BFZ', '*', 'HH*', t1, t2),
             ('NZ', 'CPWZ', '*', 'EH*', t1, t2),
             ('NZ', 'ANWZ', '*', 'EH*', t1, t2),
             ('NZ', 'DVHZ', '*', 'EH*', t1, t2),
             ('NZ', 'PRWZ', '*', 'EH*', t1, t2)]
# Note this will take a little while.
st = client.get_waveforms_bulk(bulk_info)
# Work out what data we actually have to cope with possible lost data
stations = list(set([tr.stats.station for tr in st]))

# Now we can generate the templates.  Note there are wrappers for generating
# templates from QuakeML and seisan format pick files.
for event in catalog:
    picks = []
    for pick in event.picks:
        if pick.waveform_id.station_code in stations:
            picks.append(pick)
    print picks
    # Now we can generate the template, for this tutorial we will plot the
    # templates.
    # Note, this makes more plots than I want, I think there is an issue with
    # this loop
    template = template_gen._template_gen(picks=picks, st=st, length=3.0,
                                          swin='all', prepick=0.05, plot=True)
