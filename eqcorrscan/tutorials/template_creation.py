"""
Simple tutorial detailing how to generate a series of templates from catalog\
data available online.
"""

from collections import Counter
from eqcorrscan.core import template_gen
from eqcorrscan.utils import pre_processing
from joblib import Parallel, delayed

# This import section copes with namespace changes between obspy versions
import obspy
if int(obspy.__version__.split('.')[0]) >= 1:
    from obspy.clients.fdsn import Client
    from obspy import read_events
else:
    from obspy.fdsn import Client
    from obspy import readEvents as read_events
from obspy import UTCDateTime, Stream
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
    catalog += read_events(data_stream, format="quakeml")
    data_stream.close()

# Now we should have a catalog of a four events with pick and event information
# We need the seismic data to go with this.  We are going to focus on one day
# of data for this tutorial, all events are from the 4th of January.  We will
# download that day of data for five sites nearby the earthquakes.

t1 = UTCDateTime("2016-01-04T00:00:00.000")
t2 = t1 + 86400

# Work out what stations we have picks for all events from

all_picks = []
for event in catalog:
    all_picks += [(pick.waveform_id.station_code,
                   pick.waveform_id.channel_code) for pick in event.picks]
all_picks = Counter(all_picks).most_common(5)
all_picks = [pick[0] for pick in all_picks]

bulk_info = []
for station in all_picks:
    bulk_info.append(('NZ', station[0], '*', station[1][0:2]+'*', t1, t2))

# Note this will take a little while.
print('Downloading seismic data, this may take a while')
st = client.get_waveforms_bulk(bulk_info)
# Merge the stream, it will be downloaded in chunks
st.merge(fill_value='interpolate')
# Work out what data we actually have to cope with possible lost data
stations = list(set([tr.stats.station for tr in st]))

# Pre-process the data to set frequency band and sampling rate
print('Processing the seismic data')
st = Parallel(n_jobs=10)(delayed(pre_processing.dayproc)
                         (tr=tr, lowcut=2.0, highcut=9.0, filt_order=3,
                          samp_rate=20.0, debug=0, starttime=t1)
                         for tr in st)
# Convert from list to stream
st = Stream(st)

# Now we can generate the templates.  Note there are wrappers for generating
# templates from QuakeML and seisan format pick files.
templates = []
for event in catalog:
    picks = []
    for pick in event.picks:
        if pick.waveform_id.station_code in stations:
            picks.append(pick)
    # Now we can generate the template, for this tutorial we will plot the
    # templates.
    # Note, this makes more plots than I want, I think there is an issue with
    # this loop
    templates.append(template_gen._template_gen(picks=picks, st=st, length=3.0,
                                                swin='all', prepick=0.05,
                                                plot=True))

# We now have a series of templates! Using obspys Stream.write() method we
# can save these to disk for later use.  We will do that now for use in the
# following tutorials.
for i, template in enumerate(templates):
    template.write('tutorial_template_' + str(i) + '.ms', format='MSEED')
    # Note that this will warn you about data types.  As we don't care
    # at the moment, whatever obspy chooses is fine.
