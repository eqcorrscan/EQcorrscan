"""
Test of large-scale match-filter, not to be run on CI.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import os
import unittest
import warnings

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core.event import Pick

from eqcorrscan.tutorials.get_geonet_events import get_geonet_events
from eqcorrscan.utils import catalog_utils
from eqcorrscan.core import template_gen
from eqcorrscan.utils import pre_processing


class TestGeoNetCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        client = Client('GEONET')
        t1 = UTCDateTime(2016, 9, 4)
        t2 = t1 + (86400 * 20)
        catalog = get_geonet_events(
            startdate=t1, enddate=t2, minmag=4, minlat=-49, maxlat=-35,
            minlon=175.0, maxlon=185.0)
        catalog = catalog_utils.filter_picks(
            catalog, channels=['EHZ'], top_n_picks=5)
        for event in catalog:
            extra_pick = Pick()
            extra_pick.phase_hint = 'S'
            extra_pick.time = event.picks[0].time + 10
            extra_pick.waveform_id = event.picks[0].waveform_id
            event.picks.append(extra_pick)
        cls.templates = template_gen.from_client(
            catalog=catalog, client_id='GEONET', lowcut=2.0, highcut=9.0,
            samp_rate=50.0, filt_order=4, length=3.0, prepick=0.15, swin='all',
            process_len=3600)
        # Download and process the day-long data
        bulk_info = [(tr.stats.network, tr.stats.station, '*',
                      tr.stats.channel, t1 + (4 * 3600), t1 + (5 * 3600))
                     for tr in cls.templates[0]]
        # Just downloading an hour of data
        print('Downloading data')
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        st.trim(t1 + (4 * 3600), t1 + (5 * 3600)).sort()
        # This is slow?
        print('Processing continuous data')
        cls.st = pre_processing.shortproc(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=50.0,
            debug=0, num_cores=1)
        cls.st.trim(t1 + (4 * 3600), t1 + (5 * 3600)).sort()
        cls.template_names = [str(template[0].stats.starttime)
                              for template in cls.templates]


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()