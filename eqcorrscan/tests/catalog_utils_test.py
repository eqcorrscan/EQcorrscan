"""
Functions to test the functions within the eqcorrscan.utils.catalog_utils \
submodule.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class CatalogUtilsTests(unittest.TestCase):
    def test_filter_picks(self):
        """ Test various methods of filetring picks in a catalog."""
        from obspy.clients.fdsn import Client
        from eqcorrscan.utils.catalog_utils import filter_picks
        from obspy import UTCDateTime
        client = Client(str("NCEDC"))
        t1 = UTCDateTime(2004, 9, 28)
        t2 = t1 + 86400
        catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3,
                                    minlatitude=35.7, maxlatitude=36.1,
                                    minlongitude=-120.6, maxlongitude=-120.2,
                                    includearrivals=True)
        stations = ['BMS', 'BAP', 'PAG', 'PAN', 'PBI', 'PKY', 'YEG', 'WOF']
        channels = ['SHZ', 'SHN', 'SHE', 'SH1', 'SH2']
        networks = ['NC']
        locations = ['']
        top_n_picks = 5
        filtered_catalog = filter_picks(catalog=catalog, stations=stations,
                                        channels=channels, networks=networks,
                                        locations=locations,
                                        top_n_picks=top_n_picks)
        for event in filtered_catalog:
            for pick in event.picks:
                self.assertTrue(pick.waveform_id.station_code in stations)
                self.assertTrue(pick.waveform_id.channel_code in channels)
                self.assertTrue(pick.waveform_id.network_code in networks)
                self.assertTrue(pick.waveform_id.location_code in locations)
        filtered_catalog = filter_picks(catalog=catalog,
                                        top_n_picks=top_n_picks)
        filtered_stations = []
        for event in filtered_catalog:
            for pick in event.picks:
                filtered_stations.append(pick.waveform_id.station_code)
        self.assertEqual(len(list(set(filtered_stations))), top_n_picks)

if __name__ == '__main__':
    unittest.main()
