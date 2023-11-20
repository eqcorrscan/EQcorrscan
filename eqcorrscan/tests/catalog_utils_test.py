"""
Functions to test the functions within the eqcorrscan.utils.catalog_utils \
submodule.
"""
import unittest
import pytest
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

from eqcorrscan.utils.catalog_utils import filter_picks


@pytest.mark.network
class CatalogUtilsTests(unittest.TestCase):
    @pytest.mark.flaky(reruns=2)  # Rerun the test in case of network timeout
    def setUpClass(cls):
        client = Client(str("NCEDC"))
        t1 = UTCDateTime(2004, 9, 28)
        t2 = t1 + 86400
        cls.catalog = client.get_events(
            starttime=t1, endtime=t2, minmagnitude=3,
            minlatitude=35.7, maxlatitude=36.1,
            minlongitude=-120.6, maxlongitude=-120.2,
            includearrivals=True)

    def test_filter_picks(self):
        """ Test various methods of filtering picks in a catalog."""

        stations = ['BMS', 'BAP', 'PAG', 'PAN', 'PBI', 'PKY', 'YEG', 'WOF']
        channels = ['SHZ', 'SHN', 'SHE', 'SH1', 'SH2']
        networks = ['NC']
        locations = ['']
        top_n_picks = 5
        filtered_catalog = filter_picks(
            catalog=self.catalog.copy(), stations=stations,
            channels=channels, networks=networks,
            locations=locations, top_n_picks=top_n_picks)
        for event in filtered_catalog:
            for pick in event.picks:
                self.assertTrue(pick.waveform_id.station_code in stations)
                self.assertTrue(pick.waveform_id.channel_code in channels)
                self.assertTrue(pick.waveform_id.network_code in networks)
                self.assertTrue(pick.waveform_id.location_code in locations)
        filtered_catalog = filter_picks(catalog=self.catalog.copy(),
                                        top_n_picks=top_n_picks)
        filtered_stations = []
        for event in filtered_catalog:
            for pick in event.picks:
                filtered_stations.append(pick.waveform_id.station_code)
        self.assertEqual(len(list(set(filtered_stations))), top_n_picks)

    def test_filter_phase_hints(self):
        filtered_catalog = filter_picks(
            self.catalog.copy(), phase_hints=["P"])

        phase_hints = set(p.phase_hint for ev in filtered_catalog
                          for p in ev.picks)
        self.assertEqual(phase_hints == {"P"})

    def test_filter_single_pick(self):
        filtered_catalog = filter_picks(
            self.catalog.copy(), enforce_single_pick="earliest")

        for ev in filtered_catalog:
            stations = {p.waveform_id.station_code for p in ev.picks}
            for station in stations:
                picks = [p for p in ev.picks if p.waveform_id.station_code == station]
                phase_hints = {p.phase_hint for p in ev.picks}
                for phase_hint in phase_hints:
                    self.assertEqual(
                        1, len([p for p in picks
                               if p.phase_hint == phase_hint]))
        return


if __name__ == '__main__':
    unittest.main()
