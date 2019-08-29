"""
Functions to test the functions within the eqcorrscan.utils.catalog_to_dd.py \
submodule.  Uses test data distributed with the EQcorrscan package.
"""
import unittest
import os

from collections import Counter

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Arrival

from eqcorrscan.utils.catalog_to_dd import (
    write_catalog, write_correlations, read_phase, write_event, _DTObs,
    _EventPair, _generate_event_id_mapper, _get_arrival_for_pick,
    _combined_weight, _prepare_stream, _compute_dt,
    _compute_dt_correlations, _make_event_pair, compute_differential_times,
    _filter_stream, _hypodd_phase_pick_str, _hypodd_event_str,
    _hypodd_phase_str)


class TestHelperObjects(unittest.TestCase):
    def test_dtobs(self):
        dtobs = _DTObs(station="FOZ", tt1=3.268, tt2=1.2857650,
                       weight=0.873265, phase="P")
        self.assertEqual(dtobs.ct_string, "FOZ     3.268   1.286 0.8733 P")
        self.assertEqual(dtobs.cc_string, "FOZ     1.982 0.8733 P")

    def test_event_pair(self):
        event_pair = _EventPair(event_id_1=12, event_id_2=54)
        event_pair.obs = [
            _DTObs(station="FOZ", tt1=3.268, tt2=1.2857650,
                   weight=0.873265, phase="P"),
            _DTObs(station="GCSZ", tt1=0.263, tt2=1.50,
                   weight=1.0, phase="S")]
        self.assertEqual(
            event_pair.ct_string,
            '#        12        54\nFOZ     3.268   1.286 0.8733 P\n'
            'GCSZ    0.263   1.500 1.0000 S')
        self.assertEqual(
            event_pair.cc_string,
            '#        12        54 0.0\nFOZ     1.982 0.8733 P\n'
            'GCSZ   -1.237 1.0000 S')


# TODO: Heaps of tests!
class TestCatalogMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        client = Client("GEONET")
        catalog = client.get_events(
            starttime=UTCDateTime(2019, 8, 12, 10),
            endtime=UTCDateTime(2019, 8, 13),
            latitude=-44.5,
            longitude=167.9,
            maxradius=0.2)
        stations_to_download = [sta for sta, _ in Counter(
            [p.waveform_id.station_code
             for ev in catalog for p in ev.picks]).most_common(5)]
        streams = []  # TODO: get some streams!

        cls.streams = streams
        cls.catalog = catalog

    def test_id_mapper(self):
        map_one = _generate_event_id_mapper(self.catalog)
        for event in self.catalog:
            self.assertIn(event.resource_id.id, map_one.keys())
        self.assertEqual(len(self.catalog), len(map_one))
        map_one.pop(self.catalog[10].resource_id.id)
        self.assertNotEqual(len(self.catalog), len(map_one))
        map_two = _generate_event_id_mapper(self.catalog, map_one)
        self.assertEqual(len(self.catalog), len(map_two))
        for event in self.catalog:
            self.assertIn(event.resource_id.id, map_two.keys())
        # Check that a new event mapping has been made.
        self.assertEqual(len(self.catalog) + 1,
                         map_two.pop(self.catalog[10].resource_id.id))

    def test_get_arrival_for_pick(self):
        event = self.catalog[0]
        arr = _get_arrival_for_pick(event, event.picks[0])
        self.assertEqual(arr.pick_id.get_referred_object(), event.picks[0])

    def test_weight_combination(self):
        arr1 = Arrival(time_weight=1.0)
        arr2 = Arrival(time_weight=0.5)
        self.assertEqual(_combined_weight(arr1, arr2), 0.75)
        self.assertEqual(_combined_weight(arr1, None), 1.0)
        self.assertEqual(_combined_weight(None, arr2), 0.75)

    def test_read_phase(self):
        """Function to test the phase reading function"""
        test_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                 'test_data', 'tunnel.phase')
        test_catalog = read_phase(test_file)
        self.assertEqual(len(test_catalog), 2)
        self.assertEqual(test_catalog[0].origins[0].latitude, -43.169)
        self.assertEqual(test_catalog[0].origins[0].longitude, 170.646)
        self.assertEqual(test_catalog[0].origins[0].depth, -288)
        self.assertEqual(test_catalog[0].origins[0].time,
                         UTCDateTime('2012-01-30T01:45:43.25'))
        self.assertEqual(test_catalog[1].origins[0].latitude, -43.169)
        self.assertEqual(test_catalog[1].origins[0].longitude, 170.646)
        self.assertEqual(test_catalog[1].origins[0].depth, -288)
        self.assertEqual(test_catalog[1].origins[0].time,
                         UTCDateTime('2012-01-30T06:48:43.07'))


if __name__ == '__main__':
    unittest.main()
