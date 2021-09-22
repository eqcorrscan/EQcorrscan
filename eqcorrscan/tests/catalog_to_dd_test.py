"""
Functions to test the functions within the eqcorrscan.utils.catalog_to_dd.py \
submodule.  Uses test data distributed with the EQcorrscan package.
"""
import unittest
import os
import numpy as np

from collections import Counter, namedtuple
from http.client import IncompleteRead

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

from eqcorrscan.utils.catalog_to_dd import (
    write_catalog, write_correlations, read_phase, write_event, _DTObs,
    _EventPair, _generate_event_id_mapper, _make_sparse_event,
    _prepare_stream, compute_differential_times, _filter_stream,
    _hypodd_event_str, write_phase, write_station)


class TestHelperObjects(unittest.TestCase):
    def test_dtobs(self):
        dtobs = _DTObs(station="FOZ", tt1=3.268, tt2=1.2857650,
                       weight=0.873265, phase="P")
        self.assertEqual(dtobs.ct_string, "FOZ       3.268   1.286 0.8733 P")
        self.assertEqual(dtobs.cc_string, "FOZ       1.982 0.8733 P")

    def test_event_pair(self):
        event_pair = _EventPair(event_id_1=12, event_id_2=54)
        event_pair.obs = [
            _DTObs(station="FOZ", tt1=3.268, tt2=1.2857650,
                   weight=0.873265, phase="P"),
            _DTObs(station="GCSZ", tt1=0.263, tt2=1.50,
                   weight=1.0, phase="S")]
        self.assertEqual(
            event_pair.ct_string,
            '#        12        54\nFOZ       3.268   1.286 0.8733 P\n'
            'GCSZ      0.263   1.500 1.0000 S')
        self.assertEqual(
            event_pair.cc_string,
            '#        12        54 0.0\nFOZ       1.982 0.8733 P\n'
            'GCSZ     -1.237 1.0000 S')

    def test_event_pair_spurious_phases(self):
        """ Make sure that only P and S phases are included. """
        event_pair = _EventPair(event_id_1=12, event_id_2=54)
        with self.assertRaises(AssertionError):
            event_pair.obs = [
                _DTObs(station="FOZ", tt1=3.268, tt2=1.2857650,
                       weight=0.873265, phase="P"),
                _DTObs(station="GCSZ", tt1=0.263, tt2=1.50,
                       weight=1.0, phase="S"),
                _DTObs(station="GCSZ", tt1=0.263, tt2=1.50,
                       weight=1.0, phase="IAML"),
            ]


class TestCatalogMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        starttime = UTCDateTime(2019, 8, 12, 10)
        endtime = UTCDateTime(2019, 8, 13)
        client = Client("GEONET")
        catalog = client.get_events(
            starttime=starttime, endtime=endtime, latitude=-44.5,
            longitude=167.9, maxradius=0.2)
        catalog.events.sort(
            key=lambda e: (e.preferred_origin() or e.origins[0]).time)
        StationInfo = namedtuple(
            "StationInfo", ["network", "station", "location"])
        picked_stations = [StationInfo(p.waveform_id.network_code,
                                       p.waveform_id.station_code,
                                       p.waveform_id.location_code)
                           for ev in catalog for p in ev.picks]
        stations_to_download = [
            sta for sta, _ in Counter(picked_stations).most_common(5)]
        streams = []
        for event in catalog[0:10]:  # Just get the first 10 events
            print(event.preferred_origin().time)
            stream = Stream()
            for station in stations_to_download:
                print(station)
                try:
                    stream += client.get_waveforms(
                        network=station.network, station=station.station,
                        location=station.location, channel="HH?",
                        starttime=event.preferred_origin().time - 10,
                        endtime=event.preferred_origin().time + 80)
                except IncompleteRead:
                    print(f"Could not download {station.station}")
            streams.append(stream)

        picked_stations = set(picked_stations)
        inv_bulk = [(sta.network, sta.station, sta.location, "HH?",
                     starttime, endtime) for sta in picked_stations]
        inventory = client.get_stations_bulk(inv_bulk, level="channel")

        cls.streams = streams
        cls.catalog = catalog
        cls.inventory = inventory

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

    def test_sparse_event(self):
        for event in self.catalog:
            sparse_event = _make_sparse_event(event)
            self.assertEqual(sparse_event.resource_id, event.resource_id.id)
            self.assertEqual(
                sparse_event.origin_time, event.preferred_origin().time)
            self.assertEqual(len(sparse_event.picks), len(event.picks))
            for pick in event.picks:
                matched_pick = [
                    p for p in sparse_event.picks
                    if p.seed_id == pick.waveform_id.get_seed_string()
                    and p.phase == pick.phase_hint]
                self.assertEqual(len(matched_pick), 1)
                self.assertEqual(
                    matched_pick[0].tt, pick.time - sparse_event.origin_time)

    def test_process_stream(self):
        stream = self.streams[0]
        seed_ids = {tr.id for tr in stream}
        event = self.catalog[0]
        extract_len = 10.
        sliced_stream = _prepare_stream(
            stream, event, extract_len=extract_len, pre_pick=1.2,
            seed_pick_ids=None)
        p_picks = [p for p in event.picks if p.phase_hint[0] == "P"
                   and p.waveform_id.get_seed_string() in seed_ids]
        s_picks = [p for p in event.picks if p.phase_hint[0] == "S"
                   and p.waveform_id.get_seed_string() in seed_ids]
        self.assertEqual(len(p_picks), len(sliced_stream["P"]))
        self.assertEqual(len(s_picks), len(sliced_stream["S"]))
        for stream in sliced_stream.values():
            for tr in stream:
                self.assertEqual(
                    tr.stats.endtime - tr.stats.starttime, extract_len)

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

    def test_compute_differential_times(self):
        max_sep = 8.
        diff_times, mapper = compute_differential_times(
            catalog=self.catalog, correlation=False, event_id_mapper=None,
            max_sep=max_sep, min_link=8)
        reverse_mapper = {value: key for key, value in mapper.items()}
        self.assertEqual(len(self.catalog), len(diff_times))
        for key, links in diff_times.items():
            master_id = mapper[key]
            master_event = [e for e in self.catalog
                            if e.resource_id.id == key][0]
            for link in links:
                self.assertEqual(master_id, link.event_id_1)
                linked_event = [
                    e for e in self.catalog
                    if e.resource_id.id == reverse_mapper[link.event_id_2]][0]
                dist, _, _ = gps2dist_azimuth(
                    lat1=master_event.preferred_origin().latitude,
                    lon1=master_event.preferred_origin().longitude,
                    lat2=linked_event.preferred_origin().latitude,
                    lon2=linked_event.preferred_origin().longitude)
                self.assertLess(dist / 1000, max_sep)

    def test_compute_correlation_times(self):
        shift_len = 2
        short_cat = self.catalog[0:10]
        stream_dict = {event.resource_id.id: stream
                       for event, stream in zip(short_cat, self.streams)}
        for interpolate in [True, False]:
            diff_times, mapper = compute_differential_times(
                catalog=short_cat, correlation=True, event_id_mapper=None,
                max_sep=8., min_link=0, min_cc=0.0, stream_dict=stream_dict,
                extract_len=2.0, pre_pick=0.5, shift_len=shift_len,
                interpolate=interpolate, include_master=True)
            diff_times_cat, _ = compute_differential_times(
                catalog=short_cat, correlation=False, event_id_mapper=mapper,
                include_master=True)
            self.assertEqual(len(diff_times), len(short_cat))
            for master_id, linked in diff_times.items():
                for link in linked:
                    cat_link = [pair for pair in diff_times_cat[master_id]
                                if pair.event_id_2 == link.event_id_2][0]
                    if link.event_id_2 == link.event_id_1:
                        # This is the event matched with itself, check that tt1
                        # and tt2 are the same.
                        for obs in link.obs:
                            self.assertTrue(
                                np.allclose(obs.tt1, obs.tt2, atol=0.000001))
                    for obs in link.obs:
                        cat_obs = [o for o in cat_link.obs
                                   if o.station == obs.station and
                                   o.phase == obs.phase][0]
                        self.assertEqual(obs.tt1, cat_obs.tt1)
                        self.assertLessEqual(
                            abs(obs.tt2 - cat_obs.tt2), shift_len)
                        self.assertLessEqual(obs.weight, 1.0)

    def test_compute_correlations_strange_lengths(self):
        """ Check that streams with too short data are unused. PR #424 """
        shift_len = 2
        short_cat = self.catalog[0:10]
        cat_dict = {ev.resource_id.id: ev for ev in short_cat}
        # Need to copy the stream because we are going to trim it!
        stream_dict = {event.resource_id.id: stream.copy()
                       for event, stream in zip(short_cat, self.streams)}
        # Trim to lengths shorter than required.
        shuffler = [-.5, -3, -5, -4, 1, 3, 4, 8, 5, 10]
        keys = list(stream_dict.keys())  # Need to retain order for testing
        for key, shuffle in zip(keys, shuffler):
            st = stream_dict[key]
            ev = cat_dict[key]
            for pick in ev.picks:
                if pick.phase_hint.startswith("P"):
                    _st = st.select(station=pick.waveform_id.station_code)
                    if len(_st) == 0:
                        continue
                    _st.trim(pick.time + shuffle, _st[0].stats.endtime)
        diff_times, mapper = compute_differential_times(
                catalog=short_cat, correlation=True, event_id_mapper=None,
                max_sep=8., min_link=0, min_cc=0.0, stream_dict=stream_dict,
                extract_len=2.0, pre_pick=0.5, shift_len=shift_len,
                interpolate=False, include_master=True)

        self.assertEqual(len(diff_times[keys[-1]]), 0)
        self.assertEqual(len(diff_times[keys[-2]]), 0)
        # There should still be all the events.
        self.assertEqual(len(diff_times), len(short_cat))
        for master_id, linked in diff_times.items():
            for link in linked:
                if link.event_id_2 == link.event_id_1:
                    # This is the event matched with itself, check that tt1
                    # and tt2 are the same.
                    for obs in link.obs:
                        self.assertTrue(
                            np.allclose(obs.tt1, obs.tt2, atol=0.000001))

    def test_write_catalog(self):
        # Contents checked elsewhere
        write_catalog(catalog=self.catalog, event_id_mapper=None,
                      max_sep=8., min_link=8)
        self.assertTrue(os.path.isfile("dt.ct"))
        os.remove("dt.ct")

    def test_write_correlations(self):
        # Contents checked elsewhere
        shift_len = 2
        short_cat = self.catalog[0:10]
        stream_dict = {event.resource_id.id: stream
                       for event, stream in zip(short_cat, self.streams)}
        write_correlations(
            catalog=short_cat, event_id_mapper=None,
            max_sep=8., min_link=0, min_cc=0.0, stream_dict=stream_dict,
            extract_len=2.0, pre_pick=0.5, shift_len=shift_len,
            interpolate=False)
        self.assertTrue(os.path.isfile("dt.cc"))
        os.remove('dt.cc')

    def test_write_correlations_parallel(self):
        # Contents checked elsewhere
        shift_len = 2
        short_cat = self.catalog[0:10]
        stream_dict = {event.resource_id.id: stream
                       for event, stream in zip(short_cat, self.streams)}
        write_correlations(
            catalog=short_cat, event_id_mapper=None,
            max_sep=8., min_link=0, min_cc=0.0, stream_dict=stream_dict,
            extract_len=2.0, pre_pick=0.5, shift_len=shift_len,
            interpolate=False, parallel_process=True, max_workers=2)
        self.assertTrue(os.path.isfile("dt.cc"))
        os.remove('dt.cc')

    def test_filter_stream(self):
        """ Check that original data are unchanged. """
        st = self.streams[0].copy()
        event_id = self.catalog[0].resource_id.id
        freqmin = (2.0, None)
        freqmax = (10.0, None)
        for _freqmin in freqmin:
            for _freqmax in freqmax:
                st_out = _filter_stream(st=st, event_id=event_id,
                                        highcut=_freqmax, lowcut=_freqmin)
                if _freqmax or _freqmin:
                    self.assertNotEqual(st_out[event_id], st)
                else:
                    self.assertEqual(st_out[event_id], st)

    def test_write_phase(self):
        """ This file has been tested with ph2dt """
        write_phase(self.catalog)
        test_data_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        with open(os.path.join(test_data_path, "phase.dat"), "r") as f:
            original_phase = f.read()
        with open("phase.dat") as f:
            phase = f.read()
        self.assertEqual(phase, original_phase)

    def test_write_station(self):
        """ This file has been tested with ph2dt """
        write_station(self.inventory)
        test_data_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        with open(os.path.join(test_data_path, "station.dat"), "r") as f:
            original_station = f.read()
        with open("station.dat") as f:
            station = f.read()
        self.assertEqual(station, original_station)

    def test_write_station_elevations(self):
        """ Include elevations in station.dat: PR #424. """
        write_station(self.inventory, use_elevation=True,
                      filename="station_elev.dat")
        test_data_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        # To do: create this file and test with ph2dt!
        with open(os.path.join(test_data_path, "station_elev.dat"), "r") as f:
            original_station = f.read()
        with open("station_elev.dat") as f:
            station = f.read()
        self.assertEqual(station, original_station)

    def test_write_event(self):
        # Contents checked below
        write_event(self.catalog)
        self.assertTrue(os.path.isfile("event.dat"))
        os.remove("event.dat")

    def test_event_string(self):
        event_str = _hypodd_event_str(event=self.catalog[0], event_id=1)
        self.assertEqual(event_str, '20190812  10352775  -44.5152   167.8586'
                                    '    12.0000   5.51   0.01   0.00   '
                                    '0.74          1')


if __name__ == '__main__':
    unittest.main()
