"""
Functions to test the mag_calc functions within EQcorrscan.
"""
import unittest
import pytest
import numpy as np
import os
import glob
import logging

from typing import Tuple
from http.client import IncompleteRead

from obspy.core.event import Event, Pick, WaveformStreamID
from obspy import UTCDateTime, read, Trace, read_inventory, Stream
from obspy.clients.fdsn import Client

from eqcorrscan.utils import mag_calc
from eqcorrscan.utils.mag_calc import (
    dist_calc, _sim_WA, _max_p2t, _pairwise, svd_moments, amp_pick_event,
    _snr, relative_amplitude, relative_magnitude, PAZ_WA)
from eqcorrscan.utils.clustering import svd
from eqcorrscan.helpers.mock_logger import MockLoggingHandler


class TestMagCalcMethods(unittest.TestCase):
    """Test all mag_calc functions."""
    def test_dist_calc(self):
        """
        Test the distance calculation that comes with mag_calc.
        """
        self.assertEqual(dist_calc((0, 0, 0), (0, 0, 10)), 10)
        self.assertEqual(round(dist_calc((0, 0, 0), (0, 1, 0))), 111)
        self.assertEqual(round(dist_calc((0, 0, 0), (1, 0, 0))), 111)
        self.assertEqual(round(dist_calc((45, 45, 0), (45, 45, 10))), 10)
        self.assertEqual(round(dist_calc((45, 45, 0), (45, 46, 0))), 79)
        self.assertEqual(round(dist_calc((45, 45, 0), (46, 45, 0))), 111)
        self.assertEqual(round(dist_calc((90, 90, 0), (90, 90, 10))), 10)
        self.assertEqual(round(dist_calc((90, 90, 0), (90, 89, 0))), 0)
        self.assertEqual(round(dist_calc((90, 90, 0), (89, 90, 0))), 111)
        self.assertEqual(round(dist_calc((0, 180, 0), (0, -179, 0))), 111)
        self.assertEqual(round(dist_calc((0, -180, 0), (0, 179, 0))), 111)

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_sim_WA(self):
        """ Test simulating a Wood Anderson instrument. """
        t1 = UTCDateTime("2010-09-3T16:30:00.000")
        t2 = UTCDateTime("2010-09-3T17:00:00.000")
        fdsn_client = Client('IRIS')
        st = fdsn_client.get_waveforms(
            network='NZ', station='BFZ', location='10', channel='HHZ',
            starttime=t1, endtime=t2, attach_response=True)
        inventory = fdsn_client.get_stations(
            network='NZ', station='BFZ', location='10', channel='HHZ',
            starttime=t1, endtime=t2, level="response")
        tr = st[0]
        # Test with inventory
        _sim_WA(trace=tr, inventory=inventory, water_level=10)
        # TODO: This doesn't really test the accuracy

    def test_max_p2t(self):
        """Test the minding of maximum peak-to-trough."""
        data = np.random.randn(1000)
        data[500] = 20
        delta = 0.01
        amplitude, period, time = _max_p2t(data=data, delta=delta)
        self.assertTrue(amplitude > 15)
        self.assertTrue(amplitude < 25)
        self.assertEqual(round(time), 5)

    def test_pairwise(self):
        """Test the itertools wrapper"""
        pairs = _pairwise(range(20))
        for i, pair in enumerate(pairs):
            self.assertEqual(pair[0] + 1, pair[1])
        self.assertEqual(i, 18)

    def test_SVD_mag(self):
        """Test the SVD magnitude calculator."""
        # Do the set-up
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events_processed')
        stream_files = glob.glob(os.path.join(testing_path, '*DFDPC*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        event_list = []
        for i, stream in enumerate(stream_list):
            st_list = []
            for tr in stream:
                if (tr.stats.station, tr.stats.channel) not in\
                        [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]:
                    stream.remove(tr)
                    continue
                st_list.append(i)
            event_list.append(st_list)
        event_list = np.asarray(event_list).T.tolist()
        SVectors, SValues, Uvectors, stachans = svd(stream_list=stream_list)
        M, events_out = svd_moments(u=Uvectors, s=SValues, v=SVectors,
                                    stachans=stachans, event_list=event_list)
        self.assertEqual(len(M), len(stream_list))
        self.assertEqual(len(events_out), len(stream_list))


class TestRelativeAmplitudes(unittest.TestCase):
    """Test the relative amplitude methods."""
    @classmethod
    @pytest.mark.network
    def setUpClass(cls):
        client = Client("GEONET")
        event1 = client.get_events(eventid="2016p912302")[0]
        event2 = client.get_events(eventid="3470170")[0]
        shared_chans = {p1.waveform_id.get_seed_string()
                        for p1 in event1.picks}.intersection(
            {p2.waveform_id.get_seed_string() for p2 in event2.picks})
        bulk = [(p.waveform_id.network_code, p.waveform_id.station_code,
                 p.waveform_id.location_code, p.waveform_id.channel_code,
                 p.time - 20, p.time + 60) for p in event1.picks
                if p.waveform_id.get_seed_string() in shared_chans]
        st1 = client.get_waveforms_bulk(bulk)
        st1 = st1.detrend().filter("bandpass", freqmin=2, freqmax=20)
        bulk = [(p.waveform_id.network_code, p.waveform_id.station_code,
                 p.waveform_id.location_code, p.waveform_id.channel_code,
                 p.time - 20, p.time + 60) for p in event2.picks
                if p.waveform_id.get_seed_string() in shared_chans]
        st2 = client.get_waveforms_bulk(bulk)
        st2 = st2.detrend().filter("bandpass", freqmin=2, freqmax=20)
        cls.event1 = event1
        cls.event2 = event2
        cls.st1 = st1
        cls.st2 = st2

    def test_snr(self):
        noise = np.random.randn(100)
        signal = np.random.randn(100)
        signal[50] = 100
        trace = Trace(data=np.concatenate([noise, signal]))
        trace.stats.sampling_rate = 1.
        snr = _snr(
            tr=trace,
            noise_window=(trace.stats.starttime, trace.stats.starttime + 100),
            signal_window=(trace.stats.starttime + 100, trace.stats.endtime))
        self.assertLessEqual(abs(100 - snr), 30)

    def test_scaled_event(self):
        scale_factor = 0.2
        st1 = read()
        st1.filter("bandpass", freqmin=2, freqmax=20)
        st2 = st1.copy()
        event1 = Event(picks=[
            Pick(time=tr.stats.starttime + 5, phase_hint="P",
                 waveform_id=WaveformStreamID(seed_string=tr.id))
            for tr in st1])
        event2 = event1
        for tr in st2:
            tr.data *= scale_factor
        relative_amplitudes = relative_amplitude(
            st1=st1, st2=st2, event1=event1, event2=event2)
        self.assertEqual(len(relative_amplitudes), len(st1))
        for value in relative_amplitudes.values():
            self.assertAlmostEqual(value, scale_factor)

    def test_no_suitable_picks_event1(self):
        scale_factor = 0.2
        st1 = read()
        st1.filter("bandpass", freqmin=2, freqmax=20)
        st2 = st1.copy()
        event1 = Event(picks=[
            Pick(time=tr.stats.starttime + 5, phase_hint="S",
                 waveform_id=WaveformStreamID(seed_string=tr.id))
            for tr in st1])
        event2 = event1
        for tr in st2:
            tr.data *= scale_factor
        relative_amplitudes = relative_amplitude(
            st1=st1, st2=st2, event1=event1, event2=event2)
        self.assertEqual(len(relative_amplitudes), 0)

    def test_no_suitable_picks_event2(self):
        import copy

        scale_factor = 0.2
        st1 = read()
        st1.filter("bandpass", freqmin=2, freqmax=20)
        st2 = st1.copy()
        event1 = Event(picks=[
            Pick(time=tr.stats.starttime + 5, phase_hint="P",
                 waveform_id=WaveformStreamID(seed_string=tr.id))
            for tr in st1])
        event2 = copy.deepcopy(event1)
        for pick in event2.picks:
            pick.phase_hint = "S"
        for tr in st2:
            tr.data *= scale_factor
        relative_amplitudes = relative_amplitude(
            st1=st1, st2=st2, event1=event1, event2=event2)
        self.assertEqual(len(relative_amplitudes), 0)

    def test_no_picks_event2(self):
        import copy

        scale_factor = 0.2
        st1 = read()
        st1.filter("bandpass", freqmin=2, freqmax=20)
        st2 = st1.copy()
        event1 = Event(picks=[
            Pick(time=tr.stats.starttime + 5, phase_hint="P",
                 waveform_id=WaveformStreamID(seed_string=tr.id))
            for tr in st1])
        event2 = copy.deepcopy(event1)
        event2.picks = []
        for tr in st2:
            tr.data *= scale_factor
        relative_amplitudes = relative_amplitude(
            st1=st1, st2=st2, event1=event1, event2=event2)
        self.assertEqual(len(relative_amplitudes), 0)

    def test_no_picks_event1(self):
        scale_factor = 0.2
        st1 = read()
        st1.filter("bandpass", freqmin=2, freqmax=20)
        st2 = st1.copy()
        event1 = Event()
        event2 = event1
        for tr in st2:
            tr.data *= scale_factor
        relative_amplitudes = relative_amplitude(
            st1=st1, st2=st2, event1=event1, event2=event2)
        self.assertEqual(len(relative_amplitudes), 0)

    def test_low_snr(self):
        scale_factor = 0.2
        st1 = read()
        st1[0].data += np.random.randn(st1[0].stats.npts) * st1[0].data.max()
        st2 = st1.copy()
        st2[1].data += np.random.randn(st2[1].stats.npts) * st2[1].data.max()
        event1 = Event(picks=[
            Pick(time=tr.stats.starttime + 5, phase_hint="P",
                 waveform_id=WaveformStreamID(seed_string=tr.id))
            for tr in st1])
        event2 = event1
        for tr in st2:
            tr.data *= scale_factor
        relative_amplitudes = relative_amplitude(
            st1=st1, st2=st2, event1=event1, event2=event2)
        self.assertEqual(len(relative_amplitudes), 1)
        for value in relative_amplitudes.values():
            self.assertAlmostEqual(value, scale_factor)

    def test_real_near_repeat(self):
        relative_amplitudes = relative_amplitude(
            st1=self.st1, st2=self.st2, event1=self.event1, event2=self.event2)
        for seed_id, ratio in relative_amplitudes.items():
            self.assertLess(abs(0.8 - ratio), 0.1)

    def test_real_near_repeat_magnitudes(self):
        relative_magnitudes, correlations = relative_magnitude(
            st1=self.st1, st2=self.st2, event1=self.event1, event2=self.event2,
            return_correlations=True)
        for seed_id, mag_diff in relative_magnitudes.items():
            self.assertLess(abs(mag_diff + 0.15), 0.1)

    def test_real_near_repeat_magnitudes_corr_provided(self):
        relative_magnitudes = relative_magnitude(
            st1=self.st1, st2=self.st2, event1=self.event1, event2=self.event2,
            correlations={tr.id: 1.0 for tr in self.st1})
        for seed_id, mag_diff in relative_magnitudes.items():
            self.assertLess(abs(mag_diff + 0.15), 0.1)

    def test_real_near_repeat_magnitudes_bad_corr(self):
        relative_magnitudes = relative_magnitude(
            st1=self.st1, st2=self.st2, event1=self.event1, event2=self.event2,
            correlations={tr.id: 0.2 for tr in self.st1})
        self.assertEqual(len(relative_magnitudes), 0)

    def test_real_near_repeat_magnitudes_S_picks(self):
        self.event1.picks.append(self.event1.picks[0].copy())
        self.event1.picks[-1].phase_hint = "S"
        self.event1.picks[-1].waveform_id.channel_code = "EHN"
        self.st1.append(self.st1[0].copy())
        self.st1[-1].id = self.st1[-1].id[0:-1] + "N"
        self.event2.picks.append(self.event2.picks[13].copy())
        self.event2.picks[-1].phase_hint = "S"
        self.event2.picks[-1].waveform_id.channel_code = "EHN"
        self.st2.append(
            self.st2.select(station=self.st1[0].stats.station)[0].copy())
        self.st2[-1].id = self.st2[-1].id[0:-1] + "N"
        # Check that the results are not equal when using or not using S:
        relative_magnitudes1 = relative_magnitude(
            st1=self.st1, st2=self.st2, event1=self.event1, event2=self.event2,
            use_s_picks=False)
        relative_magnitudes2 = relative_magnitude(
            st1=self.st1, st2=self.st2, event1=self.event1, event2=self.event2,
            use_s_picks=True)
        self.assertNotEqual(relative_magnitudes1, relative_magnitudes2)


class TestAmpPickEvent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        log = logging.getLogger(mag_calc.__name__)
        cls._log_handler = MockLoggingHandler(level='DEBUG')
        log.addHandler(cls._log_handler)
        cls.log_messages = cls._log_handler.messages
        client = Client("GEONET")
        cls.event = client.get_events(eventid="2019p498440")[0]
        origin_time = cls.event.preferred_origin().time
        bulk = [(
            p.waveform_id.network_code, p.waveform_id.station_code,
            p.waveform_id.location_code, p.waveform_id.channel_code,
            origin_time - 10, origin_time + 120)
            for p in cls.event.picks]
        cls.inventory = client.get_stations_bulk(bulk, level='response')
        cls.st = Stream()
        for _bulk in bulk:
            try:
                cls.st += client.get_waveforms(*_bulk)
            except IncompleteRead:
                print(f"Could not download {_bulk}")
        cls.available_stations = len({p.waveform_id.station_code
                                      for p in cls.event.picks})

    def setUp(self):
        self._log_handler.reset()

    def test_amp_pick_event(self):
        """Test the main amplitude picker."""
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy(),
            inventory=self.inventory)
        self.assertEqual(len(picked_event.picks),
                         len(self.st) + self.available_stations)

    def test_amp_pick_remove_old_picks(self):
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy(),
            inventory=self.inventory, remove_old=True)
        self.assertEqual(len(picked_event.amplitudes), self.available_stations)

    def test_amp_pick_missing_channel(self):
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy()[0:-3],
            inventory=self.inventory, remove_old=True)
        missed = False
        for warning in self.log_messages['warning']:
            if 'not found in the stream' in warning:
                missed = True
        self.assertTrue(missed)
        self.assertEqual(len(picked_event.amplitudes),
                         self.available_stations - 1)

    def test_amp_pick_not_varwin(self):
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy(),
            inventory=self.inventory, remove_old=True,
            var_wintype=False)
        self.assertEqual(len(picked_event.amplitudes), self.available_stations)

    def test_amp_pick_not_varwin_no_S(self):
        event = self.event.copy()
        for pick in event.picks:
            if pick.phase_hint.upper() == 'S':
                event.picks.remove(pick)
        picked_event = amp_pick_event(
            event=event, st=self.st.copy(), inventory=self.inventory,
            remove_old=True, var_wintype=False)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_varwin_no_S(self):
        event = self.event.copy()
        for pick in event.picks:
            if pick.phase_hint.upper() == 'S':
                event.picks.remove(pick)
        picked_event = amp_pick_event(
            event=event, st=self.st.copy(), inventory=self.inventory,
            remove_old=True, var_wintype=True)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_varwin_no_P(self):
        event = self.event.copy()
        for pick in event.picks:
            if pick.phase_hint.upper() == 'P':
                event.picks.remove(pick)
        picked_event = amp_pick_event(
            event=event, st=self.st.copy(), inventory=self.inventory,
            remove_old=True, var_wintype=True)
        self.assertEqual(len(picked_event.amplitudes), self.available_stations)

    def test_amp_pick_high_min_snr(self):
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy(),
            inventory=self.inventory, remove_old=True, var_wintype=False,
            min_snr=15)
        self.assertEqual(len(picked_event.amplitudes), 0)

    def test_amp_pick_no_prefilt(self):
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy(),
            inventory=self.inventory, remove_old=True,
            var_wintype=False, pre_filt=False)
        self.assertEqual(len(picked_event.amplitudes), 5)


class TestAmpPickAccuracy(unittest.TestCase):
    sampling_rate = 100.0  # Sampling-rate in Hz for synthetic
    length = 60.0  # Length in seconds for synthetic
    inv = read_inventory().select(
        station="RJOB", channel="EHZ", starttime=UTCDateTime(2020, 1, 1))

    def simulate_trace(
        self,
        max_amplitude: float = 1e-4,
        width: float = 0.2,
        noise: bool = False,
    ) -> Tuple[Trace, Event, float, float]:
        """ Make a wood-anderson trace and convolve it with a response. """
        from scipy.signal import ricker

        # Make dummy data in meters on Wood Anderson
        np.random.seed(42)
        if noise:
            data = np.random.randn(int(self.sampling_rate * self.length))
        else:
            data = np.zeros(int(self.sampling_rate * self.length))
        wavelet = ricker(
            int(self.sampling_rate * self.length),
            a=(width * self.sampling_rate / 4))
        wavelet /= wavelet.max()
        wavelet *= max_amplitude
        vals = sorted([(val, ind) for ind, val in enumerate(wavelet)],
                      key=lambda x: x[0])
        period = abs(vals[0][1] - vals[1][1]) / self.sampling_rate
        half_max = (wavelet.max() - wavelet.min()) / 2

        data += wavelet

        # Make dummy trace
        tr = Trace(data=data, header=dict(
            station=self.inv[0][0].code,
            network=self.inv[0].code,
            location=self.inv[0][0][0].location_code,
            channel=self.inv[0][0][0].code,
            sampling_rate=self.sampling_rate,
            starttime=UTCDateTime(2020, 1, 1)))
        tr = tr.detrend()

        # Remove Wood-Anderson response and simulate seismometer
        resp = self.inv[0][0][0].response
        paz = resp.get_paz()
        paz = {
            'poles': paz.poles,
            'zeros': paz.zeros,
            'gain': paz.normalization_factor,
            'sensitivity': resp.instrument_sensitivity.value,
        }
        tr = tr.simulate(
            paz_remove=PAZ_WA, paz_simulate=paz)

        # Make an event
        mid_point = tr.stats.starttime + 0.5 * (
                tr.stats.endtime - tr.stats.starttime)
        event = Event(picks=[
            Pick(phase_hint="P",
                 time=mid_point - 3.0,
                 waveform_id=WaveformStreamID(
                     network_code=tr.stats.network,
                     station_code=tr.stats.station,
                     location_code=tr.stats.location,
                     channel_code=tr.stats.channel)),
            Pick(phase_hint="S", time=mid_point,
                 waveform_id=WaveformStreamID(
                     network_code=tr.stats.network,
                     station_code=tr.stats.station,
                     location_code=tr.stats.location,
                     channel_code=tr.stats.channel))])
        return tr, event, period, half_max

    def test_no_filter(self):
        max_amplitude, width = 1e-3, 0.2
        tr, event, period_at_max, half_max = self.simulate_trace(
            max_amplitude=max_amplitude, width=width)
        event_out = amp_pick_event(
            event=event, st=Stream(tr), inventory=self.inv,
            pre_filt=False)
        self.assertLessEqual(
            abs(event_out.amplitudes[0].period - period_at_max), 0.05)
        self.assertLessEqual(
            abs(event_out.amplitudes[0].generic_amplitude - half_max), 1e-5)

    def test_range_of_freqs(self):
        max_amplitude, width = 1e-3, 0.2
        tr, event, period_at_max, half_max = self.simulate_trace(
            max_amplitude=max_amplitude, width=width)
        for lowcut in range(1, 5, 1):
            for highcut in range(5, 50, 2):
                if lowcut >= highcut:
                    continue
                event_out = amp_pick_event(
                    event=event.copy(), st=Stream(tr), inventory=self.inv,
                    pre_filt=True, lowcut=lowcut, highcut=highcut)
                self.assertLessEqual(
                    abs(event_out.amplitudes[0].period - period_at_max), 0.1)
                self.assertLessEqual(
                    abs(event_out.amplitudes[0].generic_amplitude - half_max),
                    1e-3)

    def test_iaspei(self):
        max_amplitude, width = 1e-3, 0.2
        tr, event, period_at_max, half_max = self.simulate_trace(
            max_amplitude=max_amplitude, width=width)
        event_out = amp_pick_event(
            event=event, st=Stream(tr), inventory=self.inv,
            pre_filt=False, iaspei_standard=True)
        self.assertLessEqual(
            abs((event_out.amplitudes[0].generic_amplitude *
                 PAZ_WA["sensitivity"]) - half_max), 1e-5)

    def test_velocity(self):
        # TODO
        return


if __name__ == '__main__':
    unittest.main()
