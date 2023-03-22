"""
A series of test functions for the core functions in EQcorrscan.
"""
import glob
import os
import shutil
import unittest
import numpy as np
import logging

from obspy import read, UTCDateTime

from eqcorrscan.utils.synth_seis import generate_synth_data
from eqcorrscan.utils.correlate import get_stream_xcorr
from eqcorrscan.core import lag_calc as lag_calc_module
from eqcorrscan.core.lag_calc import (
    _xcorr_interp, LagCalcError, _prepare_data, lag_calc, xcorr_pick_family,
    _concatenate_and_correlate)
from eqcorrscan.core.match_filter import Detection, Family, Party, Template
from eqcorrscan.helpers.mock_logger import MockLoggingHandler


class SyntheticTests(unittest.TestCase):
    """ Test lag-calc with synthetic data. """
    @classmethod
    def setUpClass(cls):
        np.random.seed(999)
        print("Setting up class")
        samp_rate = 50
        t_length = .75
        # Make some synthetic templates
        templates, data, seeds = generate_synth_data(
            nsta=5, ntemplates=5, nseeds=10, samp_rate=samp_rate,
            t_length=t_length, max_amp=10, max_lag=15, phaseout="both",
            jitter=0, noise=False, same_phase=True)
        # Rename channels
        channel_mapper = {"SYN_Z": "HHZ", "SYN_H": "HHN"}
        for tr in data:
            tr.stats.channel = channel_mapper[tr.stats.channel]
        for template in templates:
            for tr in template:
                tr.stats.channel = channel_mapper[tr.stats.channel]
        party = Party()
        t = 0
        data_start = data[0].stats.starttime
        for template, template_seeds in zip(templates, seeds):
            template_name = "template_{0}".format(t)
            detections = []
            for i, sample in enumerate(template_seeds["time"]):
                det = Detection(
                    template_name=template_name,
                    detect_time=data_start + (sample / samp_rate),
                    detect_val=template_seeds["SNR"][i] / len(data),
                    no_chans=len(data),
                    chans=[
                        (tr.stats.station, tr.stats.channel) for tr in data],
                    threshold=0.0, threshold_input=0.0, threshold_type="abs",
                    typeofdet="ccc")
                det._calculate_event(
                    template_st=template, estimate_origin=False)
                detections.append(det)
            # Make a fully formed Template
            _template = Template(
                name=template_name, st=template, lowcut=2.0, highcut=15.0,
                samp_rate=samp_rate, filt_order=4, process_length=86400,
                prepick=10. / samp_rate, event=None)
            family = Family(template=_template, detections=detections)
            party += family
            t += 1
        cls.party = party
        cls.data = data
        cls.t_length = t_length
        assert len(data) == 10

    def _prepare_data_checks(self, detect_stream_dict, shift_len, family):
        self.assertEqual(
            len(detect_stream_dict), len(family.detections))
        for detection_id, stream in detect_stream_dict.items():
            detection = [d for d in family if d.id == detection_id][0]
            self.assertEqual(len(stream), detection.no_chans)
            for tr in stream:
                self.assertAlmostEqual(
                    tr.stats.endtime - tr.stats.starttime,
                    (2 * shift_len) + self.t_length, 1)
                pick = [p for p in detection.event.picks
                        if p.waveform_id.get_seed_string() == tr.id][0]
                self.assertAlmostEqual(
                    tr.stats.starttime,
                    pick.time - (family.template.prepick + shift_len), 1)

    def test_prepare_data(self):
        shift_len = 0.2
        for family in self.party:
            detect_stream_dict = _prepare_data(
                family=family, detect_data=self.data, shift_len=shift_len)
            self._prepare_data_checks(detect_stream_dict=detect_stream_dict,
                                      family=family, shift_len=shift_len)

    def test_prepare_data_too_short(self):
        data = self.data.copy()
        data.trim(self.party[0][0].detect_time,
                  self.party[0][0].detect_time + 3)
        shift_len = 0.2
        detect_stream_dict = _prepare_data(
            family=self.party[0], detect_data=data, shift_len=shift_len)
        self.assertEqual(len(detect_stream_dict), 1)

    def test_prepare_data_masked(self):
        data = self.data.copy()
        data.cutout(self.party[0][0].detect_time,
                    self.party[0][0].detect_time + 3)
        data.merge()
        shift_len = 0.2
        detect_stream_dict = _prepare_data(
            family=self.party[0], detect_data=data, shift_len=shift_len)
        short_key = self.party[0][0].id
        for key, value in detect_stream_dict.items():
            detection = [d for d in self.party[0] if d.id == key][0]
            if key == short_key:
                self.assertNotEqual(len(value), detection.no_chans)
            else:
                self.assertEqual(len(value), detection.no_chans)

    def test_family_picking(self):
        catalog_dict = xcorr_pick_family(
            family=self.party[0], stream=self.data, shift_len=0.2, plot=False,
            export_cc=False)
        for event in catalog_dict.values():
            self.assertEqual(len(event.picks), len(self.data))
            for pick in event.picks:
                self.assertTrue("cc_max=" in pick.comments[0].text)
                self.assertAlmostEqual(
                    float(pick.comments[0].text.split("=")[-1]), 1.0, 1)

    def test_family_picking_missing_data(self):
        """ Check that this all works okay when one channel has some gaps """
        gappy_data = self.data.copy()
        sr = gappy_data[0].stats.sampling_rate
        mask = np.zeros(len(gappy_data[0].data), dtype=bool)
        gap_start = self.party[0][5].detect_time - 5
        gap_length = 3600
        gap_end = gap_start + gap_length
        gap_start = int(sr * (gap_start - UTCDateTime(0)))
        gap_end = int(sr * (gap_end - UTCDateTime(0)))
        mask[gap_start: gap_end] = np.ones(int(gap_length * sr), dtype=bool)
        gappy_data[0].data = np.ma.masked_array(
            data=gappy_data[0].data, mask=mask)
        catalog_dict = xcorr_pick_family(
            family=self.party[0], stream=gappy_data, shift_len=0.2, plot=False,
            export_cc=False)
        gap = gappy_data.split().get_gaps()
        for event in catalog_dict.values():
            if len(event.picks) != len(gappy_data):
                self.assertEqual(len(event.picks), len(gappy_data) - 1)
                # Check that the event happened to be in the gap
                self.assertTrue(gap[0][4] <= event.picks[0].time <= gap[0][5])
                # Check that there isn't a pick on the channel missing data
                self.assertNotIn('.'.join(gap[0][0:4]),
                                 [pick.waveform_id.get_seed_string()
                                  for pick in event.picks])
            for pick in event.picks:
                self.assertTrue("cc_max=" in pick.comments[0].text)
                self.assertAlmostEqual(
                    float(pick.comments[0].text.split("=")[-1]), 1.0, 1)

    def test_family_picking_with_interpolation(self):
        catalog_dict = xcorr_pick_family(
            family=self.party[0], stream=self.data, shift_len=0.2, plot=False,
            interpolate=True, export_cc=False)
        for event in catalog_dict.values():
            for pick in event.picks:
                self.assertTrue("cc_max=" in pick.comments[0].text)
                self.assertAlmostEqual(
                    float(pick.comments[0].text.split("=")[-1]), 1.0, 1)

    def test_family_picking_with_new_interpolation(self):
        catalog_dict = xcorr_pick_family(
            family=self.party[0], stream=self.data, shift_len=0.2, plot=False,
            interpolate=True, use_new_resamp_method=True, export_cc=False)
        for event in catalog_dict.values():
            for pick in event.picks:
                self.assertTrue("cc_max=" in pick.comments[0].text)
                self.assertAlmostEqual(
                    float(pick.comments[0].text.split("=")[-1]), 1.0, 1)

    def test_lag_calc_api(self):
        detections = [d for f in self.party for d in f]
        templates = [f.template.st for f in self.party]
        template_names = [f.template.name for f in self.party]
        output_cat = lag_calc(
            detections, self.data, template_names, templates,
            shift_len=0.2, min_cc=0.4, min_cc_from_mean_cc_factor=1,
            horizontal_chans=['E', 'N', '1', '2'],
            vertical_chans=['Z'], cores=1, interpolate=False,
            plot=False, export_cc=False)
        self.assertEqual(len(output_cat), len(detections))
        for event in output_cat:
            self.assertEqual(len(event.picks), len(self.data))
            for pick in event.picks:
                self.assertTrue("cc_max=" in pick.comments[0].text)
                self.assertAlmostEqual(
                    float(pick.comments[0].text.split("=")[-1]), 1.0, 1)

    def test_xcorr_pick_family_export_cc(self):
        cc_dir = 'lag_calc_cc_exported'
        xcorr_pick_family(
            family=self.party[0], stream=self.data.copy(), shift_len=0.2,
            plot=False, interpolate=False, export_cc=True, cc_dir=cc_dir)
        cc_files = glob.glob(os.path.join(cc_dir, '*.npy'))
        assert len(cc_files) == len(self.party[0])
        for fcc in cc_files:
            np.load(fcc)
        if os.path.isdir(cc_dir):
            shutil.rmtree(cc_dir)


class SimpleRealDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        st = read().detrend().filter("bandpass", freqmin=2, freqmax=20)
        cls.template = st.slice(starttime=st[0].stats.starttime + 4,
                                endtime=st[0].stats.starttime + 10).copy()
        cls.shift_len = 1.0
        cls.detect_streams = []
        for i in range(1, 8):
            cls.detect_streams.append(
                st.slice(starttime=st[0].stats.starttime + i - cls.shift_len,
                         endtime=st[0].stats.starttime + i + 6 + cls.shift_len
                         ).copy())

    def test_correlation_max_and_position(self):
        ccc, chans = _concatenate_and_correlate(
            streams=self.detect_streams, template=self.template, cores=1)
        samp_rate = self.template[0].stats.sampling_rate
        t_start, t_end = (
            self.template[0].stats.starttime, self.template[0].stats.endtime)
        for _ccc, detect_stream in zip(ccc, self.detect_streams):
            d_start, d_end = (detect_stream[0].stats.starttime,
                              detect_stream[0].stats.endtime)
            if d_start <= t_start and d_end >= t_end:
                for ccc_chan in _ccc:
                    self.assertEqual(round(ccc_chan.max(), 5), 1.0)
                    self.assertEqual(ccc_chan.argmax(),
                                     samp_rate * (t_start - d_start))
            else:
                for ccc_chan in _ccc:
                    self.assertNotEqual(round(ccc_chan.max(), 5), 1.0)

    def test_correlation_precision(self):
        """Compare to correlation function outputs"""
        ccc, chans = _concatenate_and_correlate(
            streams=self.detect_streams, template=self.template, cores=1)
        fftw_xcorr_func = get_stream_xcorr("fftw")
        for _ccc, detect_stream in zip(ccc, self.detect_streams):
            fftw_ccc, _, _ = fftw_xcorr_func(
                templates=[self.template], stream=detect_stream, stack=False)
            for chan_ccc, fftw_chan_ccc in zip(_ccc, fftw_ccc[0]):
                self.assertTrue(np.allclose(
                    chan_ccc, fftw_chan_ccc, atol=.00001))


class ShortTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        log = logging.getLogger(lag_calc_module.__name__)
        cls._log_handler = MockLoggingHandler(level='DEBUG')
        log.addHandler(cls._log_handler)
        cls.log_messages = cls._log_handler.messages

    def setUp(self):
        self._log_handler.reset()

    def test_error(self):
        with self.assertRaises(LagCalcError):
            raise LagCalcError('Generic error')
        err = LagCalcError('Generic error')
        self.assertEqual('Generic error', err.value)
        self.assertEqual('Generic error', err.__repr__())
        self.assertEqual('LagCalcError: Generic error', err.__str__())

    def test_bad_interp(self):
        ccc = np.array([-0.21483282, -0.59443731, 0.1898917, -0.67516038,
                        0.60129057, -0.71043723,  0.16709118, 0.96839009,
                        1.58283915, -0.3053663])

        _xcorr_interp(ccc, 0.1, use_new_resamp_method=True)
        self.assertEqual(len(self.log_messages['warning']), 1)
        self.assertTrue(
            'not give an accurate result' in self.log_messages['warning'][0])


if __name__ == '__main__':
    unittest.main()
