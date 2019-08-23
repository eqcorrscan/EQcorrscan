"""
A series of test functions for the core functions in EQcorrscan.
"""
import unittest
import numpy as np
import logging

from eqcorrscan.utils.synth_seis import generate_synth_data
from eqcorrscan.core.lag_calc import (
    _xcorr_interp, LagCalcError, _prepare_data, lag_calc, xcorr_pick_family)
from eqcorrscan.core.match_filter import Detection, Family, Party, Template
from eqcorrscan.helpers.mock_logger import MockLoggingHandler

np.random.seed(999)


class SyntheticTests(unittest.TestCase):
    """ Test lag-calc with synthetic data. """
    @classmethod
    def setUpClass(cls) -> None:
        samp_rate = 50
        cls.t_length = .75
        # Make some synthetic templates
        templates, data, seeds = generate_synth_data(
            nsta=5, ntemplates=5, nseeds=10, samp_rate=samp_rate,
            t_length=cls.t_length, max_amp=10, max_lag=15, phaseout="both",
            jitter=10)
        # Rename channels
        channel_mapper = {"SYN_Z": "HHZ", "SYN_H": "HHN"}
        for tr in data:
            tr.stats.channel = channel_mapper[tr.stats.channel]
        for template in templates:
            for tr in template:
                tr.stats.channel = channel_mapper[tr.stats.channel]
        cls.party = Party()
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
                    chans=[(tr.stats.station, tr.stats.channel) for tr in data],
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
            cls.party += family
            t += 1
        cls.data = data

    def test_prepare_data(self):
        shift_len = 0.2
        detect_stream_dict = _prepare_data(
            family=self.party[0], detect_data=self.data, shift_len=shift_len)
        self.assertEqual(
            len(detect_stream_dict), len(self.party[0].detections))
        for detection_id, stream in detect_stream_dict.items():
            detection = [d for d in self.party[0] if d.id == detection_id][0]
            self.assertEqual(len(stream), detection.no_chans)
            for tr in stream:
                self.assertAlmostEqual(
                    tr.stats.endtime - tr.stats.starttime,
                    (2 * shift_len) + self.t_length, 1)

    # def test_prepare_data_too_short(self):
    #
    # def test_prepare_data_masked(self):
    #
    def test_family_picking(self):
        catalog = xcorr_pick_family(
            family=self.party[0], stream=self.data, shift_len=0.2, plot=True)

    # def test_family_picking_cccsum_reduce(self):
    #
    # def test_lag_calc_api(self):


class ShortTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        log = logging.getLogger(lag_calc.__name__)
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

        _xcorr_interp(ccc, 0.1)
        self.assertEqual(len(self.log_messages['warning']), 2)
        self.assertTrue(
            'Fewer than 5 samples' in self.log_messages['warning'][0])
        self.assertTrue(
            'Residual in quadratic fit' in self.log_messages['warning'][1])


if __name__ == '__main__':
    unittest.main()
