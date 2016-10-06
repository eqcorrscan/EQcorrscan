"""
A series of test functions for the core functions in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import os
import numpy as np

from eqcorrscan.core import lag_calc
from eqcorrscan.core.template_gen import from_sfile
from eqcorrscan.core.match_filter import normxcorr2


class TestMethods(unittest.TestCase):
    def setUp(self):
        self.testing_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data', 'REA', 'TEST_')
        self.template = from_sfile(
            sfile=os.path.join(self.testing_path, '21-1412-02L.S201309'),
            lowcut=5, highcut=15, samp_rate=40, filt_order=4, length=3,
            swin='all', prepick=0.05)
        self.detection = from_sfile(
            sfile=os.path.join(self.testing_path, '21-1759-04L.S201309'),
            lowcut=5, highcut=15, samp_rate=40, filt_order=4, length=4,
            swin='all', prepick=0.55)

    def test_channel_loop(self):
        """Test the main lag_calc function"""
        i, event = lag_calc._channel_loop(detection=self.detection,
                                          template=self.template,
                                          min_cc=0.4, i=0,
                                          detection_id='Tester_01',
                                          interpolate=False)
        matched_traces = []
        detection_stachans = [(tr.stats.station, tr.stats.channel)
                              for tr in self.detection]
        picked_stachans = [(pick.waveform_id.station_code,
                            pick.waveform_id.channel_code)
                           for pick in event.picks]
        for master_tr in self.template:
            stachan = (master_tr.stats.station, master_tr.stats.channel)
            if stachan in detection_stachans:
                matched_traces.append(stachan)

        for picked_stachan in picked_stachans:
            self.assertTrue(picked_stachan in matched_traces)

    def test_interpolate(self):
        """Test channel loop with interpolation."""
        i, event = lag_calc._channel_loop(detection=self.detection,
                                          template=self.template,
                                          min_cc=0.4, i=0,
                                          detection_id='Tester_01',
                                          interpolate=True)
        matched_traces = []
        detection_stachans = [(tr.stats.station, tr.stats.channel)
                              for tr in self.detection]
        picked_stachans = [(pick.waveform_id.station_code,
                            pick.waveform_id.channel_code)
                           for pick in event.picks]
        for master_tr in self.template:
            stachan = (master_tr.stats.station, master_tr.stats.channel)
            if stachan in detection_stachans:
                matched_traces.append(stachan)

        for picked_stachan in picked_stachans:
            self.assertTrue(picked_stachan in matched_traces)

    def test_interp_normal(self):
        synth_template = np.sin(np.arange(0, 2, 0.001))
        synth_detection = synth_template[100:]
        synth_template = synth_template[0:-10]
        ccc = normxcorr2(synth_detection, synth_template)
        shift, coeff = lag_calc._xcorr_interp(ccc, 0.01)
        self.assertEqual(shift.round(), 1.0)
        self.assertEqual(coeff.round(), 1.0)

    def test_interp_few_samples(self):
        synth_template = np.sin(np.arange(0, 2, 0.001))
        synth_detection = synth_template[13:]
        synth_template = synth_template[0:-10]
        ccc = normxcorr2(synth_detection, synth_template)
        shift, coeff = lag_calc._xcorr_interp(ccc, 0.01)
        self.assertEqual(shift.round(), 0.0)
        self.assertEqual(coeff.round(), 1.0)

    def test_interp_not_enough_samples(self):
        synth_template = np.sin(np.arange(0, 2, 0.001))
        synth_detection = synth_template[11:]
        synth_template = synth_template[0:-10]
        ccc = normxcorr2(synth_detection, synth_template)[0]
        with self.assertRaises(IndexError):
            lag_calc._xcorr_interp(ccc, 0.01)


if __name__ == '__main__':
    unittest.main()
