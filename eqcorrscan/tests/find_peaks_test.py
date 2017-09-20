"""
Functions for testing the utils.findpeaks functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest

import numpy as np
import os

from eqcorrscan.utils.timer import time_func
from eqcorrscan.utils.findpeaks import find_peaks2_short, coin_trig


#TODO: Write a test for multi-peak finding and time the functions!
class TestPeakFinding(unittest.TestCase):
    def test_main_find_peaks(self):
        """Test find_peaks2_short"""
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        expected_ccc = np.load(os.path.join(testing_path, 'test_ccc.npy'))
        peaks = find_peaks2_short(arr=expected_ccc, thresh=0.2, trig_int=10,
                                  debug=0, starttime=False, samp_rate=200.0,
                                  compiled=False)
        peaks_c = find_peaks2_short(arr=expected_ccc, thresh=0.2, trig_int=10,
                                    debug=0, starttime=False, samp_rate=200.0,
                                    compiled=True)
        expected_peaks = np.load(os.path.join(testing_path, 'test_peaks.npy'))
        # Check length first as this will be a more obvious issue
        self.assertEqual(len(peaks), len(expected_peaks),
                         msg='Peaks are not the same length, has ccc been ' +
                         'updated?')
        self.assertTrue((np.array(peaks) == expected_peaks).all())
        self.assertEqual(len(peaks_c), len(expected_peaks),
                         msg='Peaks are not the same length, has ccc been ' +
                         'updated?')
        self.assertTrue((np.array(peaks_c) == expected_peaks).all())

    def test_coincidence(self):
        """Test the coincidence trigger."""
        peaks = [[(0.5, 100), (0.3, 800), (0.3, 105)],
                 [(0.4, 120), (0.7, 850)]]
        triggers = coin_trig(peaks, [('a', 'Z'), ('b', 'Z')], samp_rate=10,
                             moveout=3, min_trig=2, trig_int=1)
        self.assertEqual(triggers, [(0.45, 100)])


class TestPeakFindSpeeds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data_len = 1000000
        cls.datasets = {'random': np.random.randn(data_len),
                        'noisy': np.random.randn(data_len) ** 5,
                        'spiky': np.random.randn(data_len),
                        'clustered': np.random.randn(data_len)}
        spike_locs = np.random.randint(0, data_len, size=500)
        for spike_loc in spike_locs:
            cls.datasets['spiky'][spike_loc] *= 1000
        spike_locs = np.random.randint(0, data_len / 10, size=200)
        spike_locs = np.append(
            spike_locs, np.random.randint(
                2 * (data_len / 10), 4 * (data_len / 10), size=400))
        for spike_loc in spike_locs:
            cls.datasets['clustered'][spike_loc] *= 1000

    def test_python_speed(self):
        print("Running Python declustering")
        for key in self.datasets.keys():
            print("\tRunning timings for %s dataset" % key)
            mad_thresh = 10 * np.median(np.abs(self.datasets[key]))
            peaks = time_func(
                find_peaks2_short, "find_peaks2_short", arr=self.datasets[key],
                thresh=mad_thresh, trig_int=600, debug=0, starttime=False,
                samp_rate=100)
            print('Found %i peaks' % len(peaks))

    def test_c_speed(self):
        print("Running compiled declustering")
        for key in self.datasets.keys():
            print("\tRunning timings for %s dataset" % key)
            mad_thresh = 10 * np.median(np.abs(self.datasets[key]))
            peaks = time_func(
                find_peaks2_short, "find_peaks2_short", arr=self.datasets[key],
                thresh=mad_thresh, trig_int=600, debug=0, starttime=False,
                samp_rate=100, compiled=True)
            print('Found %i peaks' % len(peaks))


if __name__ == '__main__':
    """
    Run tests
    """
    unittest.main()
