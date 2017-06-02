"""
Functions for testing the utils.findpeaks functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestPeakFinding(unittest.TestCase):
    def test_main_find_peaks(self):
        """Test find_peaks2_short"""
        from eqcorrscan.utils.findpeaks import find_peaks2_short
        import numpy as np
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        expected_ccc = np.load(os.path.join(testing_path, 'test_ccc.npy'))
        peaks = find_peaks2_short(arr=expected_ccc, thresh=0.2, trig_int=10,
                                  debug=0, starttime=False, samp_rate=200.0)
        expected_peaks = np.load(os.path.join(testing_path, 'test_peaks.npy'))
        # Check length first as this will be a more obvious issue
        self.assertEqual(len(peaks), len(expected_peaks),
                         msg='Peaks are not the same length, has ccc been ' +
                         'updated?')
        self.assertTrue((np.array(peaks) == expected_peaks).all())

    def test_coincidence(self):
        """Test the coincidence trigger."""
        from eqcorrscan.utils.findpeaks import coin_trig
        peaks = [[(0.5, 100), (0.3, 800), (0.3, 105)],
                 [(0.4, 120), (0.7, 850)]]
        triggers = coin_trig(peaks, [('a', 'Z'), ('b', 'Z')], samp_rate=10,
                             moveout=3, min_trig=2, trig_int=1)
        self.assertEqual(triggers, [(0.45, 100)])

if __name__ == '__main__':
    """
    Run tests
    """
    unittest.main()
