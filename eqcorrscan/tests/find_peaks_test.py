"""
Functions for testing the utils.findpeaks functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestTutorialScripts(unittest.TestCase):
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
        self.assertTrue((np.array(peaks) == expected_peaks).all())

if __name__ == '__main__':
    """
    Run tests
    """
    unittest.main()
