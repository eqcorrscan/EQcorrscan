"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.core import sliding_normxcorr
import unittest
import numpy as np


class TestSliding(unittest.TestCase):
    """Test the sliding window correlation."""
    def test_sine(self):
        arr1 = np.sin(np.arange(10000)/200.0)
        arr2 = np.sin(np.arange(10000)/200.0+10)
        shift, ccs = sliding_normxcorr.\
            sliding_normxcorr(arr1.astype(np.float32),
                              arr2.astype(np.float32),
                              600)
        self.assertEqual(round(ccs, 5), 1.0)
        self.assertEqual(shift, -513)  # Interpolation rounding

    def test_fail(self):
        arr1 = np.sin(np.arange(10000)/200.0)
        arr2 = np.sin(np.arange(10000)/200.0+10)
        with self.assertRaises(IndexError):
            shift, ccs = sliding_normxcorr.\
                sliding_normxcorr(arr1.astype(np.float32),
                                  arr2.astype(np.float32),
                                  100000)

    def test_full(self):
        arr1 = np.sin(np.arange(10000)/200.0)
        arr2 = np.sin(np.arange(10000)/200.0+10)
        shift, ccs, xcorr = sliding_normxcorr.\
            sliding_normxcorr(arr1.astype(np.float32),
                              arr2.astype(np.float32),
                              10, full_xcorr=True)
        self.assertEqual(len(xcorr), 21)


if __name__ == '__main__':
    """
    Run stacking tests
    """
    unittest.main()