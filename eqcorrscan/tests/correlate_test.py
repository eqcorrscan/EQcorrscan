from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest
import numpy as np
from eqcorrscan.utils.correlate import *


class CorrelateTests(unittest.TestCase):
    def test_same_various_methods(self):
        templates = np.random.randn(5, 200)
        stream = np.random.randn(10000)
        pads = [0, 0, 0, 0, 0]
        scipy_ccc, no_chans= scipy_normxcorr(templates, stream, pads)
        fftw_ccc, no_chans = fftw_xcorr(templates, stream, pads)
        time_ccc, no_chans = time_multi_normxcorr(templates, stream, pads)
        self.assertTrue(np.allclose(scipy_ccc, fftw_ccc, atol=0.5))
        self.assertTrue(np.allclose(scipy_ccc, time_ccc, atol=0.5))
        self.assertTrue(np.allclose(fftw_ccc, time_ccc, atol=0.5))


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()