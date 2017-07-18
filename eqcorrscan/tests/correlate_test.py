from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest
import numpy as np
import time

from eqcorrscan.utils.correlate import scipy_normxcorr, fftw_xcorr_2d, \
    fftw_xcorr, time_multi_normxcorr


class CorrelateTests(unittest.TestCase):
    def test_same_various_methods(self):
        templates = np.random.randn(200, 200)
        stream = np.random.randn(10000)
        stream *= stream ** 10
        pads = np.zeros(templates.shape[0], dtype=int)
        tic = time.time()
        scipy_ccc, no_chans = scipy_normxcorr(templates, stream, pads)
        toc = time.time()
        print('Scipy took: %f seconds' % (toc-tic))
        # tic = time.time()
        # fftw_ccc, no_chans = fftw_xcorr(templates, stream, pads)
        # toc = time.time()
        # print('FFTW took: %f seconds' % (toc-tic))
        tic = time.time()
        fftw_ccc_2d, no_chans = fftw_xcorr_2d(templates, stream, pads)
        toc = time.time()
        print('FFTW-2D took: %f seconds' % (toc-tic))
        # tic = time.time()
        # time_ccc, no_chans = time_multi_normxcorr(templates, stream, pads)
        # toc = time.time()
        # print('Time-domain took: %f seconds' % (toc-tic))
        # self.assertTrue(np.allclose(scipy_ccc, fftw_ccc, atol=0.001))
        self.assertTrue(np.allclose(scipy_ccc, fftw_ccc_2d, atol=0.001))
        # self.assertTrue(np.allclose(scipy_ccc, time_ccc, atol=0.1))
        # self.assertTrue(np.allclose(fftw_ccc, time_ccc, atol=0.1))


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
