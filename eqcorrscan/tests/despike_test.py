"""Test the despiking algorithms implimented in EQcorrscan."""

import unittest
import os
import numpy as np

from obspy import Trace, read

from eqcorrscan.utils.despike import median_filter, template_remove


class DespikeTesting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        cls.spiked = read(
            os.path.join(cls.testing_path, 'random_spiked.ms'))[0]

    def test_median_filter(self):
        """Test the median filter implementation."""
        despiked = median_filter(tr=self.spiked.copy(), multiplier=2,
                                 windowlength=0.5, interp_len=0.05)
        self.assertNotEqual(despiked.data[100], 20)
        self.assertNotEqual(despiked.data[400], 40)
        self.assertNotEqual(despiked.data[450], -40)

    def test_template_remove(self):
        """Test the despiker based on correlations."""
        template = np.zeros(10)
        template[2] = 1
        template = Trace(template)
        template.stats.sampling_rate = 100
        despiked = template_remove(
            tr=self.spiked.copy(), template=template, cc_thresh=0.9,
            windowlength=0.05, interp_len=0.05)
        self.assertNotEqual(despiked.data[100], 20)
        self.assertNotEqual(despiked.data[400], 40)
        self.assertNotEqual(despiked.data[450], -40)


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
