"""
Functions to test the functions within the eqcorrscan.utils.catalogue2DD.py \
submodule.  Uses test data distributed with the EQcorrscan pacakge.
"""
import unittest


class TestCatalogueMethods(unittest.TestCase):
    def test_rounding(self):
        """Simple test to test that _cc_round gives correct result.
        """
        from eqcorrscan.utils.catalogue2DD import _cc_round
        import numpy as np
        # Use an irrational number and round it to various decimal places
        test_no = np.pi
        self.assertEqual(_cc_round(test_no, 0), '3')
        self.assertEqual(_cc_round(test_no, 1), '3.1')
        self.assertEqual(_cc_round(test_no, 2), '3.14')
        self.assertEqual(_cc_round(test_no, 3), '3.142')
        self.assertEqual(_cc_round(test_no, 4), '3.1416')
        self.assertEqual(_cc_round(test_no, 5), '3.14159')
        self.assertEqual(_cc_round(test_no, 6), '3.141593')
        self.assertEqual(_cc_round(test_no, 7), '3.1415927')
        self.assertEqual(_cc_round(test_no, 8), '3.14159265')
        self.assertEqual(_cc_round(test_no, 9), '3.141592654')
        self.assertEqual(_cc_round(test_no, 10), '3.1415926536')
        self.assertEqual(_cc_round(test_no, 11), '3.14159265359')
        self.assertEqual(_cc_round(test_no, 12), '3.141592653590')
        self.assertEqual(_cc_round(test_no, 13), '3.1415926535898')
        self.assertEqual(_cc_round(test_no, 14), '3.14159265358979')
        self.assertEqual(_cc_round(test_no, 15), '3.141592653589793')

    def test_weight_averaging(self):
        """Simple function to test _av_weight returns the correct weights.
        """
        from eqcorrscan.utils.catalogue2DD import _av_weight
        self.assertEqual(_av_weight('0', '0'), '1.0000')
        self.assertEqual(_av_weight('0', '1'), '0.8750')
        self.assertEqual(_av_weight('0', '2'), '0.7500')
        self.assertEqual(_av_weight('0', '3'), '0.6250')
        self.assertEqual(_av_weight('0', '4'), '0.5000')
        self.assertEqual(_av_weight('0', '9'), '0.5000')
        self.assertEqual(_av_weight('1', '0'), '0.8750')
        self.assertEqual(_av_weight('1', '1'), '0.7500')
        self.assertEqual(_av_weight('1', '2'), '0.6250')
        self.assertEqual(_av_weight('1', '3'), '0.5000')
        self.assertEqual(_av_weight('1', '4'), '0.3750')
        self.assertEqual(_av_weight('1', '9'), '0.3750')
        self.assertEqual(_av_weight('2', '0'), '0.7500')
        self.assertEqual(_av_weight('2', '1'), '0.6250')
        self.assertEqual(_av_weight('2', '2'), '0.5000')
        self.assertEqual(_av_weight('2', '3'), '0.3750')
        self.assertEqual(_av_weight('2', '4'), '0.2500')
        self.assertEqual(_av_weight('2', '9'), '0.2500')
        self.assertEqual(_av_weight('3', '0'), '0.6250')
        self.assertEqual(_av_weight('3', '1'), '0.5000')
        self.assertEqual(_av_weight('3', '2'), '0.3750')
        self.assertEqual(_av_weight('3', '3'), '0.2500')
        self.assertEqual(_av_weight('3', '4'), '0.1250')
        self.assertEqual(_av_weight('3', '9'), '0.1250')
        self.assertEqual(_av_weight('4', '0'), '0.5000')
        self.assertEqual(_av_weight('4', '1'), '0.3750')
        self.assertEqual(_av_weight('4', '2'), '0.2500')
        self.assertEqual(_av_weight('4', '3'), '0.1250')
        self.assertEqual(_av_weight('4', '4'), '0.0000')
        self.assertEqual(_av_weight('4', '9'), '0.0000')

    def test_readSTATION0(self):
        """Simple function to test the ability to read a test STATION0.HYP \
        file."""
        from eqcorrscan.utils.catalogue2DD import readSTATION0
        import os
        station_input_list = []
        STATION0_path = os.path.join(os.path.abspath(os.path.
                                                     dirname(__file__)),
                                     'test_data')
        station_output_list = readSTATION0(STATION0_path, station_input_list)
        # Check that the output file exists, and remove it
        self.assertTrue(os.path.isfile('station.dat'))
        os.remove('station.dat')


if __name__ == '__main__':
    unittest.main()
