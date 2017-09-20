"""
Functions for testing the utils.findpeaks functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest

from os.path import join

import numpy as np
import pytest

from eqcorrscan.utils.timer import time_func
from eqcorrscan.utils.findpeaks import find_peaks2_short, coin_trig


class TestStandardPeakFinding:
    """ Run peak finding against a standard cc array """
    trig_index = 10

    # fixtures
    @pytest.fixture
    def cc_array(self):
        """ load the test cc array case """
        return np.load(join(pytest.test_data_path, 'test_ccc.npy'))

    @pytest.fixture
    def expected_peak_array(self):
        """ load the peak array from running find peaks on cc_array """
        return np.load(join(pytest.test_data_path, 'test_peaks.npy'))

    @pytest.fixture
    def peak_array(self, cc_array):
        """ run find_peaks2_short on cc_array and return results """
        peaks = find_peaks2_short(
            arr=cc_array, thresh=0.2, trig_int=self.trig_index,
            debug=0, starttime=None, samp_rate=200.0)
        return peaks

    @pytest.fixture
    def full_peak_array(self, cc_array):
        """ run find_peaks2_short on cc_array and return results """
        peaks = find_peaks2_short(
            arr=cc_array, thresh=0.2, trig_int=self.trig_index,
            debug=0, starttime=None, samp_rate=200.0, full_peaks=True)
        return peaks

    # tests
    # def test_max_values(self, old_peak_array, cc_array):
    #     """
    #     ensure the values in the peak array are max values in
    #     the expected window.
    #
    #     This doesn't really work - for the case where peak a is at position
    #     0 and peak b is at position 11, with trig_int = 10, there may be
    #     values higher than peak b in positions 1-10, but these are ruled out
    #     by peak a... Hence test fails.
    #     """
    #     abs_array = np.abs(cc_array)
    #     for val, ind in old_peak_array:
    #         assert val == cc_array[ind]
    #         start, stop = ind - self.trig_index, ind + self.trig_index + 1
    #         window = abs_array[start: stop]
    #         assert np.max(window) == np.abs(val)

    def test_main_find_peaks(self, peak_array, expected_peak_array):
        """Test find_peaks2_short returns expected peaks """

        # Check length first as this will be a more obvious issue
        assert len(peak_array) == len(expected_peak_array), (
            'Peaks are not the same length, has ccc been updated?')
        assert (np.array(peak_array) == expected_peak_array).all()

    def test_find_all_peaks(self, full_peak_array, expected_peak_array):
        """Test finding all the peaks."""
        for peak in expected_peak_array[0]:
            assert peak in np.array(full_peak_array)
        assert len(full_peak_array) == 287


class TestCoincidenceTrigger:
    # fixtures
    @pytest.fixture
    def peaks(self):
        """" create a sample peak array """
        peaks = [[(0.5, 100), (0.3, 800), (0.3, 105)],
                 [(0.4, 120), (0.7, 850)]]
        return peaks

    # tests
    def test_coincidence(self, peaks):
        """Test the coincidence trigger."""

        triggers = coin_trig(peaks, [('a', 'Z'), ('b', 'Z')], samp_rate=10,
                             moveout=3, min_trig=2, trig_int=1)
        assert triggers, [(0.45, 100)]


class TestPeakFindSpeeds:
    @pytest.fixture
    def datasets(self):
        data_len = 1000000
        datasets = {'random': np.random.randn(data_len),
                    'noisy': np.random.randn(data_len) ** 5,
                    'spiky': np.random.randn(data_len),
                    'clustered': np.random.randn(data_len)}
        spike_locs = np.random.randint(0, data_len, size=500)
        for spike_loc in spike_locs:
            datasets['spiky'][spike_loc] *= 1000
        spike_locs = np.random.randint(0, data_len / 10, size=200)
        spike_locs = np.append(
            spike_locs, np.random.randint(
                2 * (data_len / 10), 4 * (data_len / 10), size=400))
        for spike_loc in spike_locs:
            datasets['clustered'][spike_loc] *= 1000
        return datasets

    def test_python_speed(self, datasets):
        print("Running Python declustering")
        for key in datasets.keys():
            print("\tRunning timings for %s dataset" % key)
            mad_thresh = 10 * np.median(np.abs(datasets[key]))
            peaks = time_func(
                find_peaks2_short, "find_peaks2_short", arr=datasets[key],
                thresh=mad_thresh, trig_int=600, debug=0, starttime=False,
                samp_rate=100)
            print('Found %i peaks' % len(peaks))
