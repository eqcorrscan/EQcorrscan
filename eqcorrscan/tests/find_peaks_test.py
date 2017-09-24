"""
Functions for testing the utils.findpeaks functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from os.path import join

import numpy as np
import pytest

from eqcorrscan.utils.findpeaks import (
    find_peaks2_short, coin_trig, multi_find_peaks)
from eqcorrscan.utils.timer import time_func


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
        assert len(full_peak_array) == 315


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


@pytest.mark.serial
class TestPeakFindSpeeds:
    """ test findpeaks on various themes of arrays """
    datasets_1d = []
    datasets_2d = []
    data_len = 1000000

    # fixtures that create 1D datasets
    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def random(self):
        """ simple random array """
        return np.random.randn(self.data_len)

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def noisy(self):
        """ noisy array """
        return np.random.randn(self.data_len) ** 5

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def spiky(self):
        """ array with large spikes """
        arr = np.random.randn(self.data_len) ** 5
        spike_locs = np.random.randint(0, self.data_len, size=500)
        for spike_loc in spike_locs:
            arr[spike_loc] *= 1000
        return arr

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def clustered(self):
        arr = np.random.randn(self.data_len)
        spike_locs = np.random.randint(0, self.data_len / 10, size=200)
        spike_locs = np.append(
            spike_locs, np.random.randint(
                2 * (self.data_len / 10), 4 * (self.data_len / 10), size=400))
        for spike_loc in spike_locs:
            arr[spike_loc] *= 1000
        return arr

    @pytest.fixture(scope='class', params=datasets_1d)
    def dataset_1d(self, request):
        """ parametrize the 1d datasets """
        return request.getfuncargvalue(request.param)

    # fixtures that create 2D datasets
    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_2d)
    def aggregated_1d_datasets(self, request):
        """ create a 2d numpy array of all the datasets """
        return np.array([request.getfuncargvalue(x)
                         for x in self.datasets_1d])

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_2d)
    def noisy_multi_array(self):
        """ a noisy 2d array """
        return np.random.randn(4, self.data_len) ** 5

    @pytest.fixture(scope='class', params=datasets_2d)
    def dataset_2d(self, request):
        """ parametrize the 2d datasets """
        return request.getfuncargvalue(request.param)

    # tests
    def test_python_speed(self, dataset_1d, request):
        """ test that findpeaks works on each of the arrays and print the
        time it took to run. Tests passes if no error is raised"""
        print('starting find_peak profiling on: ' + request.node.name)
        threshold = np.median(dataset_1d)  # get threshold
        peaks = time_func(
            find_peaks2_short, "find_peaks2_short", arr=dataset_1d,
            thresh=threshold, trig_int=600, debug=0, starttime=False,
            samp_rate=100)
        print('Found %i peaks' % len(peaks))

    def test_multi_find_peaks(self, dataset_2d, request):
        """ ensure the same results are returned for serial and parallel
        in multi_find_peaks """
        print('starting find_peak profiling on: ' + request.node.name)
        arr = dataset_2d
        threshold = [10 * np.median(np.abs(x)) for x in dataset_2d]
        print("Running serial loop")
        serial_peaks = time_func(
            multi_find_peaks, name="serial", arr=arr, thresh=threshold,
            trig_int=600, parallel=False)
        print("Running parallel loop")
        parallel_peaks = time_func(
            multi_find_peaks, name="parallel", arr=arr, thresh=threshold,
            trig_int=600, parallel=True)
        assert serial_peaks == parallel_peaks

    def test_noisy_timings(self, noisy_multi_array):
        threshold = [np.median(np.abs(d)) for d in noisy_multi_array]
        print("Running serial loop")
        serial_peaks = time_func(
            multi_find_peaks, name="serial", arr=noisy_multi_array,
            thresh=threshold, trig_int=600, parallel=False)
        print("Running parallel loop")
        parallel_peaks = time_func(
            multi_find_peaks, name="parallel", arr=noisy_multi_array,
            thresh=threshold, trig_int=600, parallel=True)
        assert serial_peaks == parallel_peaks
