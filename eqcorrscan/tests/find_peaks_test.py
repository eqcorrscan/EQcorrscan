"""
Functions for testing the utils.findpeaks functions
"""
from os.path import join

import numpy as np
import pytest

from obspy.core.event import Catalog, Event, Origin

from eqcorrscan.utils.findpeaks import (
    find_peaks2_short, coin_trig, multi_find_peaks, find_peaks_compiled,
    _find_peaks_c, decluster, decluster_distance_time)
from eqcorrscan.utils.timer import time_func


class TestDeclustering:
    def test_unclustered_time(self):
        """ Check that nothing is lost if there aren't any clustered events."""
        peaks = np.array([100, 65, 20, 120, 300])
        index = np.array([2000, 5000, 10, 70, 500])
        trig_int = 30
        peaks_out = decluster(peaks, index, trig_int, threshold=0)
        assert len(peaks) == len(peaks_out)
        assert peaks_out == [(300.0, 500), (120.0, 70), (100.0, 2000),
                             (65.0, 5000), (20.0, 10)]

    def test_clustered_time(self):
        """ Check that the smallest is removed. """
        peaks = np.array([100, 65, 20, 120, 300])
        index = np.array([2000, 5000, 10, 70, 500])
        trig_int = 100
        peaks_out = decluster(peaks, index, trig_int, threshold=0)
        assert len(peaks) > len(peaks_out)
        assert peaks_out == [(300.0, 500), (120.0, 70), (100.0, 2000),
                             (65.0, 5000)]

    def test_clustered_time_longlong(self):
        """ Check that the smallest is removed when longlong func is used. """
        peaks = np.array([100, 65, 20, 120, 300], dtype=np.float32)
        index = np.array([2000, 5000, 10, 70, 500])
        index = index * int(1e10)
        trig_int = 100 * int(1e10)
        peaks_out = decluster(peaks, index, trig_int, threshold=0)
        assert len(peaks) > len(peaks_out)

    def test_clustered_dist_time(self):
        peaks = np.array([100, 65, 20, 120, 300])
        index = np.array([2000, 5000, 10, 70, 500])
        trig_int = 100
        hypocentral_separation = 10.0
        catalog = Catalog([
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
        ])
        peaks_out = decluster_distance_time(
            peaks, index, trig_int, catalog, hypocentral_separation,
            threshold=0)
        assert len(peaks) > len(peaks_out)
        assert peaks_out == [(300.0, 500), (120.0, 70), (100.0, 2000),
                             (65.0, 5000)]

    def test_separated_dist(self):
        peaks = np.array([100, 65, 20, 120, 300])
        index = np.array([2000, 5000, 10, 70, 500])
        trig_int = 100
        hypocentral_separation = 10.0
        catalog = Catalog([
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=80.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
        ])
        peaks_out = decluster_distance_time(
            peaks, index, trig_int, catalog, hypocentral_separation,
            threshold=0)
        assert len(peaks) == len(peaks_out)
        assert peaks_out == [(300.0, 500), (120.0, 70), (100.0, 2000),
                             (65.0, 5000), (20.0, 10)]

    def test_separated_dist_longlong(self):
        peaks = np.array([100, 65, 20, 120, 300])
        index = np.array([2000, 5000, 10, 70, 500])
        index = index * int(1e10)
        trig_int = 100 * int(1e10)
        hypocentral_separation = 10.0
        catalog = Catalog([
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=80.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
            Event(origins=[Origin(latitude=0.0, longitude=90.0, depth=1000.)]),
        ])
        peaks_out = decluster_distance_time(
            peaks, index, trig_int, catalog, hypocentral_separation,
            threshold=0)
        assert len(peaks) == len(peaks_out)


class TestStandardPeakFinding:
    """ Run peak finding against a standard cc array """
    trig_index = 100

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
    def peak_array_py(self, cc_array):
        """ run find_peaks2_short on cc_array and return results """
        peaks = find_peaks2_short(
            arr=cc_array, thresh=0.2, trig_int=self.trig_index)
        return peaks

    @pytest.fixture
    def full_peak_array_py(self, cc_array):
        """ run find_peaks2_short on cc_array and return results """
        peaks = find_peaks2_short(
            arr=cc_array, thresh=0.2, trig_int=self.trig_index,
            full_peaks=True)
        return peaks

    @pytest.fixture
    def peak_array_c(self, cc_array):
        """ run find_peaks2_short on cc_array and return results """
        peaks = find_peaks_compiled(
            arr=cc_array, thresh=0.2, trig_int=self.trig_index)
        return peaks

    @pytest.fixture
    def full_peak_array_c(self, cc_array):
        """ run find_peaks2_short on cc_array and return results """
        peaks = find_peaks_compiled(
            arr=cc_array, thresh=0.2, trig_int=self.trig_index,
            full_peaks=True)
        return peaks

    # tests
    def test_main_find_peaks(self, peak_array_c, peak_array_py,
                             expected_peak_array):
        """Test find_peaks2_short returns expected peaks """

        # Check length first as this will be a more obvious issue
        assert len(peak_array_c) == len(expected_peak_array), (
            'Peaks are not the same length, has ccc been updated?')
        assert (np.array(peak_array_c) == expected_peak_array).all()
        assert len(peak_array_py) == len(expected_peak_array), (
            'Peaks are not the same length, has ccc been updated?')
        assert (np.array(peak_array_py) == expected_peak_array).all()

    def test_find_all_peaks(self, full_peak_array_c, full_peak_array_py,
                            expected_peak_array):
        """Test finding all the peaks."""
        for peak in expected_peak_array[0]:
            assert peak in np.array(full_peak_array_c)
        assert len(full_peak_array_c) == 69
        for peak in expected_peak_array[0]:
            assert peak in np.array(full_peak_array_py)
        assert len(full_peak_array_py) == 69


class TestEdgeCases:
    """ A selection of weird datasets to find peaks in. """
    datasets = []
    data_len = 100000
    trig_int = 600

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets)
    def start_end_peaks(self):
        """ Array with peaks at first and last index. """
        end_arr = np.arange(self.data_len / 2)
        return np.concatenate([end_arr[::-1], end_arr]), \
            [0, self.data_len - 1], (self.data_len // 2) - 20

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets)
    def not_below_threshold(self):
        """ An array that doesn't drop below the threshold. """
        arr = np.ones(self.data_len, dtype=np.float32)
        spike_locs = np.random.randint(0, self.data_len, size=500)
        threshold = 0.5
        for spike_loc in spike_locs:
            arr[spike_loc] *= 100 * np.random.random()
        return arr, spike_locs, threshold

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets)
    def nearby_peaks(self):
        """ Peaks within trig-int of one-another - ensure highest is kept."""
        arr = np.ones(self.data_len, dtype=np.float32)
        biggest_peak_loc = self.data_len // 2
        arr[biggest_peak_loc] *= 1000
        arr[biggest_peak_loc - (self.trig_int // 2)] *= 900
        arr[biggest_peak_loc + (self.trig_int // 2)] *= 900
        return arr, [biggest_peak_loc], 10

    @pytest.fixture(scope='class', params=datasets)
    def dataset(self, request):
        return request.getfixturevalue(request.param)

    def test_python_speed(self, dataset, request):
        """ test that findpeaks works on each of the arrays and print the
        time it took to run. Tests passes if no error is raised"""
        print('starting find_peak profiling on: ' + request.node.name)
        threshold = dataset[2]
        print("Threshold: {0}".format(threshold))
        # Run the prototype and check that is gets the same results!
        proto_peaks = time_func(
            find_peaks_compiled, "find_peaks_compiled", arr=dataset[0],
            thresh=threshold, trig_int=self.trig_int)
        print('Found %i peaks' % len(proto_peaks))
        for peak in proto_peaks:
            assert abs(peak[0]) > threshold
            assert np.allclose(peak[0], dataset[0][peak[1]], atol=0.001)
            assert peak[1] in dataset[1]
        assert len(proto_peaks) <= len(dataset[1])
        peaks = time_func(
            find_peaks2_short, "find_peaks2_short", arr=dataset[0],
            thresh=threshold, trig_int=self.trig_int, full_peaks=True)
        print('Found %i peaks' % len(peaks))
        assert len(peaks) <= len(dataset[1])
        for peak in peaks:
            assert abs(peak[0]) > threshold  # Assert peak is above threshold
            assert peak[0] == dataset[0][peak[1]]
            assert peak[1] in dataset[1]  # Assert peak in expected peaks
            # Assert that the correct value is given!
        assert len(proto_peaks) <= len(dataset[1])
        # The C func does a better job
        # assert np.allclose(peaks, proto_peaks, atol=0.001)


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
    data_len = 100000
    n_channels = 50  # Number of channels for noisy 2D dataset
    DEBUG = False

    # fixtures that create 1D datasets
    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def random(self):
        """ simple random array """
        arr = np.random.randn(self.data_len)
        return arr, [], np.abs(arr).max() + 10

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def noisy(self):
        """ noisy array """
        arr = np.random.randn(self.data_len) ** 5
        return arr, [], np.abs(arr).max() + 10

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def spiky(self):
        """ array with large spikes """
        arr = np.random.randn(self.data_len) ** 5
        spike_locs = np.random.randint(0, self.data_len, size=500)
        threshold = np.abs(arr).max() + 20
        for spike_loc in spike_locs:
            arr[spike_loc] *= 1000
        return arr, spike_locs, threshold

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def all_above_threshold(self):
        """ array with large spikes """
        arr = np.ones(self.data_len, dtype=float)
        spike_locs = np.random.randint(0, self.data_len, size=500)
        threshold = 0.5
        for spike_loc in spike_locs:
            arr[spike_loc] = 10 * spike_loc
            # Deliberately make each peak different height. When all peaks
            # are the same height C and Python return different (but both
            # valid) peaks).
        return arr, spike_locs, threshold

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_1d)
    def clustered(self):
        arr = np.random.randn(self.data_len)
        spike_locs = np.random.randint(0, self.data_len / 10, size=200)
        spike_locs = np.append(
            spike_locs, np.random.randint(
                2 * (self.data_len / 10), 4 * (self.data_len / 10), size=400))
        threshold = np.abs(arr).max() + 20
        for spike_loc in spike_locs:
            arr[spike_loc] *= 1000
        return arr, spike_locs, threshold

    @pytest.fixture(scope='class', params=datasets_1d)
    def dataset_1d(self, request):
        """ parametrize the 1d datasets """
        return request.getfixturevalue(request.param)

    # fixtures that create 2D datasets
    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_2d)
    def aggregated_1d_datasets(self, request):
        """ create a 2d numpy array of all the datasets """
        return np.array([request.getfixturevalue(x)[0]
                         for x in self.datasets_1d])

    @pytest.fixture(scope='class')
    @pytest.append_name(datasets_2d)
    def noisy_multi_array(self):
        """ a noisy 2d array """
        return np.random.randn(self.n_channels, self.data_len) ** 5

    @pytest.fixture(scope='class', params=datasets_2d)
    def dataset_2d(self, request):
        """ parametrize the 2d datasets """
        return request.getfixturevalue(request.param)

    # tests
    def test_python_speed(self, dataset_1d, request):
        """ test that findpeaks works on each of the arrays and print the
        time it took to run. Tests passes if no error is raised"""
        print('starting find_peak profiling on: ' + request.node.name)
        threshold = dataset_1d[2]
        print("Threshold: {0}".format(threshold))
        peaks = time_func(
            find_peaks2_short, "find_peaks2_short", arr=dataset_1d[0],
            thresh=threshold, trig_int=600)
        print('Found %i peaks' % len(peaks))
        assert len(peaks) <= len(dataset_1d[1])
        for peak in peaks:
            assert abs(peak[0]) > threshold
            assert peak[1] in dataset_1d[1]
        # Run the prototype and check that is gets the same results!
        proto_peaks = time_func(
            find_peaks_compiled, "find_peaks_compiled", arr=dataset_1d[0],
            thresh=threshold, trig_int=600)
        print('Found %i peaks' % len(proto_peaks))
        for peak in proto_peaks:
            assert abs(peak[0]) > threshold
            assert peak[1] in dataset_1d[1]
        assert len(proto_peaks) <= len(dataset_1d[1])
        # The C func does a better job
        # assert np.allclose(peaks, proto_peaks, atol=0.001)

    def test_multi_find_peaks(self, dataset_2d, request):
        """ ensure the same results are returned for serial and parallel
        in multi_find_peaks """
        print('starting find_peak profiling on: ' + request.node.name)
        arr = dataset_2d.astype(np.float32)
        threshold = [np.float32(10 * np.median(np.abs(x))) for x in arr]
        print("Running serial C loop")
        serial_c_peaks = time_func(
            multi_find_peaks, name="serial-C", arr=arr, thresh=threshold,
            trig_int=600, parallel=False, internal_func=find_peaks_compiled)
        print("Running parallel C loop")
        parallel_c_peaks = time_func(
            multi_find_peaks, name="parallel-C", arr=arr, thresh=threshold,
            trig_int=600, parallel=True, internal_func=find_peaks_compiled)

        print("Running serial Python loop")
        serial_py_peaks = time_func(
            multi_find_peaks, name="serial-Python", arr=arr, thresh=threshold,
            trig_int=600, parallel=False, internal_func=find_peaks2_short,
            full_peaks=True)
        print("Running parallel Python loop")
        parallel_py_peaks = time_func(
            multi_find_peaks, name="parallel-Python", arr=arr,
            thresh=threshold, trig_int=600, parallel=True,
            internal_func=find_peaks2_short, full_peaks=True)
        assert serial_py_peaks == parallel_py_peaks
        if not serial_c_peaks == parallel_c_peaks and self.DEBUG:
            for _serial_c_peaks, _parallel_c_peaks in zip(
                    serial_c_peaks, parallel_c_peaks):
                for peak in _serial_c_peaks:
                    if peak not in _parallel_c_peaks:
                        print("Peak in serial but not parallel: {0}".format(
                            peak))
                for peak in _parallel_c_peaks:
                    if peak not in _serial_c_peaks:
                        print("Peak in parallel but not serial: {0}".format(
                            peak))
            # Test the first step
            parallel_peak_vals, parallel_peak_indices = multi_find_peaks(
                arr=arr, thresh=threshold, cores=2,
                internal_func=find_peaks_compiled)
            parallel_sorted = []
            parallel_peaks_sorted = []
            parallel_indices_sorted = []
            for _peaks, _indices in zip(parallel_peak_vals,
                                        parallel_peak_indices):
                if len(_peaks) == 0:
                    parallel_sorted.append([])
                    continue
                _peaks_sort = sorted(
                    zip(_peaks, _indices),
                    key=lambda amplitude: abs(amplitude[0]), reverse=True)
                parallel_sorted.append(_peaks_sort)
                _arr, _ind = zip(*_peaks_sort)
                parallel_peaks_sorted.extend(_arr)
                parallel_indices_sorted.extend(_ind)
            serial_peak_vals, serial_peak_indices, serial_sorted = ([], [], [])
            for sub_arr, thresh in zip(arr, threshold):
                _peak_vals, _peak_indices = _find_peaks_c(
                    array=sub_arr, threshold=thresh)
                serial_peak_vals.append(_peak_vals)
                serial_peak_indices.append(_peak_indices)
                _peaks_sort = sorted(
                    zip(_peak_vals, _peak_indices),
                    key=lambda amplitude: abs(amplitude[0]), reverse=True)
                serial_sorted.append(_peaks_sort)
            for i in range(len(serial_peak_vals)):
                if parallel_peak_vals[i].size > 0:
                    assert np.all(parallel_peak_vals[i] == serial_peak_vals[i])
                    assert np.all(parallel_peak_indices[i] ==
                                  serial_peak_indices[i])
                    assert parallel_sorted == serial_sorted
                else:
                    assert serial_peak_vals[i].size == 0
            np.save("test_2d_array.npy", arr)
        assert serial_c_peaks == parallel_c_peaks
        for i in range(len(serial_c_peaks)):
            diff_count = 0
            for j in range(len(serial_c_peaks[i])):
                if not serial_c_peaks[i][j] in serial_py_peaks[i]:
                    diff_count += 1
                    print("Peak {0} in C but not in py".format(
                        serial_c_peaks[i][j]))
            for j in range(len(serial_py_peaks[i])):
                if not serial_py_peaks[i][j] in serial_c_peaks[i]:
                    diff_count += 1
                    print("Peak {0} in py but not in C".format(
                        serial_py_peaks[i][j]))
            assert diff_count <= 0.0001 * self.data_len
        if self.DEBUG:
            np.save("test_2d_array.npy", arr)
            np.save("test_c_peaks_serial.npy", serial_c_peaks)
            np.save("test_c_peaks_parallel.npy", parallel_c_peaks)
            np.save("test_py_peaks_serial.npy", serial_py_peaks)
            np.save("test_py_peaks_parallel.npy", parallel_py_peaks)

    def test_noisy_timings(self, noisy_multi_array):
        arr = noisy_multi_array.astype(np.float32)
        threshold = [np.float32(np.median(np.abs(d)))
                     for d in noisy_multi_array]
        print("Running serial loop")
        serial_peaks = time_func(
            multi_find_peaks, name="serial", arr=arr,
            thresh=threshold, trig_int=600, parallel=False)
        print("Running parallel loop")
        parallel_peaks = time_func(
            multi_find_peaks, name="parallel", arr=arr,
            thresh=threshold, trig_int=600, parallel=True)
        assert len(serial_peaks) == len(parallel_peaks)
        for i in range(len(serial_peaks)):
            assert serial_peaks[i] == parallel_peaks[i]
