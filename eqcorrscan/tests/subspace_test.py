"""
Functions for testing the core.subspace functions
"""
import numpy as np
import unittest
import pytest
import os
import copy

from obspy import Stream, read

from eqcorrscan.core import subspace


class SimpleSubspaceMethods(unittest.TestCase):
    """
    Tests that do not require data to be downloaded.
    """
    def test_read(self):
        """Test reading from hdf5 file"""
        detector = subspace.Detector()
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'test_data', 'Test_detector.h5')
        detector.read(path)
        self.assertEqual(detector.name, 'Tester')
        self.assertEqual(detector.multiplex, False)
        self.assertEqual(detector.lowcut, 2)
        self.assertEqual(detector.highcut, 9)
        self.assertEqual(detector.filt_order, 4)
        self.assertEqual(detector.dimension, 9)
        self.assertEqual(detector.sampling_rate, 20)

    def test_read_func(self):
        """Check that the read function works too."""
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'test_data', 'Test_detector.h5')
        detector = subspace.read_detector(path)
        _detector = subspace.Detector()
        _detector.read(path)
        self.assertEqual(detector, _detector)

    def test_energy_capture(self):
        """Check that the energy capture calc works okay"""
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'test_data', 'Test_detector.h5')
        detector = subspace.read_detector(path)
        energy = detector.energy_capture()
        self.assertTrue(0 < energy < 100)
        self.assertEqual(round(energy), 60)

    def test_partition_fail(self):
        """Check that partition fails with the expected error."""
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'test_data', 'Test_detector.h5')
        detector = subspace.read_detector(path)
        with self.assertRaises(IndexError):
            detector.partition(dimension=40)

    def test_align(self):
        """Check that alignment does as expected."""
        test_stream = Stream(read()[0])
        # Shift it
        length = 15
        st1 = test_stream.copy().trim(test_stream[0].stats.starttime + 3,
                                      test_stream[0].stats.starttime +
                                      3 + length)
        st2 = test_stream.trim(test_stream[0].stats.starttime,
                               test_stream[0].stats.starttime + length)
        aligned = subspace.align_design(design_set=[st1.copy(), st2.copy()],
                                        shift_len=5, reject=0.3,
                                        multiplex=False, plot=False)
        self.assertEqual(aligned[0][0].stats.starttime,
                         aligned[1][0].stats.starttime)

    def test_stat(self):
        """Test that the statistic calculation is the same regardless of
        system."""
        detector = subspace.Detector()
        detector.read(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   'test_data', 'subspace',
                                   'stat_test_detector.h5'))
        detector.partition(2)
        stream = read(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   'test_data', 'subspace', 'test_trace.ms'))
        st = [stream]
        fft_vars = subspace._do_ffts(detector, st, len(detector.stachans))
        stat = subspace._det_stat_freq(fft_vars[0][0], fft_vars[1][0],
                                       fft_vars[2][0], fft_vars[3],
                                       len(detector.stachans), fft_vars[4],
                                       fft_vars[5])
        self.assertEqual((stat.max().round(6) - 0.229755).round(6), 0)


@pytest.mark.superslow
@pytest.mark.network
class SubspaceTestingMethods(unittest.TestCase):
    """
    Main tests for the subspace module.
    """
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        """Set up the test templates."""
        cls.templates, cls.st = get_test_data()

    def test_write(self):
        """Test writing to an hdf5 file"""
        templates = copy.deepcopy(self.templates)
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=True, shift_len=0.8,
                           reject=0.2)
        detector.write('Test_file.h5')
        self.assertTrue(os.path.isfile('Test_file.h5'))
        os.remove('Test_file.h5')
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=0.8,
                           reject=0.2)
        detector.write('Test_file.h5')
        self.assertTrue(os.path.isfile('Test_file.h5'))
        os.remove('Test_file.h5')

    def test_create_multiplexed_unaligned(self):
        """Test subspace creation - checks that np.dot(U.T, U) is identity."""
        templates = copy.deepcopy(self.templates)
        templates = [template.select(station='TMWZ') for template in templates]
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=False, shift_len=0)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(
                identity, np.diag(np.ones(len(identity), dtype=np.float16))))
        comparison_detector = subspace.read_detector(
                os.path.join(os.path.abspath(
                    os.path.dirname(__file__)), 'test_data', 'subspace',
                    'master_detector_multi_unaligned.h5'))
        # Run to re-fresh file after SVD changes upstream
        # detector.write(os.path.join(os.path.abspath(
        #             os.path.dirname(__file__)), 'test_data', 'subspace',
        #             'master_detector_multi_unaligned.h5'))
        for key in ['name', 'sampling_rate', 'multiplex', 'lowcut', 'highcut',
                    'filt_order', 'dimension', 'stachans']:
            # print(key)
            self.assertEqual(comparison_detector.__getattribute__(key),
                             detector.__getattribute__(key))
        for key in ['data', 'u', 'v', 'sigma']:
            list_item = detector.__getattribute__(key)
            other_list = comparison_detector.__getattribute__(key)
            self.assertEqual(len(list_item), len(other_list))
            for item, other_item in zip(list_item, other_list):
                print(f"{key} is not equal")
                if not np.allclose(np.abs(item), np.abs(other_item)):
                    print(item)
                    print(other_item)
                    print("Differences:")
                    print(item - other_item)
                    print(f"Max difference: {np.max(np.abs(item - other_item))}")
                self.assertTrue(np.allclose(np.abs(item), np.abs(other_item),
                                            atol=0.001))
        # Finally check that the __eq__ method works if all the above passes.
        self.assertEqual(detector, comparison_detector)

    def test_create_nonmultiplexed_unaligned(self):
        """Test creation of a non-multiplexed detector."""
        # Test a non-multiplexed version
        detector = subspace.Detector()
        templates = copy.deepcopy(self.templates)
        templates = [template.select(station='TMWZ') for template in templates]
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=False, shift_len=0)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(
                identity, np.diag(np.ones(len(identity), dtype=np.float16))))
        comparison_detector = subspace.read_detector(
            os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                'test_data', 'subspace', 'master_detector_unaligned.h5'))
        # Run to re-fresh file after SVD changes upstream
        # detector.write(os.path.join(
        #         os.path.abspath(os.path.dirname(__file__)),
        #         'test_data', 'subspace', 'master_detector_unaligned.h5'))
        for key in ['name', 'sampling_rate', 'multiplex', 'lowcut', 'highcut',
                    'filt_order', 'dimension', 'stachans']:
            # print(key)
            self.assertEqual(comparison_detector.__getattribute__(key),
                             detector.__getattribute__(key))
        for key in ['sigma', 'v', 'u', 'data']:
            # print(key)
            list_item = detector.__getattribute__(key)
            other_list = comparison_detector.__getattribute__(key)
            self.assertEqual(len(list_item), len(other_list))
            for item, other_item in zip(list_item, other_list):
                if not np.allclose(np.abs(item), np.abs(other_item)):
                    print(f"Well fuck. {key} is different...")
                    print(item)
                    print(other_item)
                    print("Differences:")
                    print(item - other_item)
                    print(f"Max difference: {np.max(np.abs(item - other_item))}")
                self.assertTrue(np.allclose(np.abs(item), np.abs(other_item),
                                            atol=0.001))
        # Finally check that the __eq__ method works if all the above passes.
        self.assertEqual(detector, comparison_detector)

    def test_create_multiplexed_aligned(self):
        """Test subspace creation - checks that np.dot(U.T, U) is identity."""
        templates = copy.deepcopy(self.templates)
        templates = [template.select(station='TMWZ') for template in templates]
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(
            streams=templates, lowcut=2, highcut=9, filt_order=4,
            sampling_rate=20, multiplex=True, name=str('Tester'), align=True,
            shift_len=3.0, reject=0.2)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(
                identity, np.diag(np.ones(len(identity), dtype=np.float16))))
        comparison_detector = subspace.read_detector(
            os.path.join(os.path.abspath(
                os.path.dirname(__file__)), 'test_data', 'subspace',
                'master_detector_multi.h5'))
        # Run to re-fresh file after SVD changes upstream
        # detector.write(os.path.join(os.path.abspath(
        #         os.path.dirname(__file__)), 'test_data', 'subspace',
        #         'master_detector_multi.h5'))
        for key in ['name', 'sampling_rate', 'multiplex', 'lowcut', 'highcut',
                    'filt_order', 'dimension', 'stachans']:
            # print(key)
            self.assertEqual(comparison_detector.__getattribute__(key),
                             detector.__getattribute__(key))
        for key in ['data', 'u', 'v', 'sigma']:
            # print(key)
            list_item = detector.__getattribute__(key)
            other_list = comparison_detector.__getattribute__(key)
            self.assertEqual(len(list_item), len(other_list))
            for item, other_item in zip(list_item, other_list):
                if not np.allclose(np.abs(item), np.abs(other_item)):
                    print(item)
                    print(other_item)
                self.assertTrue(np.allclose(np.abs(item), np.abs(other_item),
                                            atol=0.001))
        # Finally check that the __eq__ method works if all the above passes.
        self.assertEqual(detector, comparison_detector)

    def test_create_nonmultiplexed_aligned(self):
        """Test creation of a non-multiplexed detector."""
        # Test a non-multiplexed version
        detector = subspace.Detector()
        templates = copy.deepcopy(self.templates)
        templates = [template.select(station='TMWZ') for template in templates]
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=6,
                           reject=0.2)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(
                identity, np.diag(np.ones(len(identity), dtype=np.float16))))
        comparison_detector = subspace.read_detector(
            os.path.join(os.path.abspath(os.path.dirname(__file__)),
                         'test_data', 'subspace', 'master_detector.h5'))
        # Run to re-fresh file after SVD changes upstream
        # detector.write(os.path.join(os.path.abspath(os.path.dirname(__file__)),
        #                'test_data', 'subspace', 'master_detector.h5'))
        for key in ['name', 'sampling_rate', 'multiplex', 'lowcut', 'highcut',
                    'filt_order', 'dimension', 'stachans']:
            # print(key)
            self.assertEqual(comparison_detector.__getattribute__(key),
                             detector.__getattribute__(key))
        for key in ['data', 'u', 'v', 'sigma']:
            # print(key)
            list_item = detector.__getattribute__(key)
            other_list = comparison_detector.__getattribute__(key)
            self.assertEqual(len(list_item), len(other_list))
            for item, other_item in zip(list_item, other_list):
                self.assertEqual(item.shape, other_item.shape)
                if not np.allclose(np.abs(item),
                                   np.abs(other_item)):
                    print(key)
                    print(item)
                    print(other_item)
                self.assertTrue(
                    np.allclose(np.abs(item), np.abs(other_item), atol=0.005))

    def test_refactor(self):
        """
        Test subspace refactoring, checks that np.dot(U.T, U) is identity.
        """
        templates = copy.deepcopy(self.templates)
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=False, shift_len=None)
        for dim in range(2, len(detector.v[0])):
            detector.partition(dim)
            for u in detector.data:
                identity = np.dot(u.T, u).astype(np.float16)
                self.assertTrue(np.allclose(
                    identity, np.diag(np.ones(len(identity),
                                              dtype=np.float16))))
        # Test a non-multiplexed version
        detector = subspace.Detector()
        templates = copy.deepcopy(self.templates)
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=0.2,
                           reject=0.0)
        for dim in range(2, len(detector.v[0])):
            detector.partition(dim)
            for u in detector.data:
                identity = np.dot(u.T, u).astype(np.float16)
                self.assertTrue(np.allclose(
                    identity, np.diag(np.ones(len(identity),
                                              dtype=np.float16))))

    def test_detect(self):
        """Test standard detection with known result."""

        templates = copy.deepcopy(self.templates)
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=True,
                           shift_len=4, reject=0.3,
                           no_missed=False).partition(4)
        st = self.st
        detections = detector.detect(st=st, threshold=0.2, trig_int=4)
        self.assertEqual(len(detections), 34)

    def test_not_multiplexed(self):
        """Test that a non-multiplexed detector gets the same result."""
        templates = copy.deepcopy(self.templates)
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True,
                           shift_len=4, reject=0.3).partition(4)
        st = self.st
        detections = detector.detect(st=st, threshold=0.5, trig_int=4,
                                     moveout=2, min_trig=5)
        self.assertEqual(len(detections), 16)

    def test_multi_detectors(self):
        """Test the efficient looping in subspace."""
        templates = copy.deepcopy(self.templates)
        detector1 = subspace.Detector()
        detector1.construct(streams=templates, lowcut=2, highcut=9,
                            filt_order=4, sampling_rate=20, multiplex=False,
                            name=str('Tester1'), align=True,
                            shift_len=6, reject=0.2).partition(4)
        templates = copy.deepcopy(self.templates)
        detector2 = subspace.Detector()
        detector2.construct(streams=templates[0:20], lowcut=2, highcut=9,
                            filt_order=4, sampling_rate=20, multiplex=False,
                            name=str('Tester2'), align=True,
                            shift_len=6, reject=0.2).partition(4)
        detections = subspace.subspace_detect(detectors=[detector1, detector2],
                                              stream=self.st.copy(),
                                              threshold=0.7,
                                              trig_int=10, moveout=5,
                                              min_trig=5,
                                              parallel=False, num_cores=2)
        self.assertEqual(len(detections), 6)
        detections = subspace.subspace_detect(detectors=[detector1, detector2],
                                              stream=self.st.copy(),
                                              threshold=0.7,
                                              trig_int=10, moveout=5,
                                              min_trig=5,
                                              parallel=True, num_cores=2)
        self.assertEqual(len(detections), 6)

    def partition_fail(self):
        templates = copy.deepcopy(self.templates)
        detector2 = subspace.Detector()
        with self.assertRaises(IndexError):
            detector2.construct(streams=templates[0:10], lowcut=2, highcut=9,
                                filt_order=4, sampling_rate=20,
                                multiplex=False, name=str('Tester'),
                                align=True, shift_len=6,
                                reject=0.2).partition(9)


def get_test_data():
    """
    Generate a set of waveforms from GeoNet for use in subspace testing

    :return: List of cut templates with no filters applied
    :rtype: list
    """
    from http.client import IncompleteRead
    from obspy import UTCDateTime
    from eqcorrscan.utils.catalog_utils import filter_picks
    from eqcorrscan.utils.clustering import catalog_cluster
    from obspy.clients.fdsn import Client

    client = Client("GEONET")
    cat = client.get_events(
        minlatitude=-40.98, maxlatitude=-40.85, minlongitude=175.4,
        maxlongitude=175.5, starttime=UTCDateTime(2016, 5, 11),
        endtime=UTCDateTime(2016, 5, 13))
    cat = filter_picks(catalog=cat, top_n_picks=5)
    stachans = list(set([
        (pick.waveform_id.station_code, pick.waveform_id.channel_code)
        for event in cat for pick in event.picks]))
    clusters = catalog_cluster(
        catalog=cat, thresh=2, show=False, metric='distance')
    cluster = sorted(clusters, key=lambda c: len(c))[-1]
    client = Client('GEONET')
    design_set = []
    st = Stream()
    for event in cluster:
        # This print is just in to force some output during long running test
        print("Downloading for event {0}".format(event.resource_id))
        bulk_info = []
        t1 = event.origins[0].time + 5
        t2 = t1 + 15.1
        for station, channel in stachans:
            bulk_info.append(
                ('NZ', station, '10', channel[0:2] + '?', t1, t2))
        st += client.get_waveforms_bulk(bulk=bulk_info)
    for event in cluster:
        t1 = event.origins[0].time + 5
        t2 = t1 + 15
        design_set.append(st.copy().trim(t1, t2))
    t1 = UTCDateTime(2016, 5, 11, 19)
    t2 = UTCDateTime(2016, 5, 11, 20)
    bulk_info = [('NZ', stachan[0], '10', stachan[1][0:2] + '?', t1, t2)
                 for stachan in stachans]
    st = Stream()
    for _bulk in bulk_info:
        try:
            st += client.get_waveforms(*_bulk)
        except IncompleteRead:
            print(f"Could not download {_bulk}")
    st.merge().detrend('simple').trim(starttime=t1, endtime=t2)
    return design_set, st


if __name__ == '__main__':
    unittest.main()
