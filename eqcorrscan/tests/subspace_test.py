"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


import numpy as np
import unittest
import os
import copy

from obspy import Stream, read

from eqcorrscan.core import subspace, subspace_statistic


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
        stream = read(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   'test_data', 'subspace', 'test_trace.ms'))
        tr_data = stream[0].data
        stat = subspace_statistic.det_statistic(detector.data[0].
                                                astype(np.float32),
                                                tr_data.astype(np.float32))
        self.assertEqual((stat.max().round(6) - 0.252336).round(6), 0)


class SubspaceTestingMethods(unittest.TestCase):
    """
    Main tests for the subspace module.
    """
    @classmethod
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
            self.assertTrue(np.allclose(identity,
                                        np.diag(np.ones(len(identity),
                                                        dtype=np.float16))))
        comparison_detector = \
            subspace.read_detector(os.path.
                                   join(os.path.
                                        abspath(os.path.
                                                dirname(__file__)),
                                        'test_data', 'subspace',
                                        'master_detector_multi_unaligned.h5'))
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
                if not np.allclose(item, other_item):
                    print(item)
                    print(other_item)
                self.assertTrue(np.allclose(item, other_item))
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
            self.assertTrue(np.allclose(identity,
                                        np.diag(np.ones(len(identity),
                                                        dtype=np.float16))))
        comparison_detector = \
            subspace.read_detector(os.path.
                                   join(os.path.
                                        abspath(os.path.
                                                dirname(__file__)),
                                        'test_data', 'subspace',
                                        'master_detector_unaligned.h5'))
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
                if not np.allclose(item, other_item):
                    print(item)
                    print(other_item)
                self.assertTrue(np.allclose(item, other_item))
        # Finally check that the __eq__ method works if all the above passes.
        self.assertEqual(detector, comparison_detector)

    def test_create_multiplexed_aligned(self):
        """Test subspace creation - checks that np.dot(U.T, U) is identity."""
        templates = copy.deepcopy(self.templates)
        templates = [template.select(station='TMWZ') for template in templates]
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=True, shift_len=3.0,
                           reject=0.2)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(identity,
                                        np.diag(np.ones(len(identity),
                                                        dtype=np.float16))))
        comparison_detector = \
            subspace.read_detector(os.path.
                                   join(os.path.
                                        abspath(os.path.
                                                dirname(__file__)),
                                        'test_data', 'subspace',
                                        'master_detector_multi.h5'))
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
                if not np.allclose(item, other_item):
                    print(item)
                    print(other_item)
                self.assertTrue(np.allclose(item, other_item))
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
            self.assertTrue(np.allclose(identity,
                                        np.diag(np.ones(len(identity),
                                                        dtype=np.float16))))
        comparison_detector = \
            subspace.read_detector(os.path.join(os.path.
                                                abspath(os.path.
                                                        dirname(__file__)),
                                                'test_data', 'subspace',
                                                'master_detector.h5'))
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
            # for item, other_item in zip(list_item, other_list):
            #     print(item.shape)
            #     print(other_item.shape)
            #     print('Next')
            for item, other_item in zip(list_item, other_list):
                self.assertEqual(item.shape, other_item.shape)
                if not np.allclose(item, other_item):
                    print(item)
                    print(other_item)
                self.assertTrue(np.allclose(item, other_item))
        # Finally check that the __eq__ method works if all the above passes.
        self.assertEqual(detector, comparison_detector)

    def test_refactor(self):
        """Test subspace refactoring, checks that np.dot(U.T, U) is\
         identity."""
        templates = copy.deepcopy(self.templates)
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=False, shift_len=None)
        for dim in range(2, len(detector.u[0])):
            detector.partition(dim)
            for u in detector.data:
                identity = np.dot(u.T, u).astype(np.float16)
                self.assertTrue(np.allclose(identity,
                                            np.diag(np.
                                                    ones(len(identity),
                                                         dtype=np.float16))))
        # Test a non-multiplexed version
        detector = subspace.Detector()
        templates = copy.deepcopy(self.templates)
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=0.2,
                           reject=0.0)
        for dim in range(2, len(detector.u[0])):
            detector.partition(dim)
            for u in detector.data:
                identity = np.dot(u.T, u).astype(np.float16)
                self.assertTrue(np.allclose(identity,
                                            np.diag(np.
                                                    ones(len(identity),
                                                         dtype=np.float16))))

    def test_detect(self):
        """Test standard detection with known result."""

        templates = copy.deepcopy(self.templates)
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=True,
                           shift_len=6, reject=0.2).partition(9)
        st = self.st
        detections = detector.detect(st=st, threshold=0.009, trig_int=2,
                                     debug=1)
        self.assertEqual(len(detections), 2)

    def test_not_multiplexed(self):
        """Test that a non-multiplexed detector gets the same result."""
        templates = copy.deepcopy(self.templates)
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True,
                           shift_len=6, reject=0.2).partition(9)
        st = self.st
        detections = detector.detect(st=st, threshold=0.05, trig_int=4,
                                     debug=0, moveout=2, min_trig=5)
        self.assertEqual(len(detections), 1)

    def test_multi_detectors(self):
        """Test the efficient looping in subspace."""
        templates = copy.deepcopy(self.templates)
        detector1 = subspace.Detector()
        detector1.construct(streams=templates, lowcut=2, highcut=9,
                            filt_order=4, sampling_rate=20, multiplex=False,
                            name=str('Tester1'), align=True,
                            shift_len=6, reject=0.2).partition(9)
        templates = copy.deepcopy(self.templates)
        detector2 = subspace.Detector()
        detector2.construct(streams=templates[0:20], lowcut=2, highcut=9,
                            filt_order=4, sampling_rate=20, multiplex=False,
                            name=str('Tester2'), align=True,
                            shift_len=6, reject=0.2).partition(9)
        detections = subspace.subspace_detect(detectors=[detector1, detector2],
                                              stream=self.st.copy(),
                                              threshold=0.05,
                                              trig_int=10, moveout=5,
                                              min_trig=5,
                                              parallel=False, num_cores=2)
        self.assertEqual(len(detections), 4)
        detections = subspace.subspace_detect(detectors=[detector1, detector2],
                                              stream=self.st.copy(),
                                              threshold=0.05,
                                              trig_int=10, moveout=5,
                                              min_trig=5,
                                              parallel=True, num_cores=2)
        self.assertEqual(len(detections), 4)

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
    from eqcorrscan.tutorials.get_geonet_events import get_geonet_events
    from obspy import UTCDateTime
    from eqcorrscan.utils.catalog_utils import filter_picks
    from eqcorrscan.utils.clustering import space_cluster
    from obspy.clients.fdsn import Client

    cat = get_geonet_events(minlat=-40.98, maxlat=-40.85, minlon=175.4,
                            maxlon=175.5, startdate=UTCDateTime(2016, 5, 11),
                            enddate=UTCDateTime(2016, 5, 13))
    cat = filter_picks(catalog=cat, top_n_picks=5)
    stachans = list(set([(pick.waveform_id.station_code,
                          pick.waveform_id.channel_code) for event in cat
                         for pick in event.picks]))
    clusters = space_cluster(catalog=cat, d_thresh=2, show=False)
    cluster = sorted(clusters, key=lambda c: len(c))[-1]
    client = Client('GEONET')
    design_set = []
    bulk_info = []
    for event in cluster:
        t1 = event.origins[0].time + 5
        t2 = t1 + 15
        for station, channel in stachans:
            bulk_info.append(('NZ', station, '*', channel[0:2] + '?', t1, t2))
    st = client.get_waveforms_bulk(bulk=bulk_info)
    for event in cluster:
        t1 = event.origins[0].time + 5
        t2 = t1 + 15
        design_set.append(st.copy().trim(t1, t2))
    t1 = UTCDateTime(2016, 5, 11, 19)
    t2 = UTCDateTime(2016, 5, 11, 20)
    bulk_info = [('NZ', stachan[0], '*',
                  stachan[1][0:2] + '?',
                  t1, t2) for stachan in stachans]
    st = client.get_waveforms_bulk(bulk_info)
    st.merge().detrend('simple').trim(starttime=t1, endtime=t2)
    return design_set, st

if __name__ == '__main__':
    unittest.main()
