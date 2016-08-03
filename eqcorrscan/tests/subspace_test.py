"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.core import subspace, subspace_statistic
from eqcorrscan.core.subspace import _subspace_process
import numpy as np
import unittest
from obspy import Stream
import obspy
if int(obspy.__version__.split('.')[0]) >= 1:
    from obspy.clients.fdsn import Client
else:
    from obspy.fdsn import Client
from obspy import read
import os
import copy


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

    # def test_synthetic(self):
    #     """Test a synthetic case."""
    #     self.assertEqual('This passes', 'Not yet')

    def test_write(self):
        """Test writing to an hdf5 file"""
        templates = copy.deepcopy(self.templates)
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=True, shift_len=0.2)
        detector.write('Test_file.h5')
        self.assertTrue(os.path.isfile('Test_file.h5'))
        os.remove('Test_file.h5')
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=0.2)
        detector.write('Test_file.h5')
        self.assertTrue(os.path.isfile('Test_file.h5'))
        os.remove('Test_file.h5')

    def test_create(self):
        """Test subspace creation - checks that np.dot(U.T, U) is identity."""
        templates = copy.deepcopy(self.templates)
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=True, shift_len=0.2)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(identity,
                                        np.diag(np.ones(len(identity),
                                                        dtype=np.float16))))
        # Test a non-multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=0.2)
        for u in detector.data:
            identity = np.dot(u.T, u).astype(np.float16)
            self.assertTrue(np.allclose(identity,
                                        np.diag(np.ones(len(identity),
                                                        dtype=np.float16))))

    def test_refactor(self):
        """Test subspace refactoring, checks that np.dot(U.T, U) is\
         identity."""
        templates = copy.deepcopy(self.templates)
        # Test a multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=True,
                           name=str('Tester'), align=False, shift_len=None)
        for dim in range(2, len(detector.u)):
            detector.partition(dim)
            for u in detector.data:
                identity = np.dot(u.T, u).astype(np.float16)
                self.assertTrue(np.allclose(identity,
                                            np.diag(np.ones(len(identity),
                                                            dtype=np.float16))))
        # Test a non-multiplexed version
        detector = subspace.Detector()
        detector.construct(streams=templates, lowcut=2, highcut=9,
                           filt_order=4, sampling_rate=20, multiplex=False,
                           name=str('Tester'), align=True, shift_len=0.2)
        for dim in range(2, len(detector.u)):
            detector.partition(dim)
            for u in detector.data:
                identity = np.dot(u.T, u).astype(np.float16)
                self.assertTrue(np.allclose(identity,
                                            np.diag(np.ones(len(identity),
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
        self.assertEqual(len(detections), 2)


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
                            maxlon=175.5, startdate=UTCDateTime(2016, 5, 1),
                            enddate=UTCDateTime(2016, 5, 20))
    cat = filter_picks(catalog=cat, top_n_picks=5)
    stachans = list(set([(pick.waveform_id.station_code,
                          pick.waveform_id.channel_code) for event in cat
                         for pick in event.picks]))
    clusters = space_cluster(catalog=cat, d_thresh=2, show=False)
    cluster = sorted(clusters, key=lambda c: len(c))[-1]
    client = Client('GEONET')
    design_set = []
    for event in cluster:
        t1 = event.origins[0].time
        t2 = t1 + 25
        bulk_info = []
        for station, channel in stachans:
            bulk_info.append(('NZ', station, '*', channel[0:2] + '?', t1, t2))
        st = client.get_waveforms_bulk(bulk=bulk_info)
        st.trim(t1, t2)
        design_set.append(st)
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