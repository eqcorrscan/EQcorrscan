"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.core import subspace
import numpy as np
import unittest
from obspy import UTCDateTime
import obspy
if int(obspy.__version__.split('.')[0]) >= 1:
    from obspy.clients.fdsn import Client
else:
    from obspy.fdsn import Client
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
        self.assertEqual(detector.name, 'test')
        self.assertEqual(detector.multiplex, True)
        self.assertEqual(detector.lowcut, 2)
        self.assertEqual(detector.highcut, 9)
        self.assertEqual(detector.filt_order, 3)
        self.assertEqual(detector.dimension, 3)
        self.assertEqual(detector.sampling_rate, 20)

    def test_read_func(self):
        """Check that the read function works too."""
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'test_data', 'Test_detector.h5')
        detector = subspace.read_detector(path)
        _detector = subspace.Detector()
        _detector.read(path)
        self.assertEqual(detector, _detector)

class SubspaceTestingMethods(unittest.TestCase):
    """
    Main tests for the subspace module.
    """
    @classmethod
    def setUpClass(cls):
        """Set up the test templates."""
        cls.templates = get_test_data()

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
                           shift_len=0.2).partition(4)
        t1 = UTCDateTime(2016, 5, 11, 19)
        t2 = UTCDateTime(2016, 5, 12)
        bulk_info = [('NZ', stachan[0], '*',
                      stachan[1][0] + '?' + stachan[1][-1],
                      t1, t2) for stachan in detector.stachans]
        client = Client('GEONET')
        st = client.get_waveforms_bulk(bulk_info)
        st.merge().detrend('simple').trim(starttime=t1, endtime=t2)
        for tr in st:
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        detections = detector.detect(st=st, threshold=0.53, trig_int=2, debug=1)
        self.assertEqual(len(detections), 10)


def get_test_data():
    """
    Generate a set of waveforms from GeoNet for use in subspace testing

    :return: List of cut templates with no filters applied
    :rtype: list
    """
    from eqcorrscan.tutorials.get_geonet_events import get_geonet_events
    from obspy import UTCDateTime, Catalog
    from eqcorrscan.utils.catalog_utils import filter_picks
    from eqcorrscan.utils.clustering import space_cluster
    from eqcorrscan.core import template_gen

    cat = get_geonet_events(minlat=-40.98, maxlat=-40.85, minlon=175.4,
                            maxlon=175.5, startdate=UTCDateTime(2016, 5, 1),
                            enddate=UTCDateTime(2016, 5, 20))
    cat = filter_picks(catalog=cat, top_n_picks=5)
    # Then remove events with fewer than three picks
    cat = Catalog([event for event in cat if len(event.picks) >= 3])
    # In this tutorial we will only work on one cluster, defined spatially.
    # You can work on multiple clusters, or try to whole set.
    clusters = space_cluster(catalog=cat, d_thresh=2, show=False)
    # We will work on the largest cluster
    cluster = sorted(clusters, key=lambda c: len(c))[-1]
    # This cluster contains 42 events, we will now generate simple waveform
    # templates for each of them
    templates = template_gen.from_client(catalog=cluster,
                                         client_id='GEONET',
                                         lowcut=None, highcut=None,
                                         samp_rate=100.0, filt_order=4,
                                         length=2.0, prepick=0.5,
                                         swin='all', process_len=3600,
                                         debug=0, plot=False)
    return templates

if __name__ == '__main__':
    unittest.main()