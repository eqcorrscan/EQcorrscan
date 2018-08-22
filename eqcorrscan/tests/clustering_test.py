"""
A series of test functions for the utils.clustering module in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import pytest
import numpy as np
import os
import glob
import warnings

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy import read
from obspy.clients.fdsn.header import FDSNException

from eqcorrscan.tutorials.template_creation import mktemplates
from eqcorrscan.utils.mag_calc import dist_calc
from eqcorrscan.utils.clustering import (
    cross_chan_coherence, distance_matrix, cluster, group_delays, svd, SVD,
    empirical_svd, empirical_SVD, SVD_2_stream, svd_to_stream, corr_cluster,
    dist_mat_km, catalog_cluster, space_time_cluster, remove_unclustered,
    distance_link)


@pytest.mark.network
class ClusteringTests(unittest.TestCase):
    """Main testing routines"""
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        cls.st1 = read()
        cls.st2 = cls.st1.copy()
        cls.testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        warnings.simplefilter("ignore")
        client = Client("https://earthquake.usgs.gov")
        starttime = UTCDateTime("2002-01-01")
        endtime = UTCDateTime("2002-01-02")
        cls.cat = client.get_events(starttime=starttime, endtime=endtime,
                                    minmagnitude=6)

    def test_cross_chan_coherence(self):
        """Initial test to ensure cross_chan_coherence runs."""
        cccoh, i = cross_chan_coherence(st1=self.st1.copy(),
                                        st2=self.st2.copy())
        self.assertEqual(cccoh, 1)
        self.assertEqual(i, 0)

    def test_cross_chan_coherence_shifted(self):
        """Initial test to ensure cross_chan_coherence runs."""
        cccoh, i = cross_chan_coherence(st1=self.st1.copy(),
                                        st2=self.st2.copy(), allow_shift=True)
        self.assertEqual(round(cccoh, 6), 1)
        self.assertEqual(i, 0)

    def test_inverted_coherence(self):
        """Reverse channels and ensure we get -1"""
        st2 = self.st2.copy()
        for tr in st2:
            tr.data *= -1
        cccoh, i = cross_chan_coherence(st1=self.st1.copy(), st2=st2, i=1)
        self.assertEqual(cccoh, -1)
        self.assertEqual(i, 1)

    def test_known_coherence(self):
        """Test for a real stream case"""
        st1 = read(os.path.join(self.testing_path, 'WAV', 'TEST_',
                                '2013-09-01-0410-35.DFDPC_024_00'))
        st2 = read(os.path.join(self.testing_path, 'WAV', 'TEST_',
                                '2013-09-01-2040-11.DFDPC_039_00'))
        cccoh, i = cross_chan_coherence(st1=st1, st2=st2)
        self.assertTrue(cccoh < 0.01)
        self.assertEqual(i, 0)

    def test_distance_matrix(self):
        """Test that we can create a useful distance matrix."""
        testing_path = os.path.join(self.testing_path, 'WAV', 'TEST_')
        stream_files = glob.glob(os.path.join(testing_path, '*DFDPC*'))[0:10]
        stream_list = [read(stream_file) for stream_file in stream_files]
        dist_mat = distance_matrix(stream_list=stream_list, cores=4)
        self.assertEqual(dist_mat.shape[0], len(stream_list))
        self.assertEqual(dist_mat.shape[1], len(stream_list))

    def test_unclustered(self):
        """Test clustering on unclustered data..."""
        testing_path = os.path.join(self.testing_path, 'WAV', 'TEST_')
        stream_files = glob.glob(os.path.join(testing_path, '*DFDPC*'))[0:10]
        stream_list = [(read(stream_file), i)
                       for i, stream_file in enumerate(stream_files)]
        groups = cluster(template_list=stream_list, show=False,
                         corr_thresh=0.3)
        self.assertEqual(len(groups), 10)  # They shouldn't cluster at all
        # Test setting a number of cores
        groups_2 = cluster(template_list=stream_list, show=False,
                           corr_thresh=0.3, cores=2, debug=2,
                           save_corrmat=True)
        self.assertTrue(os.path.isfile('dist_mat.npy'))
        os.remove('dist_mat.npy')
        self.assertEqual(len(groups_2), 10)  # They shouldn't cluster at all
        self.assertEqual(groups, groups_2)

    def test_clustered(self):
        """Test clustering on clustered data..."""
        testing_path = os.path.join(self.testing_path, 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [(read(stream_file), i)
                       for i, stream_file in enumerate(stream_files)]
        for stream in stream_list:
            for tr in stream[0]:
                if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                    stream[0].remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        groups = cluster(template_list=stream_list, show=False,
                         corr_thresh=0.3)
        self.assertEqual(len(groups), 9)  # They should cluster reasonably

    def test_corr_cluster(self):
        """Test the corr_cluster function."""
        samp_rate = 40
        testing_path = os.path.join(self.testing_path, 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        for stream in stream_list:
            for tr in stream:
                if not (tr.stats.station == 'GCSZ' and
                        tr.stats.channel == 'EHZ'):
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=2.0, freqmax=10.0)
                tr.resample(sampling_rate=samp_rate)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        trace_list = [stream[0] for stream in stream_list]
        output = corr_cluster(trace_list=trace_list, thresh=0.7)
        self.assertFalse(output.all())


class DistanceClusterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        client = Client("https://earthquake.usgs.gov")
        starttime = UTCDateTime("2002-01-01")
        endtime = UTCDateTime("2002-01-02")
        cls.cat = client.get_events(starttime=starttime, endtime=endtime,
                                    minmagnitude=4)
        cls.distance_threshold = 1000
        cls.time_threshold = 7200

    def test_distance_linkage(self):
        distance_cutoff = 100
        linkage = distance_link(self.cat, distance_cutoff, num_threads=1)
        for i, event in enumerate(self.cat):
            master_ori = event.preferred_origin() or event.origins[0]
            for j, other_event in enumerate(self.cat):
                slave_ori = (other_event.preferred_origin() or
                             other_event.origins[0])
                distance = dist_calc(
                    (master_ori.latitude, master_ori.longitude,
                     master_ori.depth / 1000),
                    (slave_ori.latitude, slave_ori.longitude,
                     slave_ori.depth / 1000))
                if distance < distance_cutoff:
                    assert linkage[i][j]
                else:
                    assert not linkage[i][j]

    def test_distance_linkage_parallel(self):
        distance_cutoff = 100
        linkage = distance_link(self.cat, distance_cutoff, num_threads=4)
        for i, event in enumerate(self.cat):
            master_ori = event.preferred_origin() or event.origins[0]
            for j, other_event in enumerate(self.cat):
                slave_ori = (other_event.preferred_origin() or
                             other_event.origins[0])
                distance = dist_calc(
                    (master_ori.latitude, master_ori.longitude,
                     master_ori.depth / 1000),
                    (slave_ori.latitude, slave_ori.longitude,
                     slave_ori.depth / 1000))
                if distance < distance_cutoff:
                    assert linkage[i][j]
                else:
                    assert not linkage[i][j]

    def test_remove_unclustered(self):
        distance_cutoff = 100
        cat_back = remove_unclustered(self.cat.copy(), distance_cutoff,
                                      num_threads=1)
        assert len(cat_back) < len(self.cat)
        # Check that all events have at least one event within distance_cutoff
        for i, event in enumerate(cat_back):
            master_ori = event.preferred_origin() or event.origins[0]
            distances = []
            for j, other_event in enumerate(cat_back):
                slave_ori = (other_event.preferred_origin() or
                             other_event.origins[0])
                if i == j:
                    continue
                distances.append(dist_calc(
                    (master_ori.latitude, master_ori.longitude,
                     master_ori.depth / 1000),
                    (slave_ori.latitude, slave_ori.longitude,
                     slave_ori.depth / 1000)))
            assert any(np.array(distances) < distance_cutoff)
        # Check that all events not in the catalog are correctly ignored
        for i, event in enumerate(self.cat):
            if event not in cat_back:
                master_ori = event.preferred_origin() or event.origins[0]
                distances = []
                for j, other_event in enumerate(self.cat):
                    slave_ori = (other_event.preferred_origin() or
                                 other_event.origins[0])
                    if i == j:
                        continue
                    distances.append(dist_calc(
                        (master_ori.latitude, master_ori.longitude,
                         master_ori.depth / 1000),
                        (slave_ori.latitude, slave_ori.longitude,
                         slave_ori.depth / 1000)))
                assert all(np.array(distances) > distance_cutoff)

    def test_remove_unclustered_parallel(self):
        distance_cutoff = 100
        cat_back = remove_unclustered(self.cat.copy(), distance_cutoff,
                                      num_threads=4)
        assert len(cat_back) < len(self.cat)
        for i, event in enumerate(cat_back):
            master_ori = event.preferred_origin() or event.origins[0]
            distances = []
            for j, other_event in enumerate(cat_back):
                slave_ori = (other_event.preferred_origin() or
                             other_event.origins[0])
                if i == j:
                    continue
                distances.append(dist_calc(
                    (master_ori.latitude, master_ori.longitude,
                     master_ori.depth / 1000),
                    (slave_ori.latitude, slave_ori.longitude,
                     slave_ori.depth / 1000)))
            assert any(np.array(distances) < distance_cutoff)
        # Check that all events not in the catalog are correctly ignored
        for i, event in enumerate(self.cat):
            if event not in cat_back:
                master_ori = event.preferred_origin() or event.origins[0]
                distances = []
                for j, other_event in enumerate(self.cat):
                    slave_ori = (other_event.preferred_origin() or
                                 other_event.origins[0])
                    if i == j:
                        continue
                    distances.append(dist_calc(
                        (master_ori.latitude, master_ori.longitude,
                         master_ori.depth / 1000),
                        (slave_ori.latitude, slave_ori.longitude,
                         slave_ori.depth / 1000)))
                assert all(np.array(distances) > distance_cutoff)

    def test_dist_mat_km(self):
        """Test spacial clustering."""
        dist_mat = dist_mat_km(self.cat)
        self.assertEqual(len(dist_mat), len(self.cat))
        # Diagonal should be zeros
        self.assertTrue(np.all(dist_mat.diagonal() == 0))
        # Should be symmetric
        for i in range(len(self.cat)):
            for j in range(len(self.cat)):
                self.assertEqual(dist_mat[i, j], dist_mat[j, i])
                master_ori = (
                    self.cat[i].preferred_origin() or self.cat[i].origins[0])
                slave_ori = (
                    self.cat[j].preferred_origin() or self.cat[j].origins[0])
                self.assertAlmostEqual(
                    dist_mat[i, j], dist_calc(
                        (master_ori.latitude, master_ori.longitude,
                         master_ori.depth / 1000),
                        (slave_ori.latitude, slave_ori.longitude,
                         slave_ori.depth / 1000)), 6)

    def test_space_cluster(self):
        """Test the wrapper around dist_mat_km."""
        groups = catalog_cluster(
            catalog=self.cat, thresh=self.distance_threshold,
            metric="distance", show=False)
        self.assertEqual(len([ev for group in groups for ev in group]),
                         len(self.cat))
        # Check that events within each group are within distance-threshold
        for group in groups:
            if len(group) > 1:
                master_ori = group[0].preferred_origin() or group[0].origins[0]
                for event in group[1:]:
                    slave_ori = event.preferred_origin() or event.origins[0]
                    self.assertLessEqual(dist_calc(
                        (master_ori.latitude, master_ori.longitude,
                         master_ori.depth / 1000),
                        (slave_ori.latitude, slave_ori.longitude,
                         slave_ori.depth / 1000)), self.distance_threshold)
        # Check that groups are separated by at least distance-threshold
        for i, group in enumerate(groups):
            for event in group:
                master_ori = event.preferred_origin() or event.origins[0]
                for j, other_group in enumerate(groups):
                    if i == j:
                        continue
                    for other_event in other_group:
                        slave_ori = (
                            other_event.preferred_origin() or
                            other_event.origins[0])
                        self.assertGreater(dist_calc(
                            (master_ori.latitude, master_ori.longitude,
                             master_ori.depth / 1000),
                            (slave_ori.latitude, slave_ori.longitude,
                             slave_ori.depth / 1000)), self.distance_threshold)

    def test_time_cluster(self):
        """Test the wrapper around dist_mat_time."""
        groups = catalog_cluster(
            catalog=self.cat, thresh=self.time_threshold,
            metric="time", show=False)
        self.assertEqual(len([ev for group in groups for ev in group]),
                         len(self.cat))
        # Check that events within each group are within time-threshold
        for group in groups:
            if len(group) > 1:
                master_ori = group[0].preferred_origin() or group[0].origins[0]
                for event in group[1:]:
                    slave_ori = event.preferred_origin() or event.origins[0]
                    self.assertLessEqual(
                        abs(master_ori.time - slave_ori.time),
                        self.time_threshold)
        # Check that groups are separated by at least time-threshold
        for i, group in enumerate(groups):
            master_times = []
            for master in group:
                master_ori = (
                    master.preferred_origin() or master.origins[0])
                master_times.append(master_ori.time)
            master_times.sort()
            master_median_time = master_times[0] + np.median(
                [m - master_times[0] for m in master_times])
            for j, other_group in enumerate(groups):
                if i == j:
                    continue
                slave_times = []
                for slave in other_group:
                    slave_ori = (
                        slave.preferred_origin() or slave.origins[0])
                    slave_times.append(slave_ori.time)
                slave_times.sort()
                slave_median_time = slave_times[0] + np.median(
                    [s - slave_times[0] for s in slave_times])
                self.assertGreater(
                    abs(master_median_time - slave_median_time),
                    self.time_threshold)

    def test_space_time_cluster(self):
        """Test clustering in space and time."""
        groups = space_time_cluster(
            catalog=self.cat, t_thresh=self.time_threshold,
            d_thresh=self.distance_threshold)
        self.assertEqual(len([ev for group in groups for ev in group]),
                         len(self.cat))
        # Check that events within each group are within distance-threshold and
        # time-threshold
        for group in groups:
            if len(group) > 1:
                master_ori = group[0].preferred_origin() or group[0].origins[0]
                for event in group[1:]:
                    slave_ori = event.preferred_origin() or event.origins[0]
                    self.assertLessEqual(dist_calc(
                        (master_ori.latitude, master_ori.longitude,
                         master_ori.depth / 1000),
                        (slave_ori.latitude, slave_ori.longitude,
                         slave_ori.depth / 1000)), self.distance_threshold)
                    self.assertLessEqual(abs(master_ori.time - slave_ori.time),
                                         self.time_threshold)
        # Check that groups are separated by at least distance-threshold,
        # or, by some time.
        for i, group in enumerate(groups):
            master_times = []
            for master in group:
                master_ori = (
                    master.preferred_origin() or master.origins[0])
                master_times.append(master_ori.time)
            master_times.sort()
            master_median_time = master_times[0] + np.median(
                [m - master_times[0] for m in master_times])
            # Just use the origin of one event
            master_ori = group[0].preferred_origin() or group[0].origins[0]
            for j, other_group in enumerate(groups):
                if i == j:
                    continue
                slave_times = []
                for slave in other_group:
                    slave_ori = (
                        slave.preferred_origin() or slave.origins[0])
                    slave_times.append(slave_ori.time)
                slave_times.sort()
                slave_median_time = slave_times[0] + np.median(
                    [s - slave_times[0] for s in slave_times])
                for event in other_group:
                    slave_ori = event.preferred_origin() or event.origins[0]
                    separation = dist_calc(
                        (master_ori.latitude, master_ori.longitude,
                         master_ori.depth / 1000),
                        (slave_ori.latitude, slave_ori.longitude,
                         slave_ori.depth / 1000))
                    if separation < self.distance_threshold:
                        self.assertGreater(
                            abs(master_median_time - slave_median_time),
                            self.time_threshold)
                    else:
                        self.assertGreater(separation, self.distance_threshold)


@pytest.mark.network
class ClusteringTestWarnings(unittest.TestCase):
    """Main testing routines"""
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        cls.st1 = read()
        cls.st2 = cls.st1.copy()
        cls.testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        warnings.simplefilter("always")
        client = Client("https://earthquake.usgs.gov")
        starttime = UTCDateTime("2002-01-01")
        endtime = UTCDateTime("2002-01-02")
        cls.cat = client.get_events(starttime=starttime, endtime=endtime,
                                    minmagnitude=6)

    def test_cross_chan_coherence_non_matching(self):
        """Initial test to ensure cross_chan_coherence runs."""
        st2 = self.st2.copy()
        for tr in st2:
            tr.stats.station += 'A'
        with warnings.catch_warnings(record=True) as w:
            cross_chan_coherence(st1=self.st1.copy(), st2=st2)
            self.assertEqual(len(w), 1)
            self.assertTrue('No matching channels' in str(w[0].message))

    def test_cross_chan_coherence_non_matching_sampling_rates(self):
        """Initial test to ensure cross_chan_coherence runs."""
        st2 = self.st2.copy()
        for tr in st2:
            tr.stats.sampling_rate += 20
        with warnings.catch_warnings(record=True) as w:
            cross_chan_coherence(st1=self.st1.copy(), st2=st2)
            self.assertTrue('Sampling rates do not match' in str(w[0].message))

    def test_delay_grouping(self):
        """Test grouping by delays"""
        testing_path = os.path.join(self.testing_path, 'WAV', 'TEST_')
        stream_files = glob.glob(os.path.join(testing_path, '*DFDPC*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        groups = group_delays(stream_list=stream_list)
        self.assertEqual(len(groups), 1)  # All have same lag-times

        # Make some templates
        try:
            mktemplates(plot=False)
        except FDSNException:
            warnings.warn('FDSN exception raised, is server down?')
            return
        stream_list = []
        for template_no in range(4):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            stream_list.append(template)
        groups = group_delays(stream_list=stream_list)
        self.assertEqual(len(groups), 4)  # None have same lag-times
        # Cleanup the templates
        templates = glob.glob('tutorial_template_?.ms')
        for template in templates:
            os.remove(template)

    def test_svd(self):
        """Test the svd method."""
        testing_path = os.path.join(self.testing_path, 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        for stream in stream_list:
            for tr in stream:
                if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        UVectors, SValues, SVectors, stachans = svd(stream_list=stream_list)
        self.assertEqual(len(SVectors), len(stachans))
        self.assertEqual(len(SValues), len(stachans))
        self.assertEqual(len(UVectors), len(stachans))
        for SVec in SVectors:
            self.assertEqual(len(SVec), len(stream_list))
        with warnings.catch_warnings(record=True) as w:
            SVD(stream_list=stream_list)
            self.assertEqual(len(w), 1)
            self.assertTrue('Depreciated' in str(w[0].message))

    def test_empirical_svd(self):
        """Test the empirical SVD method"""
        testing_path = os.path.join(self.testing_path, 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        for stream in stream_list:
            for tr in stream:
                if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        first_sub, second_sub = empirical_svd(stream_list=stream_list)
        self.assertEqual(len(first_sub), len(second_sub))
        self.assertEqual(len(stream_list[0]), len(first_sub))
        first_sub, second_sub = empirical_svd(stream_list=stream_list,
                                              linear=False)
        self.assertEqual(len(first_sub), len(second_sub))
        self.assertEqual(len(stream_list[0]), len(first_sub))
        with warnings.catch_warnings(record=True) as w:
            empirical_SVD(stream_list=stream_list)
            self.assertEqual(len(w), 1)
            self.assertTrue('Depreciated' in str(w[0].message))

    def test_svd_to_stream(self):
        """Test the conversion of SVD to stream."""
        samp_rate = 100
        testing_path = os.path.join(self.testing_path, 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        for stream in stream_list:
            for tr in stream:
                if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
                tr.resample(sampling_rate=samp_rate)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        SVectors, SValues, Uvectors, stachans = svd(stream_list=stream_list)
        svstreams = svd_to_stream(uvectors=SVectors, stachans=stachans, k=4,
                                  sampling_rate=samp_rate)
        self.assertEqual(len(svstreams), 4)
        with warnings.catch_warnings(record=True) as w:
            SVD_2_stream(uvectors=SVectors, stachans=stachans, k=4,
                         sampling_rate=samp_rate)
            self.assertEqual(len(w), 1)
            self.assertTrue('Depreciated' in str(w[0].message))

    def test_corr_unclustered(self):
        """Test the corr_cluster function."""
        samp_rate = 40
        testing_path = os.path.join(self.testing_path, 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        for stream in stream_list:
            for tr in stream:
                if not (tr.stats.station == 'WHAT2' and
                        tr.stats.channel == 'SH1'):
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=2.0, freqmax=10.0)
                tr.resample(sampling_rate=samp_rate)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        trace_list = [stream[0] for stream in stream_list]
        with warnings.catch_warnings(record=True) as w:
            corr_cluster(trace_list=trace_list, thresh=0.7)
            self.assertEqual(len(w), 1)
            self.assertTrue('Nothing made it past the first 0.6 threshold'
                            in str(w[0].message))


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
