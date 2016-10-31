"""
A series of test functions for the utils.clustering module in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class ClusteringTestMethods(unittest.TestCase):
    """Main testing routines"""
    def test_cross_chan_coherence(self):
        """Initial test to ensure cross_chan_coherence runs."""
        from eqcorrscan.utils.clustering import cross_chan_coherence
        from obspy import read
        st1 = read()
        st2 = read()
        cccoh, i = cross_chan_coherence(st1=st1, st2=st2)
        self.assertEqual(cccoh, 1)
        self.assertEqual(i, 0)

    def test_inverted_coherence(self):
        """Reverse channels and ensure we get -1"""
        from eqcorrscan.utils.clustering import cross_chan_coherence
        from obspy import read
        st1 = read()
        st2 = read()
        for tr in st2:
            tr.data *= -1
        cccoh, i = cross_chan_coherence(st1=st1, st2=st2, i=1)
        self.assertEqual(cccoh, -1)
        self.assertEqual(i, 1)

    def test_known_coherence(self):
        """Test for a real stream case"""
        from eqcorrscan.utils.clustering import cross_chan_coherence
        from obspy import read
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        st1 = read(os.path.join(testing_path, 'WAV', 'TEST_',
                                '2013-09-01-0410-35.DFDPC_024_00'))
        st2 = read(os.path.join(testing_path, 'WAV', 'TEST_',
                                '2013-09-01-2040-11.DFDPC_039_00'))
        cccoh, i = cross_chan_coherence(st1=st1, st2=st2)
        self.assertTrue(cccoh < 0.01)
        self.assertEqual(i, 0)

    def test_distance_matrix(self):
        """Test that we can create a useful distance matrix."""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import distance_matrix
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_')
        stream_files = glob.glob(os.path.join(testing_path, '*'))[0:10]
        stream_list = [read(stream_file) for stream_file in stream_files]
        dist_mat = distance_matrix(stream_list=stream_list, cores=4)
        self.assertEqual(dist_mat.shape[0], len(stream_list))
        self.assertEqual(dist_mat.shape[1], len(stream_list))

    def test_unclustered(self):
        """Test clustering on unclustered data..."""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import cluster
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_')
        stream_files = glob.glob(os.path.join(testing_path, '*'))[0:10]
        stream_list = [(read(stream_file), i)
                       for i, stream_file in enumerate(stream_files)]
        groups = cluster(template_list=stream_list, show=False,
                         corr_thresh=0.3)
        self.assertEqual(len(groups), 10)  # They shouldn't cluster at all

    def test_clustered(self):
        """Test clustering on clustered data..."""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import cluster
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events')
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

    def test_delay_grouping(self):
        """Test grouping by delays"""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import group_delays
        from eqcorrscan.tutorials.template_creation import mktemplates
        from obspy.clients.fdsn.header import FDSNException
        import warnings

        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
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

    def test_SVD(self):
        """Test the SVD method."""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import SVD
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events')
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
        SVectors, SValues, Uvectors, stachans = SVD(stream_list=stream_list)
        self.assertEqual(len(SVectors), len(stachans))
        self.assertEqual(len(SValues), len(stachans))
        self.assertEqual(len(Uvectors), len(stachans))
        for SVec in SVectors:
            self.assertEqual(len(SVec), len(stream_list))

    def test_empirical_SVD(self):
        """Test the empirical SVD method"""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import empirical_SVD
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events')
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
        first_sub, second_sub = empirical_SVD(stream_list=stream_list)
        self.assertEqual(len(first_sub), len(second_sub))
        self.assertEqual(len(stream_list[0]), len(first_sub))

    def test_SVD_2_stream(self):
        """Test the conversion of SVD to stream."""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import SVD, SVD_2_stream

        samp_rate = 100
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events')
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
        SVectors, SValues, Uvectors, stachans = SVD(stream_list=stream_list)
        SVstreams = SVD_2_stream(SVectors=SVectors, stachans=stachans, k=4,
                                 sampling_rate=samp_rate)
        self.assertEqual(len(SVstreams), 4)

    def test_corr_cluster(self):
        """Test the corr_cluster function."""
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import corr_cluster

        samp_rate = 100
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        for stream in stream_list:
            for tr in stream:
                if not (tr.stats.station == 'WHAT2' and
                        tr.stats.channel == 'SH1'):
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
                tr.resample(sampling_rate=samp_rate)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
        trace_list = [stream[0] for stream in stream_list]
        output = corr_cluster(trace_list=trace_list, thresh=0.7)
        self.assertFalse(output.all())

    def test_dist_mat_km(self):
        """Test spacial clustering."""
        from eqcorrscan.utils.clustering import dist_mat_km
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        client = Client("IRIS")
        starttime = UTCDateTime("2002-01-01")
        endtime = UTCDateTime("2002-01-02")
        cat = client.get_events(starttime=starttime, endtime=endtime,
                                minmagnitude=6, catalog="ISC")
        dist_mat = dist_mat_km(cat)
        self.assertEqual(len(dist_mat), len(cat))

    def test_space_cluster(self):
        """Test the wrapper around dist_mat_km."""
        from eqcorrscan.utils.clustering import space_cluster
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        client = Client("IRIS")
        starttime = UTCDateTime("2002-01-01")
        endtime = UTCDateTime("2002-01-02")
        cat = client.get_events(starttime=starttime, endtime=endtime,
                                minmagnitude=6, catalog="ISC")
        groups = space_cluster(catalog=cat, d_thresh=1000, show=False)
        self.assertEqual(len([ev for group in groups for ev in group]),
                         len(cat))

    def test_space_time_cluster(self):
        """Test clustering in space and time."""
        from eqcorrscan.utils.clustering import space_time_cluster
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        client = Client("IRIS")
        starttime = UTCDateTime("2002-01-01")
        endtime = UTCDateTime("2002-01-02")
        cat = client.get_events(starttime=starttime, endtime=endtime,
                                minmagnitude=6, catalog="ISC")
        groups = space_time_cluster(catalog=cat, t_thresh=86400, d_thresh=1000)
        self.assertEqual(len([ev for group in groups for ev in group]),
                         len(cat))


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
