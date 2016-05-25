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
        self.assertEqual(round(cccoh, 3), 0.003)
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



if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()