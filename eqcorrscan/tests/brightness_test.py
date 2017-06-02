"""
A series of test functions for the core.bright_lights module in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import numpy as np
import os
import shutil

from obspy import Trace, Stream, read
from matplotlib import path

from eqcorrscan.core.bright_lights import brightness, _read_tt, _resample_grid
from eqcorrscan.core.bright_lights import _rm_similarlags, _rms, _node_loop
from eqcorrscan.core.bright_lights import _cum_net_resp, _find_detections
from eqcorrscan.core.bright_lights import coherence, BrightnessError
from eqcorrscan.utils.synth_seis import generate_synth_data


class BrightnessTestMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testing_path = os.path.join(os.path.abspath(
            os.path.dirname(__file__)), 'test_data') + os.sep

    def test_read_tt(self):
        # Test reading S from S
        stations, nodes, lags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='S', phaseout='S')
        self.assertEqual(stations[0], 'COSA')
        self.assertEqual(len(nodes), len(lags[0]))
        # Test reading P from S
        stations, nodes, lags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='S', phaseout='P')
        self.assertEqual(stations[0], 'COSA')
        self.assertEqual(len(nodes), len(lags[0]))
        # Test reading S from P
        stations, nodes, lags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='P', phaseout='S')
        self.assertEqual(stations[0], 'COSA')
        self.assertEqual(len(nodes), len(lags[0]))
        # Test reading P from P
        stations, nodes, lags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='P', phaseout='P')
        self.assertEqual(stations[0], 'COSA')
        self.assertEqual(len(nodes), len(lags[0]))
        # Test lag_switch
        stations, nodes, lags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='P', phaseout='P',
            lags_switch=False)
        self.assertEqual(stations[0], 'COSA')
        self.assertEqual(len(nodes), len(lags[0]))

    def test_resample_grid(self):
        minlon = 168
        maxlon = 170
        minlat = -46
        maxlat = -43
        mindepth = 4
        maxdepth = 10
        stations, allnodes, alllags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='S', phaseout='S')
        corners = [(minlon, minlat),
                   (maxlon, minlat),
                   (maxlon, maxlat),
                   (minlon, maxlat)]
        corners = path.Path(corners, closed=True)
        stations, nodes, lags = _resample_grid(
            stations, allnodes, alllags, mindepth=mindepth, maxdepth=maxdepth,
            corners=corners)
        for node in nodes:
            self.assertTrue(minlon < node[0] < maxlon)
            self.assertTrue(minlat < node[1] < maxlat)
            self.assertTrue(mindepth < node[2] < maxdepth)

        for node in allnodes:
            if node not in nodes:
                self.assertFalse(
                    (minlon < node[0] < maxlon) and
                    (minlat < node[1] < maxlat) and
                    (mindepth < node[2] < maxdepth))

    def test_rm_similarlags(self):
        threshold = 2
        stations, allnodes, alllags = _read_tt(
            path=self.testing_path, stations=['COSA'], phase='S', phaseout='S')
        stations, nodes, lags = _rm_similarlags(
            stations=stations, nodes=allnodes, lags=alllags,
            threshold=threshold)
        for lag in lags:
            for _lag in lag:
                other_lags = np.array([l for l in lag if not l == _lag])
                self.assertTrue(np.all(np.abs(other_lags - _lag) > threshold))

    def test_rms(self):
        rms = _rms(np.zeros(1000) + 1)
        self.assertEqual(rms, 1)
        rms = _rms(np.random.randn(10000))
        self.assertEqual(round(rms), 1)

    def test_node_loop(self):
        stations, nodes, lags = _read_tt(path=self.testing_path,
                                         stations=['COSA', 'LABE'],
                                         phase='S', phaseout='S')
        st = Stream(Trace())
        st[0].stats.station = stations[0]
        st[0].data = np.random.randn(86400) * 3000
        st[0].data = st[0].data.astype(np.int16)
        st += Trace(np.random.randn(86400) * 3000)
        st[1].stats.station = 'LABE'
        index, energy = _node_loop(stations=stations, lags=lags[:, 1],
                                   stream=st, clip_level=4)
        self.assertEqual(index, 0)
        self.assertEqual(np.shape(energy), (1, 86400))

    def test_cum_net_resp(self):
        stations, nodes, lags = _read_tt(path=self.testing_path,
                                         stations=['COSA', 'LABE'],
                                         phase='S', phaseout='S')
        st = Stream(Trace())
        st[0].stats.station = stations[0]
        st[0].data = np.random.randn(86400) * 3000
        st[0].data = st[0].data.astype(np.int16)
        st += Trace(np.random.randn(86400) * 3000)
        st[1].stats.station = stations[1]
        if not os.path.isdir('tmp0'):
            os.makedirs('tmp0')
        index, energy_file = _node_loop(stations=stations, lags=lags[:, 0],
                                        stream=st, clip_level=4,
                                        mem_issue=True, instance=0)
        self.assertEqual(index, 0)
        self.assertTrue(type(energy_file) == str)
        _node_loop(stations=stations, lags=lags[:, 1], stream=st, clip_level=4,
                   mem_issue=True, i=1, instance=0)
        cum_net_resp, indeces = _cum_net_resp(node_lis=[0, 1], instance=0)
        self.assertEqual(len(cum_net_resp), 86400)
        self.assertEqual(len(indeces), 86400)
        shutil.rmtree('tmp0')

    def test_find_detections(self):
        stations, nodes, lags = _read_tt(path=self.testing_path,
                                         stations=['COSA', 'LABE'],
                                         phase='S', phaseout='S')
        st = Stream(Trace())
        st[0].stats.station = stations[0]
        st[0].data = np.random.randn(86400) * 3000
        st[0].data = st[0].data.astype(np.int16)
        st += Trace(np.random.randn(86400) * 3000)
        st[1].stats.station = stations[1]
        # Need to make a tmp0 folder
        os.makedirs('tmp1')
        _node_loop(stations=stations, lags=lags[:, 1], stream=st, clip_level=4,
                   mem_issue=True, instance=1)
        cum_net_resp, indeces = _cum_net_resp(node_lis=[0], instance=1)
        all_nodes = [nodes[1] for i in range(len(cum_net_resp))]
        detections = _find_detections(cum_net_resp=cum_net_resp,
                                      nodes=all_nodes, threshold=10,
                                      thresh_type='MAD', samp_rate=1,
                                      realstations=[tr.stats.station
                                                    for tr in st],
                                      length=10)
        self.assertEqual(len(detections), 0)
        detections = _find_detections(cum_net_resp=cum_net_resp,
                                      nodes=all_nodes, threshold=5,
                                      thresh_type='MAD', samp_rate=1,
                                      realstations=[tr.stats.station
                                                    for tr in st],
                                      length=10)
        self.assertTrue(len(detections) > 0)
        detections = _find_detections(cum_net_resp=cum_net_resp,
                                      nodes=all_nodes, threshold=5,
                                      thresh_type='abs', samp_rate=1,
                                      realstations=[tr.stats.station
                                                    for tr in st],
                                      length=10)
        self.assertTrue(len(detections) > 0)
        detections = _find_detections(cum_net_resp=cum_net_resp,
                                      nodes=all_nodes, threshold=5,
                                      thresh_type='RMS', samp_rate=1,
                                      realstations=[tr.stats.station
                                                    for tr in st],
                                      length=10)
        self.assertTrue(len(detections) > 0)
        shutil.rmtree('tmp1')

    def test_coherence(self):
        testing_path = os.path.join(self.testing_path + 'WAV', 'TEST_',
                                    '2013-09-01-0410-35.DFDPC_024_00')
        st = read(testing_path)
        coh = coherence(stream_in=st)
        self.assertTrue(type(coh), float)
        coh = coherence(stream_in=st, clip=(0.5, 10))
        self.assertTrue(type(coh), float)
        coh = coherence(stream_in=st, stations=[tr.stats.station
                                                for tr in st[0:-5]])
        self.assertTrue(type(coh), float)


class TestBrightnessMain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data') + os.sep
        # Test reading S from S
        cls.stations, cls.nodes, cls.lags = _read_tt(
            path=testing_path, stations=['COSA', 'LABE'],  phase='S',
            phaseout='S')
        templates, cls.st, seeds = generate_synth_data(
            nsta=2, ntemplates=2, nseeds=20, samp_rate=1, t_length=15,
            max_amp=10, max_lag=6)
        cls.st[0].stats.station = cls.stations[0]
        cls.st[1].stats.station = cls.stations[1]
        cls.st[0].data = cls.st[0].data.astype(np.int16)
        cls.st[1].data = cls.st[1].data.astype(np.int16)
        cls.st[1].stats.channel = 'HHZ'
        cls.st[0].stats.channel = 'HHZ'

    def test_brightness(self):
        st = self.st.copy()
        detections, nodes_out = brightness(
            stations=self.stations, nodes=self.nodes, lags=self.lags,
            stream=st, threshold=1.885, thresh_type='MAD', template_length=1,
            template_saveloc='.', coherence_thresh=(10, 1), cores=1)
        self.assertEqual(len(detections), 0)
        self.assertEqual(len(detections), len(nodes_out))

    def test_big_data_range(self):
        st = self.st.copy()
        for tr in st:
            tr.data *= 40000
        detections, nodes_out = brightness(
            stations=self.stations, nodes=self.nodes, lags=self.lags,
            stream=st, threshold=1.885, thresh_type='MAD', template_length=1,
            template_saveloc='.', coherence_thresh=(10, 1), cores=1)
        self.assertEqual(len(detections), 0)
        self.assertEqual(len(detections), len(nodes_out))

    def test_fail_too_many_traces(self):
        st = self.st.copy()
        for i in range(130):
            st += Trace(np.zeros(1))

        with self.assertRaises(BrightnessError):
            brightness(stations=self.stations, nodes=self.nodes,
                       lags=self.lags, stream=st, threshold=1.885,
                       thresh_type='MAD', template_length=1,
                       template_saveloc='.', coherence_thresh=(10, 1))

    def test_mem_issue(self):
        st = self.st.copy()
        detections, nodes_out = brightness(
            stations=self.stations, nodes=self.nodes, lags=self.lags,
            stream=st, threshold=1.885, thresh_type='MAD', template_length=1,
            template_saveloc='.', coherence_thresh=(10, 1), cores=1,
            mem_issue=True, instance=2)
        self.assertEqual(len(detections), 0)
        self.assertEqual(len(detections), len(nodes_out))
        shutil.rmtree('tmp2')

    def test_mem_issue_parallel(self):
        st = self.st.copy()
        detections, nodes_out = brightness(
            stations=self.stations, nodes=self.nodes, lags=self.lags,
            stream=st, threshold=1.885, thresh_type='MAD', template_length=1,
            template_saveloc='.', coherence_thresh=(10, 1), cores=3,
            mem_issue=True, instance=3)
        self.assertEqual(len(detections), 0)
        self.assertEqual(len(detections), len(nodes_out))
        shutil.rmtree('tmp3')


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
