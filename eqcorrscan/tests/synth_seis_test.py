"""
Functions to test the synth_seis generation.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import numpy as np
from obspy import Stream

from eqcorrscan.utils.synth_seis import seis_sim, SVD_sim, template_grid
from eqcorrscan.utils.synth_seis import generate_synth_data


class TestSynth(unittest.TestCase):
    def test_seis_sim(self):
        """
        Check that phases are put in the right place
        """
        for sp_samples in [10, 15, 20, 50]:
            sim_data = seis_sim(sp=sp_samples)
            # Find first non-zero sample - start of P-phase
            for i, sample in enumerate(sim_data):
                if sample != 0.0:
                    start_sample = i
                    break
            # Find maximum - start of S-phase
            max_sample = np.argmax(np.abs(sim_data))
            self.assertEqual((max_sample - start_sample) - 1, sp_samples)

    def test_fixed_lengths(self):
        for sp_samples in [10, 15, 20]:
            sim_data = seis_sim(sp=sp_samples, flength=150)
            # Find first non-zero sample - start of P-phase
            for i, sample in enumerate(sim_data):
                if sample != 0.0:
                    start_sample = i
                    break
            # Find maximum - start of S-phase
            max_sample = np.argmax(np.abs(sim_data))
            self.assertEqual((max_sample - start_sample) - 1, sp_samples)

    def test_phaseout(self):
        sp_samples = 15
        flength = 150
        sim_data = seis_sim(sp=sp_samples, flength=flength, phaseout='S')
        self.assertEqual(len(sim_data), flength)


class TestSVDSim(unittest.TestCase):
    def test_svd_sim(self):
        U, s, V, stachans = SVD_sim(sp=15, lowcut=2, highcut=8, samp_rate=20)
        self.assertEqual(V[0].shape[0], s[0].shape[0])
        self.assertEqual(V[0].shape[0], U[0].shape[1])
        self.assertEqual(len(stachans), 1)


class TestTemplateGrid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        lats = np.random.random(20) * 90.0
        lons = np.random.random(20) * 90.0
        depths = np.abs(np.random.random(20) * 40.0)
        cls.stations = ['ALPH', 'BETA', 'GAMM', 'KAPP', 'ZETA', 'BOB', 'MAGG',
                        'ALF', 'WALR', 'ALBA', 'PENG', 'BANA', 'WIGG', 'SAUS',
                        'MALC']
        cls.nodes = list(zip(lats, lons, depths))
        cls.travel_times = np.abs(
            np.random.random([len(cls.stations), 20])) * 100

    def test_template_grid(self):
        templates = template_grid(stations=self.stations, nodes=self.nodes,
                                  travel_times=self.travel_times, phase='P')
        self.assertEqual(len(templates), 20)
        templates = template_grid(stations=self.stations, nodes=self.nodes,
                                  travel_times=self.travel_times, phase='S')
        self.assertEqual(len(templates), 20)
        with self.assertRaises(IOError):
            template_grid(stations=self.stations, nodes=self.nodes,
                          travel_times=self.travel_times, phase='bob')
        templates = template_grid(stations=self.stations, nodes=self.nodes,
                                  travel_times=self.travel_times, phase='S',
                                  phaseout='S')
        self.assertEqual(len(templates), 20)
        templates = template_grid(stations=self.stations, nodes=self.nodes,
                                  travel_times=self.travel_times, phase='S',
                                  phaseout='all', flength=100)
        self.assertEqual(len(templates), 20)
        templates = template_grid(stations=self.stations, nodes=self.nodes,
                                  travel_times=self.travel_times, phase='S',
                                  phaseout='both', flength=100)
        self.assertEqual(len(templates), 20)


class TestRandomData(unittest.TestCase):
    def test_generate_synth_dataset(self):
        for debug in [0, 2, 3]:
            templates, data, seeds = generate_synth_data(
                nsta=2, ntemplates=2, nseeds=2, samp_rate=100, t_length=10,
                max_amp=10, max_lag=20, debug=debug)
            self.assertEqual(len(templates), 2)
            self.assertEqual(len(seeds[0]['time']), 2)
            self.assertTrue(isinstance(data, Stream))


if __name__ == '__main__':
    unittest.main()
