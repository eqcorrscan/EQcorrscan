"""Functions to test the parameter set-up
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest

import os

from eqcorrscan.utils.parameters import (
    _flatten_dataset_size, _inflate_dataset_size, EQcorrscanConfig)


class TestDictFlatten(unittest.TestCase):
    """Test setup and read/write parameter files"""
    def test_flatten_dict(self):
        dataset_size = {
            'data_len': 86400.0, 'n_channels': 3, 'n_stations': 10,
            'n_templates': 5, 'sampling_rate': 200.0, 'template_len': 2.0}
        flat_str = _flatten_dataset_size(dataset_size=dataset_size)
        self.assertEqual(
            flat_str,  'n_stations: 10, n_channels: 3, data_len: 86400.0, n_t'
                       'emplates: 5, template_len: 2.0, sampling_rate: 200.0')

    def test_flatten_fail(self):
        dataset_size = {
            'data_len': 86400.0, 'n_channels': 3, 'n_stations': 10,
            'n_templates': 5, 'sampling_rate': 200.0}
        with self.assertRaises(KeyError):
            _flatten_dataset_size(dataset_size)
        dataset_size.update({'template_len': 2.0, 'bob': 'wilf'})
        with self.assertRaises(KeyError):
            _flatten_dataset_size(dataset_size)

    def test_inflate(self):
        flat_str = ('n_channels: 3, data_len: 86400.0, n_templates: 5, '
                    'sampling_rate: 200.0, n_stations: 10, template_len: 2.0')
        dataset_size = _inflate_dataset_size(flat_str=flat_str)
        self.assertEqual(dataset_size, {
            'data_len': 86400.0, 'n_channels': 3, 'n_stations': 10,
            'n_templates': 5, 'sampling_rate': 200.0, 'template_len': 2.0})


class TestConfigIO(unittest.TestCase):
    def test_write_correlate_params(self):
        dataset_size = {
            'data_len': 86400.0, 'n_channels': 3, 'n_stations': 10,
            'n_templates': 5, 'sampling_rate': 200.0, 'template_len': 2.0}
        flat_str = _flatten_dataset_size(dataset_size=dataset_size)
        input_defaults = {'correlation': {flat_str: 'fftw'}}
        config = EQcorrscanConfig(defaults=input_defaults)
        self.assertEqual(config.defaults, input_defaults)
        config.write()
        self.assertTrue(os.path.isfile(os.path.join(os.path.expanduser("~"),
                                                    ".eqcorrscan.rc")))
        config_back = EQcorrscanConfig()
        self.assertEqual(config.defaults, config_back.defaults)

    def test_append_correlate_params(self):
        dataset_size = {
            'data_len': 86400.0, 'n_channels': 3, 'n_stations': 10,
            'n_templates': 5, 'sampling_rate': 200.0, 'template_len': 2.0}
        flat_str = _flatten_dataset_size(dataset_size=dataset_size)
        input_defaults = {'correlation': {flat_str: 'fftw'}}
        config = EQcorrscanConfig(defaults=input_defaults)
        config.write()
        # Add a second case
        dataset_size = {
            'data_len': 3600.0, 'n_channels': 3, 'n_stations': 10,
            'n_templates': 5, 'sampling_rate': 100.0, 'template_len': 2.0}
        flat_str2 = _flatten_dataset_size(dataset_size=dataset_size)
        input_defaults = {'correlation': {flat_str2: 'time-domain'}}
        config2 = EQcorrscanConfig(defaults=input_defaults)
        self.assertEqual(config.defaults['correlation'][flat_str],
                         config2.defaults['correlation'][flat_str])
        self.assertTrue(len(config2.defaults['correlation'].keys()) == 2)

    def tearDown(self):
        if os.path.isfile(
                os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc")):
            os.remove(
                os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc"))


if __name__ == '__main__':
    unittest.main()
