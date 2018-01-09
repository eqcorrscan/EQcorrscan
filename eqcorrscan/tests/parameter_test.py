"""Functions to test the parameter set-up
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest

import os

from eqcorrscan.utils.parameters import (
    _flatten_dataset_size, _inflate_dataset_size, EQcorrscanConfig,
    CorrelationDefaults, get_default_xcorr, ParameterError)


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
    @classmethod
    def setUpClass(cls):
        cls.config_filename = ".eqcorrscan.rc"

    def test_write_correlate_params(self):
        input_defaults = [CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr')]
        config = EQcorrscanConfig()
        config.defaults=input_defaults
        self.assertEqual(config.defaults, input_defaults)
        config.write(filename=self.config_filename)
        self.assertTrue(os.path.isfile(self.config_filename))
        config_back = EQcorrscanConfig(filename=self.config_filename)
        self.assertEqual(config.defaults, config_back.defaults)

    def test_append_correlate_params(self):
        input_defaults = [CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr')]
        config = EQcorrscanConfig()
        config.defaults = input_defaults
        config.write(filename=self.config_filename)
        # Add a second case
        input_defaults = [CorrelationDefaults(
            data_len=3600.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=100.0, template_len=2.0,
            corr_func='time_domain.multithread')]
        config2 = EQcorrscanConfig(defaults=input_defaults,
                                   filename=self.config_filename)
        self.assertTrue(len(config2.defaults) == 2)

    def test_unique_parameters(self):
        input_defaults = [CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr'),
            CorrelationDefaults(
                data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
                sampling_rate=200.0, template_len=2.0,
                corr_func='fftw.stream_xcorr'),
            CorrelationDefaults(
                data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
                sampling_rate=200.0, template_len=2.0,
                corr_func='time_domain.multithread'),
            CorrelationDefaults(
                data_len=6400.0, n_channels=3, n_stations=10, n_templates=5,
                sampling_rate=200.0, template_len=2.0,
                corr_func='time_domain.multithread')]
        config = EQcorrscanConfig()
        config.defaults = input_defaults
        self.assertEqual(len(config.defaults), 4)
        self.assertEqual(len(config.uniq().defaults), 2)

    def test_str_io(self):
        input_default = CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr')
        self.assertEqual(str(input_default), input_default.__repr__())
        self.assertEqual(
            str(input_default),
            "Correlation: {n_stations: 10, n_channels: 3, data_len: 86400.0, "
            "n_templates: 5, template_len: 2.0, sampling_rate: 200.0}: "
            "fftw.stream_xcorr")

    def test_equality(self):
        input_default = CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr')
        diff_default = CorrelationDefaults(
            data_len=6400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr')
        self.assertNotEqual(input_default, diff_default)
        self.assertNotEqual(input_default, float(2))

    def test_fail_deserialize(self):
        input_default = CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr')
        serial = input_default.serialize()
        serial = "NoneType " + serial
        with self.assertRaises(ParameterError):
            CorrelationDefaults().deserialize(serial)
        with open(self.config_filename, "w") as f:
            f.write(serial)
        config = EQcorrscanConfig(filename=self.config_filename)
        self.assertEqual(len(config.defaults), 0)

    def tearDown(self):
        if os.path.isfile(self.config_filename):
            os.remove(self.config_filename)


class TestDefaultLookup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config_filename = "eqcorrscan_test.rc"
        input_defaults = [CorrelationDefaults(
            data_len=86400.0, n_channels=3, n_stations=10, n_templates=5,
            sampling_rate=200.0, template_len=2.0,
            corr_func='fftw.stream_xcorr'),
            CorrelationDefaults(
                data_len=400.0, n_channels=3, n_stations=10, n_templates=50,
                sampling_rate=200.0, template_len=2.0,
                corr_func='time_domain.multithread'),
            CorrelationDefaults(
                data_len=900.0, n_channels=3, n_stations=10, n_templates=50,
                sampling_rate=200.0, template_len=2.0,
                corr_func='fftw.stream_xcorr'),
            CorrelationDefaults(
                data_len=3600.0, n_channels=1, n_stations=20, n_templates=200,
                sampling_rate=100.0, template_len=6.0,
                corr_func='fftw.stream_xcorr'),
            CorrelationDefaults(
                data_len=60.0, n_channels=1, n_stations=1, n_templates=50,
                sampling_rate=200.0, template_len=2.0,
                corr_func='time_domain.multiprocess'),
            CorrelationDefaults(
                data_len=6.0, n_channels=1, n_stations=1, n_templates=50,
                sampling_rate=200.0, template_len=2.0,
                corr_func='walrus.fast')]
        config = EQcorrscanConfig(defaults=input_defaults)
        config.write(filename=cls.config_filename, append=False)

    def test_set_with_empty_datasize(self):
        xcorr = get_default_xcorr(fname=self.config_filename)
        self.assertEqual(xcorr.__name__, "_fftw_stream_xcorr")

    def test_set_with_exact_match(self):
        xcorr = get_default_xcorr(
            data_len=400.0, n_channels=3, n_stations=10, n_templates=50,
            sampling_rate=200.0, template_len=2.0, fname=self.config_filename)
        self.assertEqual(xcorr.__name__, "_time_threaded_normxcorr")

    def test_fail_with_exact_match(self):
        """Default if correlation function not registered"""
        xcorr = get_default_xcorr(
            data_len=6.0, n_channels=1, n_stations=1, n_templates=50,
            sampling_rate=200.0, template_len=2.0, fname=self.config_filename)
        self.assertEqual(xcorr.__name__, "_fftw_stream_xcorr")

    def test_set_with_near_match(self):
        xcorr = get_default_xcorr(
            data_len=50.0, n_channels=1, n_stations=1, n_templates=50,
            sampling_rate=200.0, template_len=2.0, fname=self.config_filename)
        self.assertEqual(xcorr.__qualname__,
                         "_general_multiprocess.<locals>.multiproc")

    def test_fail_with_near_match(self):
        """Default if correlation function not registered"""
        xcorr = get_default_xcorr(
            data_len=5.0, n_channels=1, n_stations=1, n_templates=50,
            sampling_rate=200.0, template_len=2.0, fname=self.config_filename)
        self.assertEqual(xcorr.__name__, "_fftw_stream_xcorr")

    def test_set_with_no_match(self):
        xcorr = get_default_xcorr(
            data_len=86400.0, n_channels=3, n_stations=100, n_templates=2000,
            sampling_rate=200.0, template_len=2.0, fname=self.config_filename)
        self.assertEqual(xcorr.__name__, "_fftw_stream_xcorr")

    def test_set_without_file(self):
        with self.assertRaises(OSError):
            get_default_xcorr(
                data_len=400.0, n_channels=3, n_stations=10, n_templates=50,
                sampling_rate=200.0, template_len=2.0, fname="walrus.rc")

    @classmethod
    def tearDownClass(cls):
        if os.path.isfile(cls.config_filename):
            os.remove(cls.config_filename)


if __name__ == '__main__':
    unittest.main()
