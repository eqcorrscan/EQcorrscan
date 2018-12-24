"""
IO for eqcorrscan config files.
:copyright:
    EQcorrscan developers.
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
from collections import Counter
import numpy as np

from eqcorrscan.utils.correlate import get_stream_xcorr

CORR_KEYS = [('n_stations', int), ('n_channels', int), ('data_len', float),
             ('n_templates', int), ('template_len', float),
             ('sampling_rate', float)]


class ParameterError(Exception):
    """ Error handling for parameter functions. """

    def __init__(self, value):
        """ Raise error.
        .. rubric:: Example
        >>> ParameterError("Something happened")
        Something happened
        """
        self.value = value

    def __repr__(self):
        """ Print error value.
        .. rubric:: Example
        >>> print(ParameterError("Error").__repr__())
        Error
        """
        return self.value

    def __str__(self):
        """ Print otherwise
        .. rubric:: Example
        >>> print(ParameterError("Error"))
        Error
        """
        return self.value


"""
Note: All subclasses of EQcorrscanConfig must define at-least the following:
  - __eq__
  - __ne__
  - serialize
  - deserialize
"""


class CorrelationDefaults(object):
    """
    Holder for correlation parameter defaults.
    """
    def __init__(self, n_stations=None, n_channels=None, data_len=None,
                 n_templates=None, template_len=None, sampling_rate=None,
                 corr_func=None):
        """
        Configuration for correlation function for one dataset size
        """
        self.n_stations = n_stations
        self.n_channels = n_channels
        self.data_len = data_len
        self.n_templates = n_templates
        self.template_len = template_len
        self.sampling_rate = sampling_rate
        self.corr_func = corr_func

    def __str__(self):
        return self.serialize()

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, CorrelationDefaults):
            return False
        if self.dataset_size() == other.dataset_size():
            return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def _isempty(self):
        for value in self.__dict__.values():
            if value is not None:
                return False
        return True

    def dataset_size(self):
        d_size = self.__dict__.copy()
        d_size.pop('corr_func')
        return d_size

    def serialize(self):
        d_size = _flatten_dataset_size(self.dataset_size())
        return "Correlation: {" + d_size + "}: " + self.corr_func

    def deserialize(self, serial):
        if str(serial.lstrip()[0:12]) != str("Correlation:"):
            raise ParameterError("Config line is not for correlation.")
        serial = serial.split('}:')
        self.corr_func = serial[-1].strip()
        serial = serial[0].split("Correlation: {")[1]
        self.__dict__.update(_inflate_dataset_size(serial))
        return self


class EQcorrscanConfig(object):
    """ Holder for all EQcorrscan config values. """
    def __init__(self, defaults=[], filename=None):
        """
        Configuration parameter storage for EQcorrscan
        :type defaults: list
        :param defaults: List of default parameters
        :type filename: str
        :param filename:
            File to read defaults from. If defaults are given as an argument
            given defaults will overwrite those in this file.
        """
        self.defaults = []
        self.read(filename=filename)
        self.defaults.extend(defaults)

    def write(self, filename=None, append=True):
        """
        Write parameters to file. If filename is None will use ~/.eqcorrscan.rc
        :type filename: str
        :param filename: filename to write parameters to
        :type append: bool
        :param append: Whether to append to old parameters or overwrite
        """
        if filename is None:
            filename = os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc")
        if append:
            mode = 'a'
        else:
            mode = 'w'
        with open(filename, mode) as f:
            for default in self.defaults:
                f.write(default.serialize() + '\n')

    def read(self, filename=None):
        """
        Read parameters from a given file.
        :type filename: str
        :param filename:
            File to read from. If not specified will read from ~/.eqcorrscan.rc
        """
        if filename is not None:
            if not os.path.isfile(filename):
                raise IOError("%s not found" % filename)
        else:
            filename = os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc")
            if not os.path.isfile(filename):
                return
        with open(filename, 'r') as f:
            for line in f:
                try:
                    default = CorrelationDefaults().deserialize(line)
                except ParameterError as e:
                    print(e)
                    continue
                if default is not None:
                    self.defaults.append(default)

    def get(self, parameter_type):
        """
        Get a subset of parameters by config type.
        :type parameter_type: str
        :param parameter_type: Type of configuration parameter
        :return: list of configuration parameters
        """
        out = []
        type_map = {"Correlation": CorrelationDefaults,
                    "correlation": CorrelationDefaults,
                    "CorrelationDefaults": CorrelationDefaults}
        for default in self.defaults:
            if isinstance(default, type_map[parameter_type]):
                out.append(default)
        return out

    def uniq(self):
        """
        Remove duplicate parameters.
        """
        uniq_parameters = []
        for default in self.defaults:
            if default not in uniq_parameters:
                uniq_parameters.append(default)
        self.defaults = uniq_parameters
        return self


def _flatten_dataset_size(dataset_size):
    """
    Flatten a dict of dataset parameters to use as a key.
    :type dataset_size: dict
    :return: str
    """
    flat_str = ''
    for key in dataset_size.keys():
        if key not in dict(CORR_KEYS).keys():
            raise KeyError("Key %s is not in the standard key set" % key)
    for key, dtype in CORR_KEYS:
        if key not in dataset_size:
            raise KeyError("Key %s is missing for dataset_size" % key)
        flat_str += ': '.join([key, str(dataset_size[key])]) + ', '
    return flat_str.rstrip(', ')


def _inflate_dataset_size(flat_str):
    """
    Convert flat string to dictionary
    :type flat_str: str
    :return: dict
    """
    dataset_size = {}
    flat_list = flat_str.split(', ')
    for item in flat_list:
        dataset_size.update({item.split(': ')[0]: item.split(': ')[1]})
    for key, item in CORR_KEYS:
        if key not in dataset_size.keys():
            raise KeyError("Key %s missing from flat_str." % key)
        dataset_size.update({key: item(dataset_size[key])})
    return dataset_size


def _scalar_size(dataset_size):
    scalar_size = 1
    for key, value in dataset_size.items():
        if value is None:
            return np.inf
        scalar_size *= value
    return scalar_size


def get_default_xcorr(fname=None, n_stations=None, n_channels=None,
                      data_len=None, n_templates=None, template_len=None,
                      sampling_rate=None, scalar_gap=1000000):
    default_xcorr = get_stream_xcorr()
    config = EQcorrscanConfig(filename=fname)
    config = config.get("CorrelationDefaults")
    if len(config) == 0:
        return default_xcorr
    lookup_dataset = CorrelationDefaults(
        n_stations=n_stations, n_channels=n_channels, data_len=data_len,
        n_templates=n_templates, template_len=template_len,
        sampling_rate=sampling_rate)
    # If all are None, we should give the most common "best" option
    if lookup_dataset._isempty():
        methods = Counter([default.corr_func for default in config])
        most_common = methods.most_common(1)[0][0]
        try:
            default_xcorr = get_stream_xcorr(
                name_or_func=most_common.split('.')[0],
                concurrency=most_common.split('.')[1])
        except (ValueError, KeyError):
            pass
    else:
        # Otherwise, look for exact matches
        for default in config:
            if default.dataset_size() == lookup_dataset.dataset_size():
                try:
                    default_xcorr = get_stream_xcorr(
                        name_or_func=default.corr_func.split('.')[0],
                        concurrency=default.corr_func.split('.')[1])
                    break
                except (ValueError, KeyError):
                    pass
        else:
            # If there is no match, find the closest scalar size
            scalar_gaps = [
                abs(_scalar_size(lookup_dataset.dataset_size()) -
                    _scalar_size(default.dataset_size()))
                for default in config]
            min_gap = scalar_gaps.index(min(scalar_gaps))
            if scalar_gaps[min_gap] < scalar_gap:
                try:
                    default_xcorr = get_stream_xcorr(
                        name_or_func=config[min_gap].corr_func.split('.')[0],
                        concurrency=config[min_gap].corr_func.split('.')[1])
                except (ValueError, KeyError):
                    pass
    return default_xcorr


if __name__ == '__main__':
    import doctest
    doctest.testmod()
