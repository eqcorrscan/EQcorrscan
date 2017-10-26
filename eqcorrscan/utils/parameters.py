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
import yaml

CORR_KEYS = [('n_stations', int), ('n_channels', int), ('data_len', float),
             ('n_templates', int), ('template_len', float),
             ('sampling_rate', float)]


class EQcorrscanConfig(object):
    """

    """
    def __init__(self, defaults={}):
        """"""
        self.defaults = {}
        # Try to read from ~/.eqcorrscan.yml
        config_file = os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc")
        if os.path.isfile(config_file):
            self.read(filename=config_file)
        try:
            self.defaults['correlation'].update(defaults['correlation'])
        except KeyError:
            self.defaults.update(defaults)

    def write(self, filename=None):
        """

        :param filename:
        :return:
        """
        if filename is None:
            filename = os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc")
        with open(filename, 'a') as f:
            yaml.dump(self.defaults, f)

    def read(self, filename=None):
        """

        :param filename:
        :return:
        """
        if filename is None:
            filename = os.path.join(os.path.expanduser("~"), ".eqcorrscan.rc")
            if not os.path.isfile(filename):
                raise IOError("%s not found" % filename)
        with open(filename, 'r') as f:
            try:
                config = yaml.load(stream=f)
            except:
                config = {}
        self.defaults.update(config)


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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
