"""
Tests for the eqcorrscan.utils.picker functions.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import os

from obspy import read

from eqcorrscan.utils.picker import stalta_pick, cross_net


class PickerTests(unittest.TestCase):
    def test_stalta_picker(self):
        """
        Tests the stalta picker.
        """
        st = read()
        event = stalta_pick(stream=st, stalen=0.2, ltalen=3, trig_on=6,
                            trig_off=1.2)
        self.assertEqual(len(event.picks), 2)
        # Test a lower threshold
        event = stalta_pick(stream=st, stalen=0.2, ltalen=3, trig_on=3,
                            trig_off=1.2)
        self.assertEqual(len(event.picks), 11)
        # Test a higher threshold
        event = stalta_pick(stream=st, stalen=0.2, ltalen=3, trig_on=10,
                            trig_off=1.2)
        self.assertEqual(len(event.picks), 0)

    def test_cross_net(self):
        """
        Test the WECC style picker
        """
        st = read(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                               'test_data', 'WAV', 'TEST_',
                               '2013-09-01-0410-35.DFDPC_024_00'))
        event = cross_net(stream=st)
        self.assertEqual(len(event.picks), 24)
        event = cross_net(stream=st, env=False, debug=3)
        self.assertEqual(len(event.picks), 24)
