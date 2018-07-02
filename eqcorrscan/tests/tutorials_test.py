"""
Functions for testing the tutorials - written as a somewhat monolithic test
because we need certain behaviour.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import os
import pytest

from obspy import read

from eqcorrscan.tutorials.template_creation import mktemplates
from eqcorrscan.tutorials import match_filter, lag_calc, subspace
from eqcorrscan.core.match_filter import read_detections


slow = pytest.mark.skipif(
    not pytest.config.getoption("--runslow"),
    reason="need --runslow option to run"
)


superslow = pytest.mark.skipif(
    not pytest.config.getoption("--runsuperslow"),
    reason="need --runsuperslow option to run"
)


class TestTutorialScripts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')

    @slow
    @pytest.mark.flaky(reruns=2)
    def test_templates_and_match(self):
        """Call the template creation then the matched-filter tests."""
        mktemplates(plot=False)
        for template_no in range(4):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            expected_template = read(
                os.path.join(self.testing_path,
                             'tutorial_template_' + str(template_no) + '.ms'))
            self.assertTrue(len(template) > 1)
            # self.assertEqual(template, expected_template)
        # Run the matched-filter
        tutorial_detections = match_filter.run_tutorial(plot=False)
        # It should make 20 detections in total...
        fname = os.path.join(self.testing_path,
                             'expected_tutorial_detections.txt')
        expected_detections = read_detections(fname)

        expected_times = [detection.detect_time for detection
                          in expected_detections]
        for expected_time in expected_times:
            expected_time.precision = 3  # Lower the precision slightly
        # expected_correlations = [round(detection.detect_val, 4) for detection
        #                          in expected_detections]
        for detection in tutorial_detections:
            detection.detect_time.precision = 3
            self.assertIn(detection.detect_time, expected_times,
                          msg='Detection at %s is not in expected detections'
                          % detection.detect_time)
        if len(expected_detections) > len(tutorial_detections):
            # This is a fail but we are trying to debug
            actual_times = [tutorial_detection.detect_time
                            for tutorial_detection in tutorial_detections]
            for detection in expected_detections:
                self.assertIn(detection.detect_time, actual_times,
                              msg='Expected detection at %s was not made'
                              % detection.detect_time)
        self.assertEqual(len(tutorial_detections), 22)
        for template_no in range(4):
            if os.path.isfile('tutorial_template_' +
                              str(template_no) + '.ms'):
                os.remove('tutorial_template_' + str(template_no) + '.ms')

    @slow
    @pytest.mark.flaky(reruns=2)
    def test_lag_calc(self):
        """Test the lag calculation tutorial."""
        shift_len = 0.2
        min_mag = 4
        detections, picked_catalog, templates, template_names = \
            lag_calc.run_tutorial(min_magnitude=min_mag, shift_len=shift_len)

        self.assertEqual(len(picked_catalog), len(detections))
        self.assertEqual(len(detections), 8)
        for event, detection in zip(picked_catalog, detections):
            template = [t[0] for t in zip(templates, template_names)
                        if t[1] == detection.template_name][0]
            template_stachans = [(tr.stats.station, tr.stats.channel)
                                 for tr in template]
            for pick in event.picks:
                # First check that there is a template for the pick
                stachan = (pick.waveform_id.station_code,
                           pick.waveform_id.channel_code)
                self.assertTrue(stachan in template_stachans)
                # Now check that the pick time is within +/- shift_len of
                # The template
                tr = template.select(station=stachan[0], channel=stachan[1])[0]
                delay = tr.stats.starttime - \
                    template.sort(['starttime'])[0].stats.starttime
                re_picked_delay = pick.time - (detection.detect_time + delay)
                self.assertTrue(abs(re_picked_delay) < shift_len)

    @superslow
    @pytest.mark.flaky(reruns=2)
    def test_subspace(self):
        """Test the subspace tutorial."""
        detections = subspace.run_tutorial(plot=False)
        self.assertEqual(len(detections), 11)


if __name__ == '__main__':
    """
    Run tutorial tests
    """
    unittest.main()
