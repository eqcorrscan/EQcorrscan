"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestTutorialScripts(unittest.TestCase):
    def test_match_filter(self):
        """Test the match_filter tutorial, generates templates too."""
        from eqcorrscan.tutorials.template_creation import mktemplates
        from eqcorrscan.tutorials.match_filter import run_tutorial
        from eqcorrscan.core.match_filter import read_detections
        import os
        import glob
        from obspy import read

        # Run mktemplates first to set-up for match_filter
        mktemplates(plot=False)
        for template_no in range(4):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            self.assertTrue(len(template) > 1)
        del(template)
        # Run the matched-filter
        tutorial_detections = run_tutorial(plot=False)
        # It should make 20 detections in total...
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        fname = os.path.join(testing_path,
                             'expected_tutorial_detections.txt')
        expected_detections = read_detections(fname)

        # Annoyingly something doesn't match, event when writing out detections
        # then reading them back in and comparing them to themselves in memory.
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
        # Cleanup the templates
        templates = glob.glob('tutorial_template_?.ms')
        for template in templates:
            os.remove(template)

    def test_lag_calc(self):
        """Test the lag calculation tutorial."""
        from eqcorrscan.tutorials.lag_calc import run_tutorial

        shift_len = 0.2
        min_mag = 4
        detections, picked_catalog, templates, template_names = \
            run_tutorial(min_magnitude=min_mag, shift_len=shift_len)

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
                # Check that ccs increases
                # for comment in event.comments:
                #     if comment.text.split('=')[0] == 'detect_val':
                #         post_lag_ccs = float(comment.text.split('=')[-1])
                # for comment in detection.event.comments:
                #     if comment.text.split('=')[0] == 'detect_val':
                #         pre_lag_ccs = float(comment.text.split('=')[-1])
                # self.assertGreaterEqual(post_lag_ccs, pre_lag_ccs)

    def test_subspace(self):
        """Test the subspace tutorial."""
        from eqcorrscan.tutorials.subspace import run_tutorial

        detections = run_tutorial(plot=False)
        self.assertEqual(len(detections), 2)

if __name__ == '__main__':
    """
    Run tutorial tests
    """
    unittest.main()
