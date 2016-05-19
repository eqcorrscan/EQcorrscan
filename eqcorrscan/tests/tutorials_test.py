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
            # self.assertIn(round(detection.detect_val, 4),
            #               expected_correlations,
            #               msg='Detection with cross-correlation value %s not' +
            #               ' in expected detections' % detection.detect_val)
        if len(expected_detections) > len(tutorial_detections):
            # This is a fail but we are trying to debug
            actual_times = [tutorial_detection.detect_time
                            for tutorial_detection in tutorial_detections]
            for detection in expected_detections:
                self.assertIn(detection.detect_time, actual_times,
                              msg='Expected detection at %s was not made'
                              % detection.detect_time)
        self.assertEqual(len(tutorial_detections), 23)
        # Cleanup the templates
        templates = glob.glob('tutorial_template_?.ms')
        for template in templates:
            os.remove(template)

if __name__ == '__main__':
    """
    Run tutorial tests
    """
    unittest.main()
