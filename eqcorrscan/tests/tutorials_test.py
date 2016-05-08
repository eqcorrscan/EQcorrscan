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
        """Test the match_filter, must run after test_template_gen."""
        from eqcorrscan.tutorials.template_creation import mktemplates
        from eqcorrscan.tutorials.match_filter import run_tutorial
        import os
        import glob
        from obspy import read

        # Run mktemplates first to set-up for match_filter
        mktemplates(plot=False)
                    # network_code='GEONET',
                    # publicIDs=['2016p008122', '2016p008353'], plot=False)
        for template_no in range(2):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            self.assertTrue(len(template) > 1)
        del(template)
        # Run the matched-filter
        tutorial_detections = run_tutorial(plot=False)
        # It should make 19 detections in total...
        self.assertEqual(len(tutorial_detections), 20)
        # Cleanup the templates
        templates = glob.glob('tutorial_template_?.ms')
        for template in templates:
            os.remove(template)

if __name__ == '__main__':
    """
    Run tutorial tests
    """
    unittest.main()
