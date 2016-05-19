"""
Functions to test generating templates from SAC data.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestTemplateGeneration(unittest.TestCase):
    """Test the reading a writing of pick info."""
    def test_sac_template_gen(self):
        """Test template generation."""
        from eqcorrscan.core.template_gen import from_sac
        import glob
        import os
        import obspy

        samp_rate = 20
        length = 8

        for event in ['2014p611252', 'No_head']:
            test_files = os.path.join(os.path.abspath(os.path.
                                                      dirname(__file__)),
                                      'test_data', 'SAC', event, '*')
            sac_files = glob.glob(test_files)

            # We currently do not support SAC template generation below version
            # 1.0.0 as before this, SACIO did not fill the reference time,
            # which is needed for defining pick times.  This is usually the
            # trace start time, but isn't always...
            if int(obspy.__version__.split('.')[0]) >= 1:
                template = from_sac(sac_files, lowcut=2.0, highcut=8.0,
                                    samp_rate=samp_rate, filt_order=4,
                                    length=length,
                                    swin='all', prepick=0.1, debug=0,
                                    plot=False)
                for tr in template:
                    self.assertEqual(len(tr.data), length * samp_rate)
            else:
                with self.assertRaises(NotImplementedError):
                    template = from_sac(sac_files, lowcut=2.0, highcut=8.0,
                                        samp_rate=samp_rate, filt_order=4,
                                        length=length,
                                        swin='all', prepick=0.1, debug=0,
                                        plot=False)

    def test_tutorial_template_gen(self):
        """Test template generation from tutorial, uses from_client method.

        Checks that the tutorial generates the templates we expect it to!
        """
        from obspy import read
        from eqcorrscan.tutorials.template_creation import mktemplates
        import os
        import numpy as np

        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        mktemplates(plot=False)
        for template_no in range(4):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            expected_template = read(os.path.join(testing_path,
                                                  'tutorial_template_' +
                                                  str(template_no) + '.ms'))
            for tr in template:
                expected_tr = expected_template.select(station=tr.stats.
                                                       station,
                                                       channel=tr.stats.
                                                       channel)[0]
                self.assertTrue((expected_tr.data.astype(np.float32) ==
                                 tr.data.astype(np.float32)).all())
            del(template)
            os.remove('tutorial_template_' + str(template_no) + '.ms')

if __name__ == '__main__':
    unittest.main()
