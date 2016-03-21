"""
Functions to test the seismo_logs functions for reading reftek logging \
information.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestLogReadMethods(unittest.TestCase):
    def test_timing(self):
        """Test to check whether timing information can be read correctly.
        """
        from eqcorrscan.utils import seismo_logs
        import datetime
        import os

        log_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                'test_data', 'LABE_testlog')
        startdate = datetime.date(2012, 10, 26)
        phase_errors = seismo_logs.rt_time_log(log_path, startdate)
        self.assertEqual(len(phase_errors), 24)

    def test_location(self):
        """Check that locations can be read in from the logfile.
        """
        from eqcorrscan.utils import seismo_logs
        import os

        log_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                'test_data', 'LABE_testlog')
        locations = seismo_logs.rt_location_log(log_path)
        self.assertEqual(len(locations), 24)


if __name__ == '__main__':
    unittest.main()
