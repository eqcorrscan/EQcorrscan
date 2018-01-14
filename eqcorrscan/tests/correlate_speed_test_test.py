"""
Testing of the command line script for correlation speed-testing.

This needs to test command-line calls, python calls and that it doesn't break
schtuff.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import shutil

from eqcorrscan.scripts.correlation_speeds import (
    generate_dataset, run_correlation, run_profiling)


if __name__ == '__main__':
    unittest.main()
