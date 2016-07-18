"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.core import subspace
import unittest


class SubspaceTestingMethods(unittest.TestCase):
    """
    Main tests for the subspace module.
    """
    def test_synthetic(self):
        """Test a synthetic case."""
        self.assertEqual('This passes', 'Not yet')

if __name__ == '__main__':
    unittest.main()