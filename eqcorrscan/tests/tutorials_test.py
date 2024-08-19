"""
Functions for testing the tutorials - written as a somewhat monolithic test
because we need certain behaviour.
"""
import unittest
import os
import pytest

from eqcorrscan.tutorials import subspace


class TestTutorialScripts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')

    @pytest.mark.superslow
    def test_subspace(self):
        """Test the subspace tutorial."""
        print("Running subspace")
        detections = subspace.run_tutorial(plot=False, cores=1, verbose=True)
        print("Subspace ran")
        if not len(detections) == 8:
            for detection in detections:
                print(detection)
        self.assertEqual(len(detections), 8)


if __name__ == '__main__':
    """
    Run tutorial tests
    """
    unittest.main()
