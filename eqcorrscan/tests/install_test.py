"""Script to test if all dependencies are installed and running for the \
EQcorrscan package.
"""
import unittest


class TestImport(unittest.TestCase):
    def test_import(self):
        import sys
        sys.path.insert(0, '/usr/lib/pyshared/python2.7')
        # Insert path for travis
        i = 0
        try:
            import cv2  # NOQA
        except:
            print("You have not properly installed: cv2")
            i += 1
        try:
            import joblib  # NOQA
        except:
            print("You have not properly installed: joblib")
            i += 1
        try:
            import numpy  # NOQA
        except:
            print("You have not properly installed: numpy")
            i += 1
        try:
            import matplotlib.pyplot  # NOQA
        except:
            print("You have not properly installed: matplotlib")
            i += 1
        try:
            import scipy  # NOQA
        except:
            print("You have not properly installed: scipy")
            i += 1
        try:
            from obspy import read  # NOQA
        except:
            print("You have not properly installed: obspy")
            i += 1
        self.assertEqual(i, 0)

if __name__ == '__main__':
    unittest.main()
