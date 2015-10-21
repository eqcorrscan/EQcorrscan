#!/usr/bin/python
"""
Script to test if all dependencies are installed and running for the EQcorrscan
package.
"""

def test_import():
    import sys
    sys.path.insert(0,'/usr/lib/pyshared/python2.7') # Insert path for travis
    i=0
    try:
        import cv2
    except:
        print("You have not properly installed: cv2")
        i+=1
    try:
        import joblib
    except:
        print("You have not properly installed: joblib")
        i+=1
    try:
        import numpy
    except:
        print("You have not properly installed: numpy")
        i+=1
    try:
        import matplotlib.pyplot
    except:
        print("You have not properly installed: matplotlib")
        i+=1
    try:
        import scipy
    except:
        print("You have not properly installed: scipy")
        i+=1
    try:
        from obspy import read
    except:
        print("You have not properly installed: obspy")
        i+=1
    if not i == 0:
        return False
    return True


def run():
    """
    Where we call all the available tests from
    """
    from eqcorrscan.utils import Sfile_util
    assert test_import() == True
    assert Sfile_util.test_rw() == True

if __name__ == '__main__':
    run()
