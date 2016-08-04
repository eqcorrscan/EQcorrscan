"""
Function for computing a sliding window normalised cross-correlation
function using openCV and numpy.  Written because obspy's xcorr gives different
results on different systems.  Designed only for small windows.

:copyright:
    Calum Chamberlain, Chet Hopp.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
from eqcorrscan.core import sliding_xcorrlib


def sliding_xcorr(arr1, arr2, shift_len, full_xcorr=False):
    """
    Calculate the sliding normalized cross-correlation of a pair of arrays.

    Uses the openCV match_template routine for correlation, and implements a \
    cythonised sliding-window approach. shift_len must be less than half the \
    shortest array length to provide reasonable support. Shifts +/- the \
    shift_len

    :type arr1: np.ndarray
    :param arr1: First array, shifts will be relative to this array.  Must \
        be np.float32 datatype
    :type arr2: np.ndarry
    :param arr2: Second array to slide through arr 1.  Must be np.float32 \
        datatype.
    :type shift_len: int
    :param shift_len: Length in samples to compute the sliding window for.
    :type full_xcorr: int
    :param full_xcorr: Whether to return the full cross-correlation or not. \
        If False, only the shift for the maximum correlation and the \
        correlation value will be returned, if True, will also return the \
        correlation vector. 0=False, 1=True

    :rtype: float
    :return: Shift (in samples) for maximum correlation
    :type: float
    :return: Maximum correlation.
    :rtype: np.ndarray
    :return: Full cross-correlation array if full_xcorr=True.
    """
    # De-mean
    arr1 -= arr1.mean()
    arr2 -= arr2.mean()
    # Normalize
    arr1 /= arr1.max()
    arr2 /= arr2.max()
    # Make contiguous
    arr1 = np.ascontiguousarray(arr1, np.float32)
    arr2 = np.ascontiguousarray(arr2, np.float32)
    # Check that shift_len is less than half the shortest array length
    if shift_len * 2 > np.min([len(arr1), len(arr2)]):
        raise IndexError('Shift_len is too large, not enough support.')
    xcorr = np.empty((shift_len * 2) + 1, dtype=np.float32, order='C')
    print(len(arr1))
    print(len(arr2))
    res = sliding_xcorrlib.xcorr(arr1, arr2, xcorr, shift_len, len(arr1),
                                 len(arr2))
    if res:
        raise MemoryError('Internal C error %s' % res)
    shift = xcorr.argmax() - shift_len
    corr = xcorr.max()
    if full_xcorr:
        return shift, corr, xcorr
    else:
        return shift, corr
