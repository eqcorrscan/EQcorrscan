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
cimport numpy as np
import cython
from cv2 import matchTemplate, TM_CCOEFF_NORMED

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def sliding_normxcorr(np.ndarray[DTYPE_t, ndim=1] arr1,
                      np.ndarray[DTYPE_t, ndim=1] arr2,
                      int shift_len,
                      int full_xcorr=False):
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

    .. Note:: This is not fast and is a bit silly - should be written in pure \
        c++ and wrapped.
    """
    cdef int i
    cdef int shift
    cdef DTYPE_t corr
    cdef int arr1_len = arr1.shape[0]
    cdef int arr2_len = arr2.shape[0]
    cdef int min_len
    cdef int arr1_start
    cdef int arr1_end
    cdef int arr2_start
    cdef int arr2_end
    cdef np.ndarray[DTYPE_t, ndim=1] xcorr = np.empty((shift_len * 2) + 1,
                                                      dtype=DTYPE)
    if arr1_len <= arr2_len:
        min_len = arr1_len
    else:
        min_len = arr2_len
    # De-mean
    arr1 -= arr1.mean()
    arr2 -= arr2.mean()
    # Normalize
    arr1 /= arr1.max()
    arr2 /= arr2.max()
    # Check that shift_len is less than half the shortest array length
    if shift_len * 2 > min_len:
        raise IndexError('Shift_len is too large, not enough support.')
    # Do the first half, reduces compute a little
    arr1_start = 0
    arr2_end = arr2_len
    for i in range(shift_len):
        arr1_end = arr1_len - (shift_len - i)
        arr2_start = arr2_end - (arr1_end - arr1_start)
        xcorr[i] = matchTemplate(arr1[arr1_start: arr1_end],
                                 arr2[arr2_start: arr2_end],
                                 method=TM_CCOEFF_NORMED)[0][0]
    arr1_end = arr1_len
    arr2_start = 0
    for i in range(shift_len + 1):
        arr2_end = arr2_len - i
        arr1_start = arr1_end - (arr2_end - arr2_start)
        xcorr[i + shift_len] = matchTemplate(arr1[arr1_start: arr1_end],
                                             arr2[arr2_start: arr2_end],
                                             method=TM_CCOEFF_NORMED)[0][0]

    shift = xcorr.argmax() - shift_len
    corr = xcorr.max()
    if full_xcorr == 1:
        return shift, corr, xcorr
    else:
        return shift, corr


