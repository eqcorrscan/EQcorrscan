"""
Internal loop for subspace detection statistic calculation. Testing speedups.

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

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def det_statistic(np.ndarray[DTYPE_t, ndim=2] detector,
                  np.ndarray[DTYPE_t, ndim=1] data):
    """
    Base function to calculate the subspace detection statistic.

    Calculates for a given subspace detector and data stream. \
    The statistic is calculated by \
    projecting the data onto the N dimensional subspace defined by the given \
    detector following the equation: :math:'\\gamma = y^TUU^Ty' where y is \
    the data stream, U is the subspace detector and :math:'\\gamma' is the \
    detection statistic from 0 to 1.

    :type detector: np.ndarray
    :param detector: U matrix from singular value decomposition
    :type data: np.ndarry
    :param data: Data to detect within

    :returns: Detection statistic from 0-1
    :rtype: np.ndarray
    """
    cdef int i
    cdef int datamax = data.shape[0]
    cdef int ulen = detector.shape[1]
    cdef int umax = detector.shape[0]
    cdef int imax = datamax - ulen + 1
    cdef np.ndarray[DTYPE_t, ndim=1] stats = np.zeros(imax, dtype=DTYPE)
    # Check that there will not be an empty window
    cdef np.ndarray[DTYPE_t, ndim=1] _data = \
        np.concatenate([data, np.zeros((umax * imax) -
                                       datamax, dtype=DTYPE)])
    cdef np.ndarray[DTYPE_t, ndim=2] uut = np.dot(detector, detector.T)
    # Actual loop after static typing
    for i in range(imax):
        stats[i] = np.dot(_data[i:i + umax],
                          np.dot(uut, _data[i:i + umax].T))
    # Cope with case of errored internal loop
    if np.all(np.isnan(stats)):
        return np.zeros(imax, dtype=DTYPE)
    else:
        return stats

