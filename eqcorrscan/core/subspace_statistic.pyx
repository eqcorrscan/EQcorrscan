# cython: linetrace=True
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
from scipy.linalg.blas import sdot, sgemv

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def det_statistic(float[:,:] detector,
                  float[:] data,
                  size_t inc):
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
    :type inc: int
    :param inc: Step size during stat loop

    :returns: Detection statistic from 0-1
    :rtype: np.ndarray
    """
    cdef size_t i, datamax = data.shape[0]
    cdef size_t ulen = detector.shape[0]
    cdef size_t stat_imax = (datamax // inc) - (ulen // inc) + 1
    cdef size_t dat_imax = (datamax - ulen) + 1
    cdef float[:] stats = np.zeros(dat_imax, dtype=DTYPE)
    cdef float[:,:] uut = np.dot(detector, detector.T)
    # Actual loop after static typing
    for i in range(0, dat_imax, inc):
        xp = np.dot(data[i:i+ulen].T, np.dot(uut, data[i:i+ulen]))
        xt = np.dot(data[i:i+ulen].T, data[i:i+ulen])
        stats[i] = (xp / xt)
    # Downsample stats
    stats = stats[::inc]
    # Cope with case of errored internal loop
    if np.all(np.isnan(stats)):
        return np.zeros(stat_imax, dtype=DTYPE)
    else:
        return np.array(stats)

