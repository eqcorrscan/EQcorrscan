"""
Functions to stream-line normalisation of cross-correlation functions.

The specific intent of these functions is to reduce memory costs associated
with keeping large, redundant arrays in memory.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ctypes

import numpy as np

from eqcorrscan.utils.libnames import _load_cdll
from eqcorrscan.core.match_filter import MatchFilterError


# def normalise_corr(ccc, image, norm_sum, template_length):
#     """
#     Normalise a cross-correlation using bottleneck and scipy functions
#
#     :type ccc: np.ndarray
#     :param ccc: Non-normalised 2D array of cross-correlations
#     :type image: np.darray
#     :param image: Long data through which templates were scanned, 1D
#     :type norm_sum: np.ndarray
#     :param norm_sum: 1D array of normalised sums from templates
#     :type template_length: int
#     :param template_length: Number of samples in template
#
#     :return: np.ndarray, 2D same shape as ccc, but normalised.
#
#     .. Note::
#         Included for comparison only, bottleneck is not installed by default.
#         This routine is memory inefficient, but fast for small cases.
#     """
#     import bottleneck
#
#     image_mean_array = bottleneck.move_mean(
#         image, template_length)[template_length - 1:].astype(np.float32)
#     # CPU bound
#     image_std_array = bottleneck.move_std(
#         image, template_length)[template_length - 1:].astype(np.float32)
#
#     image_std_array[image_std_array == 0] = np.inf
#
#     normed = (ccc - norm_sum * image_mean_array) / image_std_array
#
#     return normed


def multi_norm(ccc, image, norm_sum, template_length):
    """
    Normalise a cross-correlation using memory efficient C routine.

    :type ccc: np.ndarray
    :param ccc: Non-normalised 2D array of cross-correlations
    :type image: np.darray
    :param image: Long data through which templates were scanned, 1D
    :type norm_sum: np.ndarray
    :param norm_sum: 1D array of normalised sums from templates
    :type template_length: int
    :param template_length: Number of samples in template

    :return: np.ndarray, 2D same shape as ccc, but normalised.
    """

    from future.utils import native_str

    norman = _load_cdll('norm')

    norman.multi_normalise.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int, ctypes.c_int]
    norman.multi_normalise.restype = ctypes.c_int

    ccc_len = ccc.shape[1]
    n = ccc.shape[0]
    image = np.ascontiguousarray(image, np.float32)
    norm_sum = np.ascontiguousarray(norm_sum.flatten(), np.float32)
    ccc = np.ascontiguousarray(ccc.T.flatten(), np.float32)
    ret = norman.multi_normalise(
        ccc, ccc_len, image, norm_sum, template_length, n
    )
    ccc[np.isnan(ccc)] = 0.0
    if ret != 0:
        raise MemoryError()
    if np.any(ccc > 1.001):
        raise MatchFilterError('Normalisation error in C code')
    ccc[ccc > 1.0] = 1.0
    ccc = ccc.reshape(ccc_len, n).T
    return np.round(ccc, 6)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
