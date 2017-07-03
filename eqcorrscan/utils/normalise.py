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


# def normalise_corr(ccc, image, norm_sum, template_length):
#     """
#
#     :param ccc:
#     :param image:
#     :param norm_sum:
#     :param template_length:
#
#     :return:
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


def norm_compiled(ccc, image, norm_sum, template_length):
    """

    :param ccc:
    :param image:
    :param norm_sum:
    :param template_length:
    :return:
    """
    from future.utils import native_str

    norman = _load_cdll('norm')

    norman.normalise.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                               flags=native_str('C_CONTIGUOUS')),
        ctypes.c_float, ctypes.c_int]
    norman.normalise.restype = ctypes.c_int

    ccc_len = ccc.shape[1]
    image = np.ascontiguousarray(image, np.float32)

    for i in range(ccc.shape[0]):
        normalised = np.ascontiguousarray(ccc[i], np.float32)
        ret = norman.normalise(
            normalised, ccc_len, image, norm_sum[i], template_length)
        if ret:
            raise MemoryError()
        ccc[i] = normalised
    return ccc
