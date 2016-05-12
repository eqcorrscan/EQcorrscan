#!/usr/bin/python
r"""This script contains functions relevant to executing subspace detection \
for earthquake catalogs. The overarching function calls _channel_loop and \
_template_loop inner functions from match_filter in the same way as \
the normal matched filtering workflow.

We recommend that you read Harris' detailed report on subspace detection \
theory which can be found here: https://e-reports-ext.llnl.gov/pdf/335299.pdf

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
import warnings


def det_statistic(detector, data):
    r"""Base function to calculate the detection statistic for a given \
    subspace detector and data stream. The statistic is calculated by \
    projecting the data onto the N dimensional subspace defined by the given \
    detector following the equation: :math:'\\gamma = y^TUU^Ty' where y is \
    the data stream, U is the subspace detector and :math:'\\gamma' is the \
    detection statistic from 0 to 1.
    """
    day_stats = []
    for i in range(len(data) - len(detector[0]) + 1):
        y = data[i:i + len(detector[0])]
        day_stats.append(y.T.dot(detector.T).dot(detector).dot(y))
    #XXX TODO: Figure out why det_statistic not always giving between 0 and 1
    day_stats = np.asarray(day_stats)
    return day_stats


def subspace_detect(detector_names, detector_list, st, threshold,
                    threshold_type, trig_int, plotvar, plotdir='.', cores=1,
                    tempdir=False, debug=0, plot_format='jpg',
                    output_cat=False, extract_detections=False):
    r"""Overseer function to handle subspace detection. Modelled after \
    match_filter.match_filter().
    """


def plot_e_fraction():
