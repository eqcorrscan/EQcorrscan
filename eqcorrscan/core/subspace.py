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
    #XXX TODO: Complete function to loop though one day for one channel/temp

    #XXX TODO: Figure out why det_statistic not giving between 0 and 1


def plot_e_fraction():
