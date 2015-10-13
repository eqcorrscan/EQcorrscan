"""
Module written to generate lag-times and correlation values for the
cross-corrleation of a template event with a detected event.
Written by Calum Chamberlain of Victoria University, Wellington - 2015

NOTE - unfinished!

Copyright 2015 Calum Chamberlain

This file is part of EQcorrscan.

    EQcorrscan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EQcorrscan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EQcorrscan.  If not, see <http://www.gnu.org/licenses/>.


"""
import numpy as np

def laggen(template, image):
    """
    Function to generate the lag and correlation between a template
    trace and an image trace based on their normalised cross-correlation value.

    :type template: obspy.Trace
    :param template: an obspy Trace object containing the waveform of the template
    :type image: obspy.Trace
    :param image: an obspy Trace object containing the waveform of the image

    :returns: ccc

    .. rubric:: Note
    Image can not be shorter than template, but may be longer.
    """
    from core.match_filter import normxcorr2
    # Compute the correlation between the image and the template
    ccc=normxcorr2(template.data, image.data)
    lag=template.stats.sampling_rate*np.argmax(ccc)
    correlation=np.max(ccc)
    return lag, correlation
