"""
Module written to generate lag-times and correlation values for the
cross-corrleation of a template event with a detected event.
Written by Calum Chamberlain of Victoria University, Wellington - 2015
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
