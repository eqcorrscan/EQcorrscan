"""
Test function to do a **very** basic simulation of a seismogram to see how well
a simple model would fit with real data.

Written as part of the EQcorrscan project, but not yet distributed.

All copyright and ownership of this script belongs to Calum Chamberlain.

"""
import numpy as np
def seis_sim(SP, amp_ratio=1.5):
    """
    Function to generate a simulated seismogram from a given S-P time.
    Will generate spikes separated by a given S-P time, which are then convolved
    with a decaying sine function.  The P-phase is simulated by a positive
    spike of value 1, the S-arrival is simulated by a decaying boxcar of
    maximum amplitude 1.5.  These ampitude ratios can be altered by changing
    the amp_ratio, which is the ratio S amplitude:P amplitude.

    Note, in testing this can achieve 0.3 or greater cross-correlations with data

    :type SP: int
    :param SP: S-P time in samples
    :type amp_ratio: float
    :param amp_raio: S:P amplitude ratio

    :returns: np.ndarray
    """
    print amp_ratio
    synth=np.zeros(SP+10+100) # Make the array begin 10 samples before the P\
            # and at least 100 samples after the S arrival
    synth[10]=1.0
    # The length of the decaying S-phase should depend on the SP time,\
            # Some basic estimations suggest this should be atleast 10 samples\
            # and that the coda should be about 1/10 of the SP time
    S_length=int(10+SP/3.0)
    synth[10+SP:10+SP+S_length]=np.arange(amp_ratio,0, -(amp_ratio/S_length))
    # What we actually want, or what appears better is to have a series of\
            # individual spikes, of alternating polarity...
    for i in xrange(10+SP+1, 10+SP+S_length, 2):
        synth[i]=0
    for i in xrange(10+SP, 10+SP+S_length,3):
        synth[i]*=-1
    # Generate a rough damped sine wave to convolve with the model spikes
    sine_x=np.arange(0, 10.0, 0.5)
    damped_sine=np.exp(-sine_x) * np.sin(2 * np.pi * sine_x)
    # Convolve the spike model with the damped sine!
    synth=np.convolve(synth,damped_sine)
    # Normalize snyth
    synth=synth/np.max(synth)
    return synth

def SVD_sim(SP):
    """
    Function to generate a basis vectors of a set of simulated seismograms with
    a range of S-P amplitude ratios.

    :type SP: int
    :param SP: S-P time in samples

    :returns: nd.ndarray, set of output basis vectors
    """
    from obspy import Stream, Trace
    # Scan through a range of amplitude ratios
    synthetics=[Stream(Trace(seis_sim(SP, a))) for a in np.arange(-10, 10, 0.1)]
    # We have a list of obspy Trace objects, we can pass this to EQcorrscan's\
            # SVD functions
    from utils import clustering
    S,u,V = clustering.SVD(synthetics)
    return S, u, V
