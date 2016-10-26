"""
Early development functions to do **very** basic simulations of seismograms \
to be used as general matched-filter templates and see how well a simple \
model would fit with real data.  Mostly used in EQcorrscan for testing.

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

import numpy as np
import warnings

from obspy import Stream, Trace, UTCDateTime

from eqcorrscan.utils import clustering


def seis_sim(sp, amp_ratio=1.5, flength=False, phaseout='all'):
    """
    Generate a simulated seismogram from a given S-P time.

    Will generate spikes separated by a given S-P time, which are then
    convolved with a decaying sine function.  The P-phase is simulated by a
    positive spike of value 1, the S-arrival is simulated by a decaying
    boxcar of maximum amplitude 1.5.  These amplitude ratios can be altered by
    changing the amp_ratio, which is the ratio S amplitude:P amplitude.

    .. note::
        In testing this can achieve 0.3 or greater cross-correlations with
        data.

    :type sp: int
    :param sp: S-P time in samples
    :type amp_ratio: float
    :param amp_ratio: S:P amplitude ratio
    :type flength: int
    :param flength: Fixed length in samples, defaults to False
    :type phaseout: str
    :param phaseout: Either 'P', 'S' or 'all', controls which phases to cut \
        around, defaults to 'all'. Can only be used with 'P' or 'S' options \
        if flength is set.

    :returns: Simulated data.
    :rtype: :class:`numpy.ndarray`
    """
    if flength and 2.5 * sp < flength and 100 < flength:
        additional_length = flength
    elif 2.5 * sp < 100.0:
        additional_length = 100
    else:
        additional_length = 2.5 * sp
    synth = np.zeros(int(sp + 10 + additional_length))
    # Make the array begin 10 samples before the P
    # and at least 2.5 times the S-P samples after the S arrival
    synth[10] = 1.0  # P-spike fixed at 10 samples from start of window
    # The length of the decaying S-phase should depend on the SP time,\
    # Some basic estimations suggest this should be atleast 10 samples\
    # and that the coda should be about 1/10 of the SP time
    S_length = 10 + int(sp // 3)
    S_spikes = np.arange(amp_ratio, 0, -(amp_ratio / S_length))
    # What we actually want, or what appears better is to have a series of\
    # individual spikes, of alternating polarity...
    for i in range(len(S_spikes)):
        if i in np.arange(1, len(S_spikes), 2):
            S_spikes[i] = 0
        if i in np.arange(2, len(S_spikes), 4):
            S_spikes[i] *= -1
    # Put these spikes into the synthetic
    synth[10 + sp:10 + sp + len(S_spikes)] = S_spikes
    # Generate a rough damped sine wave to convolve with the model spikes
    sine_x = np.arange(0, 10.0, 0.5)
    damped_sine = np.exp(-sine_x) * np.sin(2 * np.pi * sine_x)
    # Convolve the spike model with the damped sine!
    synth = np.convolve(synth, damped_sine)
    # Normalize snyth
    synth = synth / np.max(np.abs(synth))
    if not flength:
        return synth
    else:
        if phaseout in ['all', 'P']:
            synth = synth[0:flength]
        elif phaseout == 'S':
            synth = synth[sp:]
            if len(synth) < flength:
                # If this is too short, pad
                synth = np.append(synth, np.zeros(flength - len(synth)))
            else:
                synth = synth[0:flength]
        return synth


def SVD_sim(sp, lowcut, highcut, samp_rate,
            amp_range=np.arange(-10, 10, 0.01)):
    """
    Generate basis vectors of a set of simulated seismograms.

    Inputs should have a range of S-P amplitude ratios, in theory to simulate \
    a range of focal mechanisms.

    :type sp: int
    :param sp: S-P time in seconds - will be converted to samples according \
        to samp_rate.
    :type lowcut: float
    :param lowcut: Low-cut for bandpass filter in Hz
    :type highcut: float
    :param highcut: High-cut for bandpass filter in Hz
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz
    :type amp_range: numpy.ndarray
    :param amp_range: Amplitude ratio range to generate synthetics for.

    :returns: set of output basis vectors
    :rtype: :class:`numpy.ndarray`
    """
    # Convert SP to samples
    sp = int(sp * samp_rate)
    # Scan through a range of amplitude ratios
    synthetics = [Stream(Trace(seis_sim(sp, a))) for a in amp_range]
    for st in synthetics:
        for tr in st:
            tr.stats.station = 'SYNTH'
            tr.stats.channel = 'SH1'
            tr.stats.sampling_rate = samp_rate
            tr.filter('bandpass', freqmin=lowcut, freqmax=highcut)
    # We have a list of obspy Trace objects, we can pass this to EQcorrscan's
    # SVD functions
    V, s, U, stachans = clustering.SVD(synthetics)
    return V, s, U, stachans


def template_grid(stations, nodes, travel_times, phase, PS_ratio=1.68,
                  samp_rate=100, flength=False, phaseout='all'):
    """
    Generate a group of synthetic seismograms for a grid of sources.

    Used to simulate phase arrivals from a grid of known sources in a
    three-dimensional model.  Lags must be known and supplied, these can be
    generated from the bright_lights function: read_tt, and resampled to fit
    the desired grid dimensions and spacing using other functions therein.
    These synthetic seismograms are very simple models of seismograms using
    the seis_sim function herein. These approximate body-wave P and S first
    arrivals as spikes convolved with damped sine waves.

    :type stations: list
    :param stations: List of the station names
    :type nodes: list
    :param nodes: List of node locations in (lon,lat,depth)
    :type travel_times: numpy.ndarray
    :param travel_times: Array of travel times where travel_times[i][:] \
        refers to the travel times for station=stations[i], and \
        travel_times[i][j] refers to stations[i] for nodes[j]
    :type phase: str
    :param phase: Can be either 'P' or 'S'
    :type PS_ratio: float
    :param PS_ratio: P/S velocity ratio, defaults to 1.68
    :type samp_rate: float
    :param samp_rate: Desired sample rate in Hz, defaults to 100.0
    :type flength: int
    :param flength: Length of template in samples, defaults to False
    :type phaseout: str
    :param phaseout: Either 'S', 'P', 'all' or 'both', determines which \
        phases to clip around.  'all' Encompasses both phases in one channel, \
        but will return nothing if the flength is not long enough, 'both' \
        will return two channels for each stations, one SYN_Z with the \
        synthetic P-phase, and one SYN_H with the synthetic S-phase.

    :returns: List of :class:`obspy.core.stream.Stream`
    """
    if phase not in ['S', 'P']:
        raise IOError('Phase is neither P nor S')
    # Initialize empty list for templates
    templates = []
    # Loop through the nodes, for every node generate a template!
    for i, node in enumerate(nodes):
        st = []  # Empty list to be filled with synthetics
        # Loop through stations
        for j, station in enumerate(stations):
            tr = Trace()
            tr.stats.sampling_rate = samp_rate
            tr.stats.station = station
            tr.stats.channel = 'SYN'
            tt = travel_times[j][i]
            if phase == 'P':
                # If the input travel-time is the P-wave travel-time
                SP_time = (tt * PS_ratio) - tt
                if phaseout == 'S':
                    tr.stats.starttime += tt + SP_time
                else:
                    tr.stats.starttime += tt
            elif phase == 'S':
                # If the input travel-time is the S-wave travel-time
                SP_time = tt - (tt / PS_ratio)
                if phaseout == 'S':
                    tr.stats.starttime += tt
                else:
                    tr.stats.starttime += tt - SP_time
            else:
                raise IOError('Input grid is not P or S')
            # Set start-time of trace to be travel-time for P-wave
            # Check that the template length is long enough to include the SP
            if flength and SP_time * samp_rate < flength - 11 \
               and phaseout == 'all':
                tr.data = seis_sim(sp=int(SP_time * samp_rate), amp_ratio=1.5,
                                   flength=flength, phaseout=phaseout)
                st.append(tr)
            elif flength and phaseout == 'all':
                warnings.warn('Cannot make a bulk synthetic with this fixed ' +
                              'length for station ' + station)
            elif phaseout == 'all':
                tr.data = seis_sim(sp=int(SP_time * samp_rate), amp_ratio=1.5,
                                   flength=flength, phaseout=phaseout)
                st.append(tr)
            elif phaseout in ['P', 'S']:
                tr.data = seis_sim(sp=int(SP_time * samp_rate), amp_ratio=1.5,
                                   flength=flength, phaseout=phaseout)
                st.append(tr)
            elif phaseout == 'both':
                for _phaseout in ['P', 'S']:
                    _tr = tr.copy()
                    _tr.data = seis_sim(sp=int(SP_time * samp_rate),
                                        amp_ratio=1.5, flength=flength,
                                        phaseout=_phaseout)
                    if _phaseout == 'P':
                        _tr.stats.channel = 'SYN_Z'
                        # starttime defaults to S-time
                        _tr.stats.starttime = _tr.stats.starttime - SP_time
                    elif _phaseout == 'S':
                        _tr.stats.channel = 'SYN_H'
                    st.append(_tr)
        templates.append(Stream(st))
        # Stream(st).plot(size=(800,600))
    return templates


def generate_synth_data(nsta, ntemplates, nseeds, samp_rate, t_length,
                        max_amp, max_lag, debug=0):
    """
    Generate a synthetic dataset to be used for testing.

    This will generate both templates and data to scan through.
    Templates will be generated using the utils.synth_seis functions.
    The day of data will be random noise, with random signal-to-noise
    ratio copies of the templates randomly seeded throughout the day.
    It also returns the seed times and signal-to-noise ratios used.

    :type nsta: int
    :param nsta: Number of stations to generate data for < 15.
    :type ntemplates: int
    :param ntemplates: Number of templates to generate, will be generated \
        with random arrival times.
    :type nseeds: int
    :param nseeds: Number of copies of the template to seed within the \
        day of noisy data for each template.
    :type samp_rate: float
    :param samp_rate: Sampling rate to use in Hz
    :type t_length: float
    :param t_length: Length of templates in seconds.
    :type max_amp: float
    :param max_amp: Maximum signal-to-noise ratio of seeds.
    :param max_lag: Maximum lag time in seconds (randomised).
    :type max_lag: float
    :type debug: int
    :param debug: Debug level, bigger the number, the more plotting/output.

    :returns: Templates: List of :class:`obspy.core.stream.Stream`
    :rtype: list
    :returns: Data: :class:`obspy.core.stream.Stream` of seeded noisy data
    :rtype: :class:`obspy.core.stream.Stream`
    :returns: Seeds: dictionary of seed SNR and time with time in samples.
    :rtype: dict
    """
    # Generate random arrival times
    t_times = np.abs(np.random.random([nsta, ntemplates])) * max_lag
    # Generate random node locations - these do not matter as they are only
    # used for naming
    lats = np.random.random(ntemplates) * 90.0
    lons = np.random.random(ntemplates) * 90.0
    depths = np.abs(np.random.random(ntemplates) * 40.0)
    nodes = zip(lats, lons, depths)
    # Generating a 5x3 array to make 3 templates
    stations = ['ALPH', 'BETA', 'GAMM', 'KAPP', 'ZETA', 'BOB', 'MAGG',
                'ALF', 'WALR', 'ALBA', 'PENG', 'BANA', 'WIGG', 'SAUS',
                'MALC']
    if debug > 1:
        print(nodes)
        print(t_times)
        print(stations[0:nsta])
    templates = template_grid(stations=stations[0:nsta], nodes=nodes,
                              travel_times=t_times, phase='S',
                              samp_rate=samp_rate,
                              flength=int(t_length * samp_rate))
    if debug > 2:
        for template in templates:
            print(template)
            template.plot(size=(800, 600), equal_scale=False)
    # Now we want to create a day of synthetic data
    seeds = []
    data = templates[0].copy()  # Copy a template to get the correct length
    # and stats for data, we will overwrite the data on this copy
    for tr in data:
        tr.data = np.zeros(86400 * int(samp_rate))
        # Set all the traces to have a day of zeros
        tr.stats.starttime = UTCDateTime(0)
    for i, template in enumerate(templates):
        impulses = np.zeros(86400 * int(samp_rate))
        # Generate a series of impulses for seeding
        # Need three seperate impulse traces for each of the three templates,
        # all will be convolved within the data though.
        impulse_times = np.random.randint(86400 * int(samp_rate),
                                          size=nseeds)
        impulse_amplitudes = np.random.randn(nseeds) * max_amp
        # Generate amplitudes up to maximum amplitude in a normal distribution
        seeds.append({'SNR': impulse_amplitudes,
                      'time': impulse_times})
        for j in range(nseeds):
            impulses[impulse_times[j]] = impulse_amplitudes[j]
        # We now have one vector of impulses, we need nsta numbers of them,
        # shifted with the appropriate lags
        mintime = min([template_tr.stats.starttime
                       for template_tr in template])
        for j, template_tr in enumerate(template):
            offset = int((template_tr.stats.starttime - mintime) * samp_rate)
            pad = np.zeros(offset)
            tr_impulses = np.append(pad, impulses)[0:len(impulses)]
            # Convolve this with the template trace to give the daylong seeds
            data[j].data += np.convolve(tr_impulses,
                                        template_tr.data)[0:len(impulses)]
        if debug > 2:
            data.plot(starttime=UTCDateTime(0) +
                      impulse_times[0] / samp_rate - 3,
                      endtime=UTCDateTime(0) +
                      impulse_times[0] / samp_rate + 15)
    # Add the noise
    for tr in data:
        noise = np.random.randn(86400 * int(samp_rate))
        tr.data += noise / max(noise)
        if debug > 2:
            tr.plot()

    return templates, data, seeds


if __name__ == "__main__":
    import doctest
    doctest.testmod()
