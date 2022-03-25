"""
Functions to aid magnitude estimation.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import logging
import eqcorrscan  # Used to get version number
import os
import glob
import matplotlib.pyplot as plt
import itertools
import copy
import random
import pickle
import math

from inspect import currentframe
from scipy.signal import iirfilter, sosfreqz
from collections import Counter
from obspy import Stream, Trace
from obspy.signal.invsim import simulate_seismometer as seis_sim
from obspy.core.event import (
    Amplitude, Pick, WaveformStreamID, Origin, ResourceIdentifier)
from obspy.geodetics import degrees2kilometers


Logger = logging.getLogger(__name__)


# Magnitude - frequency funcs

def calc_max_curv(magnitudes, bin_size=0.5, plotvar=False):
    """
    Calculate the magnitude of completeness using the maximum curvature method.

    :type magnitudes: list or numpy array
    :param magnitudes:
        List of magnitudes from which to compute the maximum curvature which
        will give an estimate of the magnitude of completeness given the
        assumption of a power-law scaling.
    :type bin_size: float
    :param bin_size:
        Width of magnitude bins used to compute the non-cumulative distribution
    :type plotvar: bool
    :param plotvar: Turn plotting on and off

    :rtype: float
    :return: Magnitude at maximum curvature

    .. Note:: Should be used as a guide, often under-estimates Mc.

    .. rubric:: Example

    >>> import numpy as np
    >>> mags = np.arange(3, 6, .1)
    >>> N = 10 ** (5 - 1 * mags)
    >>> magnitudes = [0, 2, 3, 2.5, 2.2, 1.0]  # Some below completeness
    >>> for mag, n in zip(mags, N):
    ...     magnitudes.extend([mag for _ in range(int(n))])
    >>> calc_max_curv(magnitudes, plotvar=False)
    3.0
    """
    min_bin, max_bin = int(min(magnitudes)), int(max(magnitudes) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(magnitudes, bins)
    grad = (df[1:] - df[0:-1]) / bin_size
    # Need to find the second order derivative
    curvature = (grad[1:] - grad[0:-1]) / bin_size
    max_curv = bins[np.argmax(np.abs(curvature))] + bin_size
    if plotvar:
        fig, ax = plt.subplots()
        ax.scatter(bins[:-1] + bin_size / 2, df, color="k",
                   label="Magnitudes")
        ax.axvline(x=max_curv, color="red", label="Maximum curvature")
        ax1 = ax.twinx()
        ax1.plot(bins[:-1] + bin_size / 2, np.cumsum(df[::-1])[::-1],
                 color="k", label="Cumulative distribution")
        ax1.scatter(bins[1:-1], grad, color="r", label="Gradient")
        ax2 = ax.twinx()
        ax2.scatter(bins[1:-2] + bin_size, curvature, color="blue",
                    label="Curvature")
        # Code borrowed from https://matplotlib.org/3.1.1/gallery/ticks_and_
        # spines/multiple_yaxis_with_spines.html#sphx-glr-gallery-ticks-and-
        # spines-multiple-yaxis-with-spines-py
        ax2.spines["right"].set_position(("axes", 1.2))
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        for sp in ax2.spines.values():
            sp.set_visible(False)
        ax2.spines["right"].set_visible(True)

        ax.set_ylabel("N earthquakes in bin")
        ax.set_xlabel("Magnitude")
        ax1.set_ylabel("Cumulative events and gradient")
        ax2.set_ylabel("Curvature")
        fig.legend()
        fig.show()
    return float(max_curv)


def calc_b_value(magnitudes, completeness, max_mag=None, plotvar=True):
    """
    Calculate the b-value for a range of completeness magnitudes.

    Calculates a power-law fit to given magnitudes for each completeness
    magnitude.  Plots the b-values and residuals for the fitted catalogue
    against the completeness values. Computes fits using numpy.polyfit,
    which uses a least-squares technique.

    :type magnitudes: list
    :param magnitudes: Magnitudes to compute the b-value for.
    :type completeness: list
    :param completeness: list of completeness values to compute b-values for.
    :type max_mag: float
    :param max_mag: Maximum magnitude to attempt to fit in magnitudes.
    :type plotvar: bool
    :param plotvar: Turn plotting on or off.

    :rtype: list
    :return:
        List of tuples of (completeness, b-value, residual, number of
        magnitudes used)

    .. Note::
        High "residuals" indicate better fit. Residuals are calculated
        according to the Wiemer & Wyss 2000, Minimum Magnitude of Completeness
        in Earthquake Catalogs: Examples from Alaska, the Western United
        States, and Japan, BSSA.

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.mag_calc import calc_b_value
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime('2012-03-26T00:00:00')
    >>> t2 = t1 + (3 * 86400)
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
    >>> magnitudes = [event.magnitudes[0].mag for event in catalog]
    >>> b_values = calc_b_value(magnitudes, completeness=np.arange(3, 7, 0.2),
    ...                         plotvar=False)
    >>> round(b_values[4][1], 1)
    1.0
    >>> # We can set a maximum magnitude:
    >>> b_values = calc_b_value(magnitudes, completeness=np.arange(3, 7, 0.2),
    ...                         plotvar=False, max_mag=5)
    >>> round(b_values[4][1], 1)
    1.0
    """
    b_values = []
    # Calculate the cdf for all magnitudes
    counts = Counter(magnitudes)
    cdf = np.zeros(len(counts))
    mag_steps = np.zeros(len(counts))
    for i, magnitude in enumerate(sorted(counts.keys(), reverse=True)):
        mag_steps[i] = magnitude
        if i > 0:
            cdf[i] = cdf[i - 1] + counts[magnitude]
        else:
            cdf[i] = counts[magnitude]

    if not max_mag:
        max_mag = max(magnitudes)
    for m_c in completeness:
        if m_c >= max_mag or m_c >= max(magnitudes):
            Logger.warning('Not computing completeness at %s, above max_mag' %
                           str(m_c))
            break
        complete_mags = []
        complete_freq = []
        for i, mag in enumerate(mag_steps):
            if mag >= m_c <= max_mag:
                complete_mags.append(mag)
                complete_freq.append(np.log10(cdf[i]))
        if len(complete_mags) < 4:
            Logger.warning('Not computing completeness above ' + str(m_c) +
                           ', fewer than 4 events')
            break
        fit = np.polyfit(complete_mags, complete_freq, 1, full=True)
        # Calculate the residuals according to the Wiemer & Wys 2000 definition
        predicted_freqs = [fit[0][1] - abs(fit[0][0] * M)
                           for M in complete_mags]
        r = 100 - ((np.sum([abs(complete_freq[i] - predicted_freqs[i])
                           for i in range(len(complete_freq))]) * 100) /
                   np.sum(complete_freq))
        b_values.append((m_c, abs(fit[0][0]), r, len(complete_mags)))
    if plotvar:
        fig, ax1 = plt.subplots()
        b_vals = ax1.scatter(list(zip(*b_values))[0], list(zip(*b_values))[1],
                             c='k')
        resid = ax1.scatter(list(zip(*b_values))[0],
                            [100 - b for b in list(zip(*b_values))[2]], c='r')
        ax1.set_ylabel('b-value and residual')
        plt.xlabel('Completeness magnitude')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Number of events used in fit')
        n_ev = ax2.scatter(list(zip(*b_values))[0], list(zip(*b_values))[3],
                           c='g')
        fig.legend((b_vals, resid, n_ev),
                   ('b-values', 'residuals', 'number of events'),
                   'lower right')
        ax1.set_title('Possible completeness values')
        plt.show()
    return b_values


# Helpers for local magnitude estimation
# Note Wood anderson sensitivity is 2080 as per Uhrhammer & Collins 1990
PAZ_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
          'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}


def dist_calc(loc1, loc2):
    """
    Function to calculate the distance in km between two points.

    Uses the
    `haversine formula <https://en.wikipedia.org/wiki/Haversine_formula>`_
    to calculate great circle distance at the Earth's surface, then uses
    trig to include depth.

    :type loc1: tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: tuple
    :param loc2: Tuple of lat, lon, depth (in decimal degrees and km)

    :returns: Distance between points in km.
    :rtype: float
    """
    from eqcorrscan.utils.libnames import _load_cdll
    import ctypes

    utilslib = _load_cdll('libutils')

    utilslib.dist_calc.argtypes = [
        ctypes.c_float, ctypes.c_float, ctypes.c_float,
        ctypes.c_float, ctypes.c_float, ctypes.c_float]
    utilslib.dist_calc.restype = ctypes.c_float

    dist = utilslib.dist_calc(
        float(math.radians(loc1[0])), float(math.radians(loc1[1])),
        float(loc1[2]),
        float(math.radians(loc2[0])), float(math.radians(loc2[1])),
        float(loc2[2]))
    return dist


def _sim_WA(trace, inventory, water_level, velocity=False):
    """
    Remove the instrument response from a trace and simulate a Wood-Anderson.

    Returns a de-meaned, de-trended, Wood Anderson simulated trace in
    its place.

    Works in-place on data and will destroy your original data, copy the
    trace before giving it to this function!

    :type trace: obspy.core.trace.Trace
    :param trace:
        A standard obspy trace, generally should be given without
        pre-filtering, if given with pre-filtering for use with
        amplitude determination for magnitudes you will need to
        worry about how you cope with the response of this filter
        yourself.
    :type inventory: obspy.core.inventory.Inventory
    :param inventory:
        Inventory containing response information for the stations in st.
    :type water_level: float
    :param water_level: Water level for the simulation.
    :type velocity: bool
    :param velocity:
        Whether to return a velocity trace or not - velocity is non-standard
        for Wood-Anderson instruments, but institutes that use seiscomp3 or
        Antelope require picks in velocity.

    :returns: Trace of Wood-Anderson simulated data
    :rtype: :class:`obspy.core.trace.Trace`
    """
    assert isinstance(trace, Trace)
    paz_wa = copy.deepcopy(PAZ_WA)
    # Need to make a copy because we might edit it.
    if velocity:
        paz_wa['zeros'] = [0 + 0j, 0 + 0j]
    # De-trend data
    trace.detrend('simple')
    # Remove response to Velocity
    try:
        trace.remove_response(
            inventory=inventory, output="VEL", water_level=water_level)
    except Exception:
        Logger.error(f"No response for {trace.id} at {trace.stats.starttime}")
        return None
    # Simulate Wood Anderson
    trace.data = seis_sim(trace.data, trace.stats.sampling_rate,
                          paz_remove=None, paz_simulate=paz_wa,
                          water_level=water_level)
    return trace


def _max_p2t(data, delta, return_peak_trough=False):
    """
    Finds the maximum peak-to-trough amplitude and period.

    Originally designed to be used to calculate magnitudes (by
    taking half of the peak-to-trough amplitude as the peak amplitude).

    :type data: numpy.ndarray
    :param data: waveform trace to find the peak-to-trough in.
    :type delta: float
    :param delta: Sampling interval in seconds
    :type return_peak_trough: bool
    :param return_peak_trough:
        Optionally return the peak and trough

    :returns:
        tuple of (amplitude, period, time) with amplitude in the same
        scale as given in the input data, and period in seconds, and time in
        seconds from the start of the data window.
    :rtype: tuple
    """
    turning_points = []  # A list of tuples of (amplitude, sample)
    for i in range(1, len(data) - 1):
        if (data[i] < data[i - 1] and data[i] < data[i + 1]) or\
           (data[i] > data[i - 1] and data[i] > data[i + 1]):
            turning_points.append((data[i], i))
    if len(turning_points) >= 1:
        amplitudes = np.empty([len(turning_points) - 1],)
        half_periods = np.empty([len(turning_points) - 1],)
    else:
        Logger.warning(
            'Turning points has length: ' + str(len(turning_points)) +
            ' data have length: ' + str(len(data)))
        return 0.0, 0.0, 0.0
    for i in range(1, len(turning_points)):
        half_periods[i - 1] = (delta * (turning_points[i][1] -
                                        turning_points[i - 1][1]))
        amplitudes[i - 1] = np.abs(turning_points[i][0] -
                                   turning_points[i - 1][0])
    amplitude = np.max(amplitudes)
    period = 2 * half_periods[np.argmax(amplitudes)]
    delay = delta * turning_points[np.argmax(amplitudes)][1]
    if not return_peak_trough:
        return amplitude, period, delay
    max_position = np.argmax(amplitudes)
    peak = max(
        t[0] for t in turning_points[max_position: max_position + 2])
    trough = min(
        t[0] for t in turning_points[max_position: max_position + 2])
    return amplitude, period, delay, peak, trough


def _pairwise(iterable):
    """
    Wrapper on itertools for SVD_magnitude.
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


# Helpers for relative magnitude calculation

def _get_pick_for_station(event, station, channel, use_s_picks):
    """
    Get the first reported pick for a given station.

    :type event: `obspy.core.event.Event`
    :param event: Event with at least picks
    :type station: str
    :param station: Station to get pick for
    :type use_s_picks: bool
    :param use_s_picks: Whether to allow S-picks to be returned

    :rtype: `obspy.core.event.Pick`
    :return: First reported pick for station
    """
    picks = [p for p in event.picks if p.waveform_id.station_code == station
             and p.waveform_id.channel_code == channel]
    if len(picks) == 0:
        Logger.info("No pick for {0}".format(station))
        return None
    picks.sort(key=lambda p: p.time)
    for pick in picks:
        if pick.phase_hint and pick.phase_hint[0].upper() == 'S'\
                and not use_s_picks:
            continue
        return pick
    Logger.info("No suitable pick found for {0}".format(station))
    return None


def _snr(tr, noise_window, signal_window):
    """
    Compute ratio of maximum signal amplitude to rms noise amplitude.

    :param tr: Trace to compute signal-to-noise ratio for
    :param noise_window: (start, end) of window to use for noise
    :param signal_window: (start, end) of window to use for signal

    :return: Signal-to-noise ratio, noise amplitude
    """
    from eqcorrscan.core.template_gen import _rms

    noise_amp = _rms(
        tr.slice(starttime=noise_window[0], endtime=noise_window[1]).data)
    if np.isnan(noise_amp):
        Logger.warning("Could not calculate noise with this data, setting "
                       "to 1")
        noise_amp = 1.0
    try:
        signal_amp = tr.slice(
            starttime=signal_window[0], endtime=signal_window[1]).data.max()
    except ValueError as e:
        Logger.error(e)
        return np.nan
    return signal_amp / noise_amp


def _get_signal_and_noise(stream, event, seed_id, noise_window,
                          signal_window, use_s_picks):
    """
    Get noise and signal RMS-amplitudes and signal standard deviation for an
    event on a specific channel.

    (Until v.0.4.3, this function calculated noise amplitude as the RMS
    amplitude of the noise window and signal amplitude as the maximum amplitude
    in the signal window. This was changed to only RMS amplitudes to align it
    with the methodology in Schaff & Richards 2014-paper.)
    """
    from eqcorrscan.core.template_gen import _rms

    station = seed_id.split('.')[1]
    channel = seed_id.split('.')[3]
    pick = _get_pick_for_station(
        event=event, station=station, channel=channel, use_s_picks=use_s_picks)
    if pick is None:
        Logger.error("No pick for {0}".format(station))
        return None, None, None
    tr = stream.select(id=seed_id).merge()
    if len(tr) == 0:
        return None, None, None
    tr = tr[0]
    noise = tr.slice(
        starttime=pick.time + noise_window[0],
        endtime=pick.time + noise_window[1]).data
    noise_amp = _rms(noise)
    if np.isnan(noise_amp):
        noise_amp = None
    signal = tr.slice(
        starttime=pick.time + signal_window[0],
        endtime=pick.time + signal_window[1]).data
    if len(signal) == 0:
        Logger.debug("No signal data between {0} and {1}".format(
            pick.time + signal_window[0], pick.time + signal_window[1]))
        Logger.debug(tr)
        return noise_amp, None, None
    signal_amp = _rms(signal)
    return noise_amp, signal_amp, signal.std()


def relative_amplitude(st1, st2, event1, event2, noise_window=(-20, -1),
                       signal_window=(-.5, 20), min_snr=1.5,
                       use_s_picks=False):
    """
    Compute the relative amplitudes between two streams.

    Uses standard deviation of amplitudes within trace. Relative amplitudes are
    computed as:

    .. math::

       \\frac{std(tr2)}{std(tr1)}

    where tr1 is a trace from st1 and tr2 is a matching (seed ids match) trace
    from st2.  The standard deviation of the amplitudes is computed in the
    signal window given. If the ratio of amplitudes between the signal window
    and the noise window is below `min_snr` then no result is returned for that
    trace. The SNR here is defined as the ratio of RMS-amplitudes of signal
    and noise (equal to ratio of L2-norms of signal and noise, but normalized
    for signal length). The Windows are computed relative to the first pick
    for that station.

    If one stream has insufficient data to estimate noise amplitude, the noise
    amplitude of the other will be used.

    :type st1: `obspy.core.stream.Stream`
    :param st1: Stream for event1
    :type st2: `obspy.core.stream.Stream`
    :param st2: Stream for event2
    :type event1: `obspy.core.event.Event`
    :param event1: Event with picks (nothing else is needed)
    :type event2: `obspy.core.event.Event`
    :param event2: Event with picks (nothing else is needed)
    :type noise_window: tuple of float
    :param noise_window:
        Start and end of noise window in seconds relative to pick
    :type signal_window: tuple of float
    :param signal_window:
        Start and end of signal window in seconds relative to pick
    :type min_snr: float
    :param min_snr: Minimum signal-to-noise ratio allowed to make a measurement
    :type use_s_picks: bool
    :param use_s_picks:
        Whether to allow relative amplitude estimates to be made from S-picks.
        Note that noise and signal windows are relative to pick-times, so using
        an S-pick might result in a noise window including P-energy.

    :rtype: dict, dict, dict
    :return:
        Dictionary of relative amplitudes keyed by seed-id
        Dictionary of signal-to-noise ratios for st1
        Dictionary of signal-to-noise ratios for st2
    """
    # keep input safe
    event1 = event1.copy()
    # sort out S-picks if not to be used
    if not use_s_picks:
        event1.picks = [p for p in event1.picks if p.phase_hint[0] != "S"]
        st1 = Stream(
            [tr for tr in st1.copy() if (tr.stats.station, tr.stats.channel) in
             [(p.waveform_id.station_code, p.waveform_id.channel_code)
              for p in event1.picks]])
    seed_ids = {tr.id for tr in st1}.intersection({tr.id for tr in st2})
    amplitudes = {}
    snrs_1 = {}
    snrs_2 = {}
    for seed_id in seed_ids:
        noise1, signal1, std1 = _get_signal_and_noise(
            stream=st1, event=event1, signal_window=signal_window,
            noise_window=noise_window, use_s_picks=use_s_picks,
            seed_id=seed_id)
        noise2, signal2, std2 = _get_signal_and_noise(
            stream=st2, event=event2, signal_window=signal_window,
            noise_window=noise_window, use_s_picks=use_s_picks,
            seed_id=seed_id)
        noise1 = noise1 or noise2
        noise2 = noise2 or noise1
        if noise1 is None or noise2 is None:
            Logger.info("Insufficient data for noise to be estimated for "
                        "{0}".format(seed_id))
            continue
        if signal1 is None or signal2 is None:
            Logger.info("No signal data found for {0}".format(seed_id))
            continue
        snr1 = np.nan_to_num(signal1 / noise1)
        snr2 = np.nan_to_num(signal2 / noise2)
        if snr1 < min_snr or snr2 < min_snr:
            Logger.info("SNR (event1: {0:.2f}, event2: {1:.2f} too low "
                        "for {2}".format(snr1, snr2, seed_id))
            continue
        ratio = std2 / std1
        Logger.debug("Channel: {0} Relative amplitude: {1:.2f}".format(
            seed_id, ratio))
        amplitudes.update({seed_id: ratio})
        snrs_1.update({seed_id: snr1})
        snrs_2.update({seed_id: snr2})
    return amplitudes, snrs_1, snrs_2


# Magnitude estimation functions

def relative_magnitude(st1, st2, event1, event2, noise_window=(-20, -1),
                       signal_window=(-.5, 20), min_snr=5.0, min_cc=0.7,
                       use_s_picks=False, correlations=None, shift=.2,
                       return_correlations=False, correct_mag_bias=True):
    """
    Compute the relative magnitudes between two events.

    See :func:`eqcorrscan.utils.mag_calc.relative_amplitude` for information
    on how relative amplitudes are calculated. To compute relative magnitudes
    from relative amplitudes this function can weight the amplitude ratios by
    the cross-correlation of the two events. The relation used is similar to
    Schaff and Richards (2014), equation 4 and is:

    .. math::

        \\Delta m = \\log{\\frac{std(tr2)}{std(tr1)}} + \\log{
            \\frac{(1+\\frac{1}{snr_x^2})}{1+\\frac{1}{snr_y^2}}\\times CC}

    If you decide to use this function you should definitely read the paper
    to understand what you can use this for and cite the paper!

    :type st1: `obspy.core.stream.Stream`
    :param st1: Stream for event1
    :type st2: `obspy.core.stream.Stream`
    :param st2: Stream for event2
    :type event1: `obspy.core.event.Event`
    :param event1: Event with picks (nothing else is needed)
    :type event2: `obspy.core.event.Event`
    :param event2: Event with picks (nothing else is needed)
    :type noise_window: tuple of float
    :param noise_window:
        Start and end of noise window in seconds relative to pick
    :type signal_window: tuple of float
    :param signal_window:
        Start and end of signal window in seconds relative to pick
    :type min_snr: float
    :param min_snr: Minimum signal-to-noise ratio allowed to make a measurement
    :type min_cc: float
    :param min_cc:
        Minimum inter-event correlation (between -1 and 1) allowed to make a
        measurement.
    :type use_s_picks: bool
    :param use_s_picks:
        Whether to allow relative amplitude estimates to be made from S-picks.
        Note that noise and signal windows are relative to pick-times, so using
        an S-pick might result in a noise window including P-energy.
    :type correlations: dict
    :param correlations:
        Pre-computed dictionary of correlations keyed by seed-id. If None
        (default) then correlations will be computed for the provided data in
        the `signal_window`.
    :type shift: float
    :param shift:
        Shift length for correlations in seconds - maximum correlation within
        a window between +/- shift of the P-pick will be used to weight the
        magnitude.
    :type return_correlations: bool
    :param return_correlations:
        If true will also return maximum correlations as a dictionary.
    :type correct_mag_bias: bool
    :param correct_mag_bias:
        Whether to correct for the magnitude-bias introduced by cc<1 and the
        presence of noise (i.e., SNR << ∞). Without bias-correction, the
        relative magnitudes are simple L2-norm-ratio relative magnitudes.

    :rtype: dict
    :return: Dictionary of relative magnitudes keyed by seed-id
    """
    import math
    from obspy.signal.cross_correlation import correlate

    relative_magnitudes = {}
    compute_correlations = False
    if correlations is None:
        correlations = {}
        compute_correlations = True
    relative_amplitudes, snrs_1, snrs_2 = relative_amplitude(
        st1=st1, st2=st2, event1=event1, event2=event2,
        noise_window=noise_window, signal_window=signal_window,
        min_snr=min_snr, use_s_picks=use_s_picks)
    for seed_id, amplitude_ratio in relative_amplitudes.items():
        tr1 = st1.select(id=seed_id)[0]
        tr2 = st2.select(id=seed_id)[0]
        pick1 = _get_pick_for_station(
            event=event1, station=tr1.stats.station, channel=tr1.stats.channel,
            use_s_picks=use_s_picks)
        pick2 = _get_pick_for_station(
            event=event2, station=tr2.stats.station, channel=tr2.stats.channel,
            use_s_picks=use_s_picks)
        if compute_correlations:
            cc = correlate(
                tr1.slice(
                    starttime=pick1.time + signal_window[0],
                    endtime=pick1.time + signal_window[1]),
                tr2.slice(
                    starttime=pick2.time + signal_window[0],
                    endtime=pick2.time + signal_window[1]),
                shift=int(shift * tr1.stats.sampling_rate))
            cc = cc.max()
            correlations.update({seed_id: cc})
        else:
            cc = correlations.get(seed_id, 0.0)
        if cc < min_cc:
            Logger.info(
                f"Correlation of {cc} less than {min_cc} for {seed_id}, "
                "skipping.")
            continue
        snr_x = snrs_1[seed_id]
        snr_y = snrs_2[seed_id]
        if not correct_mag_bias:
            cc = snr_x = snr_y = 1.0
        # Correct for CC and SNR-bias and add to relative_magnitudes
        # This is equation 10 from Schaff & Richards 2014:
        rel_mag = math.log10(amplitude_ratio) + math.log10(
            math.sqrt((1 + 1 / snr_y**2) / (1 + 1 / snr_x**2)) * cc)
        Logger.info(f"Channel: {seed_id} Magnitude change {rel_mag:.2f}")
        relative_magnitudes.update({seed_id: rel_mag})
    if return_correlations:
        return relative_magnitudes, correlations
    return relative_magnitudes


def amp_pick_event(event, st, inventory, chans=('Z',), var_wintype=True,
                   winlen=0.9, pre_pick=0.2, pre_filt=True, lowcut=1.0,
                   highcut=20.0, corners=4, min_snr=1.0, plot=False,
                   remove_old=False, ps_multiplier=0.34, velocity=False,
                   water_level=0, iaspei_standard=False):
    """
    Pick amplitudes for local magnitude for a single event.

    Looks for maximum peak-to-trough amplitude for a channel in a stream, and
    picks this amplitude and period.  There are a few things it does
    internally to stabilise the result:

        1. Applies a given filter to the data using obspy's bandpass filter.
        The filter applied is a time-domain digital SOS filter.
        This is often necessary for small magnitude earthquakes.  To correct
        for this filter later the gain of the filter at the period of the
        maximum amplitude is retrieved using scipy's sosfreqz, and used to
        divide the resulting picked amplitude.

        2. Picks the peak-to-trough amplitude, but records half of this to
        cope with possible DC offsets.

        3. The maximum amplitude within the given window is picked. Care must
        be taken to avoid including surface waves in the window;

        4. A variable window-length is used by default that takes into account
        P-S times if available, this is in an effort to include only the
        body waves.  When P-S times are not available the ps_multiplier
        variable is used, which defaults to 0.34 x hypocentral distance.

    :type event: obspy.core.event.event.Event
    :param event: Event to pick
    :type st: obspy.core.stream.Stream
    :param st: Stream associated with event
    :type inventory: obspy.core.inventory.Inventory
    :param inventory:
        Inventory containing response information for the stations in st.
    :type chans: tuple
    :param chans:
        Tuple of the components to pick on, e.g. (Z, 1, 2, N, E)
    :type var_wintype: bool
    :param var_wintype:
        If True, the winlen will be multiplied by the P-S time if both P and
        S picks are available, otherwise it will be multiplied by the
        hypocentral distance*ps_multiplier, defaults to True
    :type winlen: float
    :param winlen:
        Length of window, see above parameter, if var_wintype is False then
        this will be in seconds, otherwise it is the multiplier to the
        p-s time, defaults to 0.9.
    :type pre_pick: float
    :param pre_pick:
        Time before the s-pick to start the cut window, defaults to 0.2.
    :type pre_filt: bool
    :param pre_filt: To apply a pre-filter or not, defaults to True
    :type lowcut: float
    :param lowcut: Lowcut in Hz for the pre-filter, defaults to 1.0
    :type highcut: float
    :param highcut: Highcut in Hz for the pre-filter, defaults to 20.0
    :type corners: int
    :param corners: Number of corners to use in the pre-filter
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio to allow a pick - see note below on
        signal-to-noise ratio calculation.
    :type plot: bool
    :param plot: Turn plotting on or off.
    :type remove_old: bool
    :param remove_old:
        If True, will remove old amplitudes and associated picks from event
        and overwrite with new picks. Defaults to False.
    :type ps_multiplier: float
    :param ps_multiplier:
        A p-s time multiplier of hypocentral distance - defaults to 0.34,
        based on p-s ratio of 1.68 and an S-velocity 0f 1.5km/s, deliberately
        chosen to be quite slow.
    :type velocity: bool
    :param velocity:
        Whether to make the pick in velocity space or not. Original definition
        of local magnitude used displacement of Wood-Anderson, MLv in seiscomp
        and Antelope uses a velocity measurement. *velocity and iaspei_standard
        are mutually exclusive*.
    :type water_level: float
    :param water_level:
        Water-level for seismometer simulation, see
        https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html
    :type iaspei_standard: bool
    :param iaspei_standard:
        Whether to output amplitude in IASPEI standard IAML (wood-anderson
        static amplification of 1), or AML with wood-anderson static
        amplification of 2080. Note: Units are SI (and specified in the
        amplitude)

    :returns: Picked event
    :rtype: :class:`obspy.core.event.Event`

    .. Note::
        Signal-to-noise ratio is calculated using the filtered data by
        dividing the maximum amplitude in the signal window (pick window)
        by the normalized noise amplitude (taken from the whole window
        supplied).

    .. Note::
        With `iaspei_standard=False`, picks will be returned in SI units
        (m or m/s), with the standard Wood-Anderson sensitivity of 2080 applied
        such that the measurements reflect the amplitude measured on a Wood
        Anderson instrument, as per the original local magnitude definitions
        of Richter and others.
    """
    if iaspei_standard and velocity:
        raise NotImplementedError("Velocity is not IASPEI standard for IAML.")
    try:
        event_origin = event.preferred_origin() or event.origins[0]
    except IndexError:
        event_origin = Origin()
    depth = event_origin.depth
    if depth is None:
        Logger.warning("No depth for the event, setting to 0 km")
        depth = 0

    # Remove amplitudes and picks for those amplitudes - this is not always
    # safe: picks may not be exclusively linked to amplitudes - hence the
    # default is *not* to do this.
    if remove_old and event.amplitudes:
        removal_ids = {amp.pick_id for amp in event.amplitudes}
        event.picks = [
            p for p in event.picks if p.resource_id not in removal_ids]
        event.amplitudes = []

    # We just want to look at P and S picks.
    picks = [p for p in event.picks
             if p.phase_hint and p.phase_hint[0].upper() in ("P", "S")]
    if len(picks) == 0:
        Logger.warning('No P or S picks found')
        return event

    st = st.copy().merge()  # merge the data, just in case! Work on a copy.
    # For each station cut the window
    for sta in {p.waveform_id.station_code for p in picks}:
        for chan in chans:
            Logger.info(f'Working on {sta} {chan}')
            tr = st.select(station=sta, component=chan)
            if not tr:
                Logger.warning(f'{sta} {chan} not found in the stream.')
                continue
            tr = tr.merge()[0]
            # Apply the pre-filter
            if pre_filt:
                tr = tr.split().detrend('simple').merge(fill_value=0)[0]
                tr.filter('bandpass', freqmin=lowcut, freqmax=highcut,
                          corners=corners)
            tr = _sim_WA(tr, inventory, water_level=water_level,
                         velocity=velocity)
            if tr is None:  # None returned when no matching response is found
                continue

            # Get the distance from an appropriate arrival
            sta_picks = [p for p in picks if p.waveform_id.station_code == sta]
            distances = []
            for pick in sta_picks:
                distances += [
                    a.distance for a in event_origin.arrivals
                    if a.pick_id == pick.resource_id and
                    a.distance is not None]
            if len(distances) == 0:
                Logger.error(f"Arrivals for station: {sta} do not contain "
                             "distances. Have you located this event?")
                hypo_dist = None
            else:
                # They should all be the same, but take the mean to be sure...
                distance = sum(distances) / len(distances)
                hypo_dist = np.sqrt(
                    np.square(degrees2kilometers(distance)) +
                    np.square(depth / 1000))

            # Get the earliest P and S picks on this station
            phase_picks = {"P": None, "S": None}
            for _hint in phase_picks.keys():
                _picks = sorted(
                    [p for p in sta_picks if p.phase_hint[0].upper() == _hint],
                    key=lambda p: p.time)
                if len(_picks) > 0:
                    phase_picks[_hint] = _picks[0]
            p_pick = phase_picks["P"]
            s_pick = phase_picks["S"]
            # Get the window size.
            if var_wintype:
                if p_pick and s_pick:
                    p_time, s_time = p_pick.time, s_pick.time
                elif s_pick and hypo_dist:
                    s_time = s_pick.time
                    p_time = s_time - (hypo_dist * ps_multiplier)
                elif p_pick and hypo_dist:
                    p_time = p_pick.time
                    s_time = p_time + (hypo_dist * ps_multiplier)
                elif (s_pick or p_pick) and hypo_dist is None:
                    Logger.error(
                        "No hypocentral distance and no matching P and S "
                        f"picks for {sta}, skipping.")
                    continue
                else:
                    raise NotImplementedError(
                        "No p or s picks - you should not have been able to "
                        "get here")
                trim_start = s_time - pre_pick
                trim_end = s_time + (s_time - p_time) * winlen
                # Work out the window length based on p-s time or distance
            else:  # Fixed window-length
                if s_pick:
                    s_time = s_pick.time
                elif p_pick and hypo_dist:
                    # In this case, there is no S-pick and the window length is
                    # fixed we need to calculate an expected S_pick based on
                    # the hypocentral distance, this will be quite hand-wavey
                    # as we are not using any kind of velocity model.
                    s_time = p_pick.time + hypo_dist * ps_multiplier
                else:
                    Logger.warning(
                        "No s-pick or hypocentral distance to predict "
                        f"s-arrival for station {sta}, skipping")
                    continue
                trim_start = s_time - pre_pick
                trim_end = s_time + winlen
            tr = tr.trim(trim_start, trim_end)
            if len(tr.data) <= 10:
                Logger.warning(f'Insufficient data for {sta}')
                continue
            # Get the amplitude
            try:
                amplitude, period, delay, peak, trough = _max_p2t(
                    tr.data, tr.stats.delta, return_peak_trough=True)
            except ValueError as e:
                Logger.error(e)
                Logger.error(f'No amplitude picked for tr {tr.id}')
                continue
            # Calculate the normalized noise amplitude
            snr = amplitude / np.sqrt(np.mean(np.square(tr.data)))
            if amplitude == 0.0:
                continue
            if snr < min_snr:
                Logger.info(
                    f'Signal to noise ratio of {snr} is below threshold.')
                continue
            if plot:
                plt.plot(np.arange(len(tr.data)), tr.data, 'k')
                plt.scatter(tr.stats.sampling_rate * delay, peak)
                plt.scatter(tr.stats.sampling_rate * (delay + period / 2),
                            trough)
                plt.show()
            Logger.info(f'Amplitude picked: {amplitude}')
            Logger.info(f'Signal-to-noise ratio is: {snr}')
            # Note, amplitude should be in meters at the moment!
            # Remove the pre-filter response
            if pre_filt:
                # Generate poles and zeros for the filter we used earlier.
                # We need to get the gain for the digital SOS filter used by
                # obspy.
                sos = iirfilter(
                    corners, [lowcut / (0.5 * tr.stats.sampling_rate),
                              highcut / (0.5 * tr.stats.sampling_rate)],
                    btype='band', ftype='butter', output='sos')
                _, gain = sosfreqz(sos, worN=[1 / period],
                                   fs=tr.stats.sampling_rate)
                gain = np.abs(gain[0])  # Convert from complex to real.
                if gain < 1e-2:
                    Logger.warning(
                        f"Pick made outside stable pass-band of filter "
                        f"on {tr.id}, rejecting")
                    continue
                amplitude /= gain
                Logger.debug(f"Removed filter gain: {gain}")
            # Write out the half amplitude, approximately the peak amplitude as
            # used directly in magnitude calculations
            amplitude *= 0.5
            # Documentation standards
            module = _sim_WA.__module__
            fname = currentframe().f_code.co_name
            # This is here to ensure that if the function name changes this
            # is still correct
            method_id = ResourceIdentifier(
                id=f"{module}.{fname}",
                prefix=f"smi:eqcorrscan{eqcorrscan.__version__}")
            filter_id = ResourceIdentifier(
                id=f"{module}._sim_WA",
                prefix=f"smi:eqcorrscan{eqcorrscan.__version__}")
            if iaspei_standard:
                # Remove wood-anderson amplification
                units, phase_hint, amplitude_type = (
                    "m", "IAML", "IAML")
                # amplitude *= 10 ** 9  # *THIS IS NOT SUPPORTED BY QML*
                amplitude /= PAZ_WA["sensitivity"]  # Remove WA sensitivity
                # Set the filter ID to state that sensitivity was removed
                filter_id = ResourceIdentifier(
                    id=f"{module}._sim_WA.WA_sensitivity_removed",
                    prefix=f"smi:eqcorrscan{eqcorrscan.__version__}")
            else:  # Not IAML, use SI units.
                if velocity:
                    units, phase_hint, amplitude_type = (
                        "m/s", "AML", "AML")
                else:
                    units, phase_hint, amplitude_type = (
                        "m", "AML", "AML")
            if tr.stats.channel.endswith("Z"):
                magnitude_hint = "MLv"
                # MLv is ML picked on the vertical channel
            else:
                magnitude_hint = "ML"
            # Append an amplitude reading to the event
            _waveform_id = WaveformStreamID(
                station_code=tr.stats.station, channel_code=tr.stats.channel,
                network_code=tr.stats.network)
            pick = Pick(
                waveform_id=_waveform_id, phase_hint=phase_hint,
                polarity='undecidable', time=tr.stats.starttime + delay,
                evaluation_mode='automatic',
                method_id=method_id, filter_id=filter_id)
            event.picks.append(pick)
            event.amplitudes.append(Amplitude(
                generic_amplitude=amplitude, period=period,
                pick_id=pick.resource_id, waveform_id=pick.waveform_id,
                unit=units, magnitude_hint=magnitude_hint,
                type=amplitude_type, category='point', method_id=method_id,
                filter_id=filter_id))
    return event


def svd_moments(u, s, v, stachans, event_list, n_svs=2):
    """
    Calculate relative moments/amplitudes using singular-value decomposition.

    Convert basis vectors calculated by singular value decomposition (see the
    SVD functions in clustering) into relative moments.

    For more information see the paper by
    `Rubinstein & Ellsworth (2010).
    <http://www.bssaonline.org/content/100/5A/1952.short>`_

    :type u: list
    :param u:
        List of the :class:`numpy.ndarray` input basis vectors from the SVD,
        one array for each channel used.
    :type s: list
    :param s:
        List of the :class:`numpy.ndarray` of singular values, one array for
        each channel.
    :type v: list
    :param v:
        List of :class:`numpy.ndarray` of output basis vectors from SVD, one
        array per channel.
    :type stachans: list
    :param stachans: List of station.channel input
    :type event_list: list
    :param event_list:
        List of events for which you have data, such that event_list[i]
        corresponds to stachans[i], U[i] etc. and event_list[i][j] corresponds
        to event j in U[i].  These are a series of indexes that map the basis
        vectors to their relative events and channels - if you have every
        channel for every event generating these is trivial (see example).
    :type n_svs: int
    :param n_svs: Number of singular values to use, defaults to 4.

    :returns: M, array of relative moments
    :rtype: :class:`numpy.ndarray`
    :returns: events_out, list of events that relate to M (in order), \
        does not include the magnitude information in the events, see note.
    :rtype: :class:`obspy.core.event.event.Event`

    .. note:: M is an array of relative moments (or amplitudes), these cannot
        be directly compared to true moments without calibration.

    .. note:: When comparing this method with the method used for creation
        of subspace detectors (Harris 2006) it is important to note that the
        input `design set` matrix in Harris contains waveforms as columns,
        whereas in Rubinstein & Ellsworth it contains waveforms as rows
        (i.e. the transpose of the Harris data matrix). The U and V matrices
        are therefore swapped between the two approaches. This is accounted
        for in EQcorrscan but may lead to confusion when reviewing the code.
        Here we use the Harris approach.

    .. rubric:: Example

    >>> from eqcorrscan.utils.mag_calc import svd_moments
    >>> from obspy import read
    >>> import glob
    >>> import os
    >>> from eqcorrscan.utils.clustering import svd
    >>> import numpy as np
    >>> # Do the set-up
    >>> testing_path = 'eqcorrscan/tests/test_data/similar_events_processed'
    >>> stream_files = glob.glob(os.path.join(testing_path, '*'))
    >>> stream_list = [read(stream_file) for stream_file in stream_files]
    >>> event_list = []
    >>> remove_list = [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]
    >>> for i, stream in enumerate(stream_list):
    ...     st_list = []
    ...     for tr in stream:
    ...         if (tr.stats.station, tr.stats.channel) not in remove_list:
    ...             stream.remove(tr)
    ...             continue
    ...         st_list.append(i)
    ...     event_list.append(st_list) # doctest: +SKIP
    >>> event_list = np.asarray(event_list).T.tolist()
    >>> SVec, SVal, U, stachans = svd(stream_list=stream_list) # doctest: +SKIP
    ['GCSZ.EHZ', 'WV04.SHZ', 'WHAT2.SH1']
    >>> M, events_out = svd_moments(u=U, s=SVal, v=SVec, stachans=stachans,
    ...                             event_list=event_list) # doctest: +SKIP

    """
    Logger.critical(
        "Proceed with caution: this function is experimental and somewhat"
        " stochastic - you should run this multiple times to ensure you get"
        " a stable result.")
    # Define maximum number of events, will be the width of K
    K_width = max([max(ev_list) for ev_list in event_list]) + 1
    # Sometimes the randomisation generates a singular matrix - rather than
    # attempting to regulerize this matrix I propose undertaking the
    # randomisation step a further time
    if len(stachans) == 1:
        Logger.critical('Only provided data from one station-channel - '
                        'will not try to invert')
        return u[0][:, 0], event_list[0]
    for i, stachan in enumerate(stachans):
        k = []  # Small kernel matrix for one station - channel
        # Copy the relevant vectors so as not to destroy them
        # Here we'll swap into the Rubinstein U and V matrices
        U_working = copy.deepcopy(v[i].T)
        V_working = copy.deepcopy(u[i])
        s_working = copy.deepcopy(s[i].T)
        ev_list = event_list[i]
        if len(ev_list) > len(U_working):
            Logger.error('U is : ' + str(U_working.shape))
            Logger.error('ev_list is len %s' % str(len(ev_list)))
            f_dump = open('mag_calc_U_working.pkl', 'wb')
            pickle.dump(U_working, f_dump)
            f_dump.close()
            raise IOError('More events than represented in U')
        # Set all non-important singular values to zero
        s_working[n_svs:len(s_working)] = 0
        s_working = np.diag(s_working)
        # Convert to numpy matrices
        U_working = np.matrix(U_working)
        V_working = np.matrix(V_working)
        s_working = np.matrix(s_working)
        SVD_weights = U_working[:, 0]
        # If all the weights are negative take the abs
        if np.all(SVD_weights < 0):
            Logger.warning('All weights are negative - flipping them')
            SVD_weights = np.abs(SVD_weights)
        SVD_weights = np.array(SVD_weights).reshape(-1).tolist()
        # Shuffle the SVD_weights prior to pairing - will give one of multiple
        # pairwise options - see p1956 of Rubinstein & Ellsworth 2010
        # We need to keep the real indexes though, otherwise, if there are
        # multiple events with the same weight we will end up with multiple
        # -1 values
        random_SVD_weights = np.copy(SVD_weights)
        # Tack on the indexes
        random_SVD_weights = random_SVD_weights.tolist()
        random_SVD_weights = [(random_SVD_weights[_i], _i)
                              for _i in range(len(random_SVD_weights))]
        random.shuffle(random_SVD_weights)
        # Add the first element to the end so all elements will be paired twice
        random_SVD_weights.append(random_SVD_weights[0])
        # Take pairs of all the SVD_weights (each weight appears in 2 pairs)
        pairs = []
        for pair in _pairwise(random_SVD_weights):
            pairs.append(pair)
        # Deciding values for each place in kernel matrix using the pairs
        for pairsIndex in range(len(pairs)):
            # We will normalize by the minimum weight
            _weights = list(zip(*list(pairs[pairsIndex])))[0]
            _indeces = list(zip(*list(pairs[pairsIndex])))[1]
            min_weight = min(np.abs(_weights))
            max_weight = max(np.abs(_weights))
            min_index = _indeces[np.argmin(np.abs(_weights))]
            max_index = _indeces[np.argmax(np.abs(_weights))]
            row = []
            # Working out values for each row of kernel matrix
            for j in range(len(SVD_weights)):
                if j == max_index:
                    result = -1
                elif j == min_index:
                    normalised = max_weight / min_weight
                    result = float(normalised)
                else:
                    result = 0
                row.append(result)
            # Add each row to the K matrix
            k.append(row)
        # k is now a square matrix, we need to flesh it out to be K_width
        k_filled = np.zeros([len(k), K_width])
        for j in range(len(k)):
            for l, ev in enumerate(ev_list):
                k_filled[j, ev] = k[j][l]
        if 'K' not in locals():
            K = k_filled
        else:
            K = np.concatenate([K, k_filled])
    # Remove any empty rows
    K_nonempty = []
    events_out = []
    for i in range(0, K_width):
        if not np.all(K[:, i] == 0):
            K_nonempty.append(K[:, i])
            events_out.append(i)
    K = np.array(K_nonempty).T
    K = K.tolist()
    K_width = len(K[0])
    # Add an extra row to K, so average moment = 1
    K.append(np.ones(K_width) * (1. / K_width))
    Logger.debug("Created Kernel matrix: ")
    del row
    Logger.debug('\n'.join([''.join([str(round(float(item), 3)).ljust(6)
                 for item in row]) for row in K]))
    Krounded = np.around(K, decimals=4)
    # Create a weighting matrix to put emphasis on the final row.
    W = np.matrix(np.identity(len(K)))
    # the final element of W = the number of stations*number of events
    W[-1, -1] = len(K) - 1
    # Make K into a matrix
    K = np.matrix(K)
    ############

    # Solve using the weighted least squares equation, K.T is K transpose
    Kinv = np.array(np.linalg.inv(K.T * W * K) * K.T * W)

    # M are the relative moments of the events
    M = Kinv[:, -1]
    # XXX TODO This still needs an outlier removal step
    return M, events_out


if __name__ == "__main__":
    import doctest
    doctest.testmod()
