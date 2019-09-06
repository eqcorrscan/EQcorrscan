"""Functions to calculate local magnitudes automatically, and to calcualte \
relative moments for near-repeating earthquakes using singular-value \
decomposition techniques.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import logging
import os
import glob
import matplotlib.pyplot as plt
import datetime as dt
import itertools
import sys
import copy
import random
import pickle
import math

from scipy.signal import iirfilter
from collections import Counter
from obspy.signal.invsim import simulate_seismometer as seis_sim
from obspy.signal.invsim import evalresp, paz_2_amplitude_value_of_freq_resp
from obspy import UTCDateTime
from obspy.core.event import Amplitude, Pick, WaveformStreamID
from obspy.geodetics import degrees2kilometers

from eqcorrscan.core.match_filter.matched_filter import MatchFilterError
from eqcorrscan.utils.catalog_utils import _get_origin


Logger = logging.getLogger(__name__)


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


def calc_max_curv(magnitudes, plotvar=False):
    """
    Calculate the magnitude of completeness using the maximum curvature method.

    :type magnitudes: list
    :param magnitudes:
        List of magnitudes from which to compute the maximum curvature which
        will give an estimate of the magnitude of completeness given the
        assumption of a power-law scaling.
    :type plotvar: bool
    :param plotvar: Turn plotting on and off

    :rtype: float
    :return: Magnitude at maximum curvature

    .. Note:: Should be used as a guide, often under-estimates Mc.

    .. rubric:: Example

    >>> import numpy as np
    >>> mags = []
    >>> for mag in np.arange(2.5,3, 0.1):
    ...     mags.extend([mag] * int(20000 - 10 * mag))
    >>> for mag in np.arange(3,7, 0.1):
    ...     mags.extend([mag] * int(10 ** (7 - 1 * mag)))
    >>> calc_max_curv(mags, plotvar=False)
    3.0
    """
    counts = Counter(magnitudes)
    df = np.zeros(len(counts))
    mag_steps = np.zeros(len(counts))
    grad = np.zeros(len(counts) - 1)
    grad_points = grad.copy()
    for i, magnitude in enumerate(sorted(counts.keys(), reverse=True)):
        mag_steps[i] = magnitude
        if i > 0:
            df[i] = counts[magnitude] + df[i - 1]
        else:
            df[i] = counts[magnitude]
    for i, val in enumerate(df):
        if i > 0:
            grad[i - 1] = (val - df[i - 1]) / (mag_steps[i] - mag_steps[i - 1])
            grad_points[i - 1] = mag_steps[i] - ((mag_steps[i] -
                                                  mag_steps[i - 1]) / 2.0)
    # Need to find the second order derivative
    curvature = np.zeros(len(grad) - 1)
    curvature_points = curvature.copy()
    for i, _grad in enumerate(grad):
        if i > 0:
            curvature[i - 1] = (_grad - grad[i - 1]) / (grad_points[i] -
                                                        grad_points[i - 1])
            curvature_points[i - 1] = grad_points[i] - ((grad_points[i] -
                                                         grad_points[i - 1]) /
                                                        2.0)
    if plotvar:
        plt.scatter(mag_steps, df, c='k', label='Magnitude function')
        plt.plot(mag_steps, df, c='k')
        plt.scatter(grad_points, grad, c='r', label='Gradient')
        plt.plot(grad_points, grad, c='r')
        plt.scatter(curvature_points, curvature, c='g', label='Curvature')
        plt.plot(curvature_points, curvature, c='g')
        plt.legend()
        plt.show()
    return curvature_points[np.argmax(abs(curvature))]


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
    :return: List of tuples of (completeness, b-value, residual,\
        number of magnitudes used)

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
    >>> round(b_values[4][1])
    1.0
    >>> # We can set a maximum magnitude:
    >>> b_values = calc_b_value(magnitudes, completeness=np.arange(3, 7, 0.2),
    ...                         plotvar=False, max_mag=5)
    >>> round(b_values[4][1])
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
        b_values.append((m_c, abs(fit[0][0]), r, str(len(complete_mags))))
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


def _sim_WA(trace, PAZ, seedresp, water_level, velocity=False):
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
    :type PAZ: dict
    :param PAZ:
        Dictionary containing lists of poles and zeros, the gain and
        the sensitivity. If unset will expect seedresp.
    :type seedresp: dict
    :param seedresp: Seed response information - if unset will expect PAZ.
    :type water_level: int
    :param water_level: Water level for the simulation.
    :type velocity: bool
    :param velocity:
        Whether to return a velocity trace or not - velocity is non-standard
        for Wood-Anderson instruments, but institutes that use seiscomp3 or
        Antelope require picks in velocity.

    :returns: Trace of Wood-Anderson simulated data
    :rtype: :class:`obspy.core.trace.Trace`
    """
    # Note Wood anderson sensitivity is 2080 as per Uhrhammer & Collins 1990
    PAZ_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
              'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    if velocity:
        PAZ_WA['zeros'] = [0 + 0j, 0 + 0j]
    # De-trend data
    trace.detrend('simple')
    # Simulate Wood Anderson
    if PAZ:
        trace.data = seis_sim(trace.data, trace.stats.sampling_rate,
                              paz_remove=PAZ, paz_simulate=PAZ_WA,
                              water_level=water_level,
                              remove_sensitivity=True)
    elif seedresp:
        trace.data = seis_sim(trace.data, trace.stats.sampling_rate,
                              paz_remove=None, paz_simulate=PAZ_WA,
                              water_level=water_level, seedresp=seedresp)
    else:
        Logger.warning('No response given to remove, will just simulate WA')
        trace.data = seis_sim(trace.data, trace.stats.sampling_rate,
                              paz_remove=None, paz_simulate=PAZ_WA,
                              water_level=water_level)
    return trace


def _max_p2t(data, delta):
    """
    Finds the maximum peak-to-trough amplitude and period.
    Originally designed to be used to calculate magnitudes (by \
    taking half of the peak-to-trough amplitude as the peak amplitude).

    :type data: numpy.ndarray
    :param data: waveform trace to find the peak-to-trough in.
    :type delta: float
    :param delta: Sampling interval in seconds

    :returns: tuple of (amplitude, period, time) with amplitude in the same \
        scale as given in the input data, and period in seconds, and time in \
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
    return amplitude, period, delta * turning_points[np.argmax(amplitudes)][1]


def _GSE2_PAZ_read(gsefile):
    """
    Read the instrument response information from a GSE Poles and Zeros file.

    Formatted for files generated by the SEISAN program RESP.

    Format must be CAL2, not coded for any other format at the moment,
    contact the authors to add others in.

    :type gsefile: string
    :param gsefile: Path to GSE file

    :returns: Dictionary of poles, zeros, gain and sensitivity
    :rtype: dict
    """
    with open(gsefile, 'r') as f:
        # First line should start with CAL2
        header = f.readline()
        if not header[0:4] == 'CAL2':
            raise IOError('Unknown format for GSE file, only coded for CAL2')
        station = header.split()[1]
        channel = header.split()[2]
        sensor = header.split()[3]
        date = dt.datetime.strptime(header.split()[7], '%Y/%m/%d')
        header = f.readline()
        if not header[0:4] == 'PAZ2':
            raise IOError('Unknown format for GSE file, only coded for PAZ2')
        gain = float(header.split()[3])  # Measured in nm/counts
        kpoles = int(header.split()[4])
        kzeros = int(header.split()[5])
        poles = []
        for i in range(kpoles):
            pole = f.readline()
            poles.append(complex(float(pole.split()[0]),
                                 float(pole.split()[1])))
        zeros = []
        for i in range(kzeros):
            zero = f.readline()
            zeros.append(complex(float(zero.split()[0]),
                                 float(zero.split()[1])))
        # Have Poles and Zeros, but need Gain and Sensitivity
        # Gain should be in the DIG2 line:
        for line in f:
            if line[0:4] == 'DIG2':
                sensitivity = float(line.split()[2])
                # measured in counts/muVolt
    PAZ = {'poles': poles, 'zeros': zeros, 'gain': gain,
           'sensitivity': sensitivity}
    return PAZ, date, station, channel, sensor


def _find_resp(station, channel, network, time, delta, directory):
    """
    Helper function to find the response information.

    Works for a given station and channel at a given time and return a
    dictionary of poles and zeros, gain and sensitivity.

    :type station: str
    :param station: Station name (as in the response files)
    :type channel: str
    :param channel: Channel name (as in the response files)
    :type network: str
    :param network: Network to scan for, can be a wildcard
    :type time: datetime.datetime
    :param time: Date-time to look for repsonse information
    :type delta: float
    :param delta: Sample interval in seconds
    :type directory: str
    :param directory: Directory to scan for response information

    :returns: dictionary of response information
    :rtype: dict
    """
    possible_respfiles = glob.glob(directory + os.path.sep + 'RESP.' +
                                   network + '.' + station +
                                   '.*.' + channel)  # GeoNet RESP naming
    possible_respfiles += glob.glob(directory + os.path.sep + 'RESP.' +
                                    network + '.' + channel +
                                    '.' + station)  # RDseed RESP naming
    possible_respfiles += glob.glob(directory + os.path.sep + 'RESP.' +
                                    station + '.' + network)
    # WIZARD resp naming
    # GSE format, station needs to be 5 characters padded with _, channel is 4
    # characters padded with _
    station = str(station)
    channel = str(channel)
    possible_respfiles += glob.glob(directory + os.path.sep +
                                    station.ljust(5, str('_')) +
                                    channel[0:len(channel) - 1].
                                    ljust(3, str('_')) +
                                    channel[-1] + '.*_GSE')
    PAZ = []
    seedresp = []
    for respfile in possible_respfiles:
        Logger.info('Reading response from: ' + respfile)
        if respfile.split(os.path.sep)[-1][0:4] == 'RESP':
            # Read from a resp file
            seedresp = {'filename': respfile, 'date': UTCDateTime(time),
                        'units': 'DIS', 'network': network, 'station': station,
                        'channel': channel, 'location': '*'}
            try:
                # Attempt to evaluate the response for this information, if not
                # then this is not the correct response info!
                freq_resp, freqs = evalresp(
                    delta, 100, seedresp['filename'], seedresp['date'],
                    units=seedresp['units'], freq=True,
                    network=seedresp['network'], station=seedresp['station'],
                    channel=seedresp['channel'])
            except Exception as e:
                Logger.warning('Issues with RESP file: {0}'.format(e))
                seedresp = []
                continue
        elif respfile[-3:] == 'GSE':
            PAZ, pazdate, pazstation, pazchannel, pazsensor =\
                _GSE2_PAZ_read(respfile)
            # check that the date is good!
            if pazdate >= time and pazchannel != channel and\
               pazstation != station:
                Logger.warning('Issue with GSE file: date: ' +
                               str(pazdate) + ' channel: ' + pazchannel +
                               ' station: ' + pazstation)
                PAZ = []
        else:
            continue
        # Check that PAZ are for the correct station, channel and date
        if PAZ or seedresp:
            break
    if PAZ:
        return PAZ
    elif seedresp:
        return seedresp


def _pairwise(iterable):
    """
    Wrapper on itertools for SVD_magnitude.
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    if sys.version_info.major == 2:
        return itertools.izip(a, b)
    else:
        return zip(a, b)


def _snr(tr, noise_window, signal_window):
    """
    Compute ratio of maximum signal amplitude to rms noise amplitude.

    :type tr: `obspy.core.Trace`
    :param tr: Trace to compute signal-to-noise ratio for
    :type noise_window: tuple of UTCDateTime
    :param noise_window: (start, end) of window to use for noise
    :type signal_window: tuple of UTCDateTime
    :param signal_window: (start, end) of window to use for signal

    :rtype: float
    :return: Signal-to-noise ratio.
    """
    from eqcorrscan.core.template_gen import _rms

    noise_amp = _rms(
        tr.slice(starttime=noise_window[0], endtime=noise_window[1]).data)
    if np.isnan(noise_amp):
        Logger.info("Noise amplitude is nan, setting to 1")
        noise_amp = 1.
    try:
        signal_amp = tr.slice(
            starttime=signal_window[0], endtime=signal_window[1]).data.max()
    except ValueError as e:
        Logger.error(e)
        return np.nan
    return signal_amp / noise_amp


def relative_amplitude(st1, st2, event1, event2, noise_window=(-20, -1),
                       signal_window=(-.5, 20), min_snr=5.0,
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
    trace. Windows are computed relative to the first pick for that station.

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

    :rtype: dict
    :return: Dictionary of relative amplitudes keyed by seed-id
    """
    from obspy import Stream

    amplitudes = {}
    for tr1 in st1:
        pick1 = _get_pick_for_station(
            event=event1, station=tr1.stats.station, use_s_picks=use_s_picks)
        if pick1 is None:
            continue
        snr1 = _snr(
            tr1, (pick1.time + noise_window[0], pick1.time + noise_window[1]),
            (pick1.time + signal_window[0], pick1.time + signal_window[1]))
        if np.isnan(snr1) or snr1 <= min_snr:
            Logger.info(
                "SNR of {0} is below min_snr ({1}) for {2} in st1".format(
                    snr1, min_snr, tr1.id))
            continue
        tr2 = [tr for tr in st2 if tr.id == tr1.id]
        if len(tr2) == 0:
            Logger.info("No matched traces for {0}".format(tr1.id))
            continue
        tr2 = Stream(tr2).merge()[0]
        pick2 = _get_pick_for_station(
            event=event2, station=tr2.stats.station, use_s_picks=use_s_picks)
        if pick2 is None:
            continue
        snr2 = _snr(
            tr2, (pick2.time + noise_window[0], pick2.time + noise_window[1]),
            (pick2.time + signal_window[0], pick2.time + signal_window[1]))
        if np.isnan(snr2) or snr2 <= min_snr:
            Logger.info(
                "SNR of {0} is below min_snr ({1}) for {2} in st2".format(
                    snr2, min_snr, tr2.id))
            continue
        # If we get here, actually compute the ratio in the signal windows
        amp1 = tr1.slice(starttime=pick1.time + signal_window[0],
                         endtime=pick1.time + signal_window[1]).data.std()
        amp2 = tr2.slice(starttime=pick2.time + signal_window[0],
                         endtime=pick2.time + signal_window[1]).data.std()
        Logger.debug("Channel: {0} Relative amplitude: {1:.2f}".format(
            tr1.id, amp2/amp1))
        amplitudes.update({tr1.id: amp2 / amp1})
    return amplitudes


def _get_pick_for_station(event, station, use_s_picks):
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
    picks = [p for p in event.picks if p.waveform_id.station_code == station]
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


def relative_magnitude(st1, st2, event1, event2, noise_window=(-20, -1),
                       signal_window=(-.5, 20), min_snr=5.0, min_cc=0.7,
                       use_s_picks=False, correlations=None, shift=.2,
                       return_correlations=False):
    """
    Compute the relative magnitudes between two events.

    See :func:`eqcorrscan.utils.mag_calc.relative_amplitude` for information
    on how relative amplitudes are calculated. To compute relative magnitudes
    from relative amplitudes this function weights the amplitude ratios by
    the cross-correlation of the two events. The relation used is similar to
    Schaff and Richards (2014) and is:

    .. math::

        \\Delta m = \\log{\\frac{std(tr2)}{std(tr1)}} \\times CC

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
    relative_amplitudes = relative_amplitude(
        st1=st1, st2=st2, event1=event1, event2=event2,
        noise_window=noise_window, signal_window=signal_window,
        min_snr=min_snr, use_s_picks=use_s_picks)
    for seed_id, amplitude_ratio in relative_amplitudes.items():
        tr1 = st1.select(id=seed_id)[0]
        tr2 = st2.select(id=seed_id)[0]
        pick1 = _get_pick_for_station(
            event=event1, station=tr1.stats.station, use_s_picks=use_s_picks)
        pick2 = _get_pick_for_station(
            event=event2, station=tr2.stats.station, use_s_picks=use_s_picks)
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
            continue
        # Weight and add to relative_magnitudes
        rel_mag = math.log10(amplitude_ratio) * cc
        Logger.debug("Channel: {0} Magnitude change {1:.2f}".format(
            tr1.id, rel_mag))
        relative_magnitudes.update({seed_id: rel_mag})
    if return_correlations:
        return relative_magnitudes, correlations
    return relative_magnitudes


def amp_pick_event(event, st, respdir, chans=['Z'], var_wintype=True,
                   winlen=0.9, pre_pick=0.2, pre_filt=True, lowcut=1.0,
                   highcut=20.0, corners=4, min_snr=1.0, plot=False,
                   remove_old=False, ps_multiplier=0.34, velocity=False):
    """
    Pick amplitudes for local magnitude for a single event.

    Looks for maximum peak-to-trough amplitude for a channel in a stream, and
    picks this amplitude and period.  There are a few things it does
    internally to stabilise the result:

        1. Applies a given filter to the data - very necessary for small
        magnitude earthquakes;

        2. Keeps track of the poles and zeros of this filter and removes them
        from the picked amplitude;

        3. Picks the peak-to-trough amplitude, but records half of this: the
        specification for the local magnitude is to use a peak amplitude on
        a horizontal, however, with modern digital seismometers, the peak
        amplitude often has an additional, DC-shift applied to it, to
        stabilise this, and to remove possible issues with de-meaning data
        recorded during the wave-train of an event (e.g. the mean may not be
        the same as it would be for longer durations), we use half the
        peak-to-trough amplitude;

        4. Despite the original definition of local magnitude requiring the
        use of a horizontal channel, more recent work has shown that the
        vertical channels give more consistent magnitude estimations between
        stations, due to a reduction in site-amplification effects, we
        therefore use the vertical channels by default, but allow the user
        to chose which channels they deem appropriate;

        5. We do not specify that the maximum amplitude should be the
        S-phase: The original definition holds that the maximum body-wave
        amplitude should be used - while this is often the S-phase, we do not
        discriminate against the P-phase.  We do note that, unless the user
        takes care when assigning winlen and filters, they may end up with
        amplitude picks for surface waves;

        6. We use a variable window-length by default that takes into account
        P-S times if available, this is in an effort to include only the
        body waves.  When P-S times are not available we us the ps_multiplier
        variable, which defaults to 0.34 x hypocentral distance.

    :type event: obspy.core.event.event.Event
    :param event: Event to pick
    :type st: obspy.core.stream.Stream
    :param st: Stream associated with event
    :type respdir: str
    :param respdir: Path to the response information directory
    :type chans: list
    :param chans:
        List of the channels to pick on, defaults to ['Z'] - should just be
        the orientations, e.g. Z, 1, 2, N, E
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
        If True, will remove old amplitude picks from event and overwrite
        with new picks. Defaults to False.
    :type ps_multiplier: float
    :param ps_multiplier:
        A p-s time multiplier of hypocentral distance - defaults to 0.34,
        based on p-s ratio of 1.68 and an S-velocity 0f 1.5km/s, deliberately
        chosen to be quite slow.
    :type velocity: bool
    :param velocity:
        Whether to make the pick in velocity space or not. Original definition
        of local magnitude used displacement of Wood-Anderson, MLv in seiscomp
        and Antelope uses a velocity measurement.

    :returns: Picked event
    :rtype: :class:`obspy.core.event.Event`

    .. Note::
        Signal-to-noise ratio is calculated using the filtered data by
        dividing the maximum amplitude in the signal window (pick window)
        by the normalized noise amplitude (taken from the whole window
        supplied).

    .. Warning::
        Works in place on data - will filter and remove response from data,
        you are recommended to give this function a copy of the data if you
        are using it in a loop.
    """
    # Convert these picks into a lists
    stations = []  # List of stations
    channels = []  # List of channels
    picktimes = []  # List of pick times
    picktypes = []  # List of pick types
    picks_out = []
    try:
        depth = _get_origin(event).depth
    except MatchFilterError:
        depth = 0
    if remove_old and event.amplitudes:
        for amp in event.amplitudes:
            # Find the pick and remove it too
            pick = [p for p in event.picks if p.resource_id == amp.pick_id][0]
            event.picks.remove(pick)
        event.amplitudes = []
    for pick in event.picks:
        if pick.phase_hint in ['P', 'S']:
            picks_out.append(pick)  # Need to be able to remove this if there
            # isn't data for a station!
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code)
            picktimes.append(pick.time)
            picktypes.append(pick.phase_hint)
    if len(picktypes) == 0:
        Logger.warning('No P or S picks found')
    st.merge()  # merge the data, just in case!
    # For each station cut the window
    uniq_stas = list(set(stations))
    for sta in uniq_stas:
        for chan in chans:
            Logger.info('Working on ' + sta + ' ' + chan)
            tr = st.select(station=sta, channel='*' + chan)
            if not tr:
                Logger.warning(
                    'There is no station and channel match in the wavefile!')
                continue
            else:
                tr = tr[0]
            # Apply the pre-filter
            if pre_filt:
                try:
                    tr.split().detrend('simple').merge(fill_value=0)
                except Exception as e:
                    Logger.warning(
                        'Some issue splitting this one: {0}'.format(e))
                    dummy = tr.split()
                    dummy.detrend('simple')
                    tr = dummy.merge(fill_value=0)
                try:
                    tr.filter('bandpass', freqmin=lowcut, freqmax=highcut,
                              corners=corners)
                except NotImplementedError:
                    Logger.error(
                        'For some reason trace is not continuous: {0}'.format(
                            tr))
                    continue
            # Find the response information
            resp_info = _find_resp(
                tr.stats.station, tr.stats.channel, tr.stats.network,
                tr.stats.starttime, tr.stats.delta, respdir)
            PAZ = []
            seedresp = []
            if resp_info and 'gain' in resp_info:
                PAZ = resp_info
            elif resp_info:
                seedresp = resp_info
            # Simulate a Wood Anderson Seismograph
            if PAZ and len(tr.data) > 10:
                # Set ten data points to be the minimum to pass
                tr = _sim_WA(tr, PAZ, None, 10, velocity=velocity)
            elif seedresp and len(tr.data) > 10:
                tr = _sim_WA(tr, None, seedresp, 10, velocity=velocity)
            elif len(tr.data) > 10:
                Logger.warning('No PAZ for ' + tr.stats.station + ' ' +
                               tr.stats.channel + ' at time: ' +
                               str(tr.stats.starttime))
                continue
            sta_picks = [i for i in range(len(stations))
                         if stations[i] == sta]
            pick_id = event.picks[sta_picks[0]].resource_id
            arrival = [arrival for arrival in event.origins[0].arrivals
                       if arrival.pick_id == pick_id][0]
            hypo_dist = np.sqrt(
                np.square(degrees2kilometers(arrival.distance)) +
                np.square(depth / 1000))
            if var_wintype and hypo_dist:
                if 'S' in [picktypes[i] for i in sta_picks] and\
                   'P' in [picktypes[i] for i in sta_picks]:
                    # If there is an S-pick we can use this :D
                    s_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'S']
                    s_pick = min(s_pick)
                    p_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'P']
                    p_pick = min(p_pick)
                    try:
                        tr.trim(starttime=s_pick - pre_pick,
                                endtime=s_pick + (s_pick - p_pick) * winlen)
                    except ValueError:
                        continue
                elif 'S' in [picktypes[i] for i in sta_picks]:
                    s_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'S']
                    s_pick = min(s_pick)
                    p_modelled = s_pick - (hypo_dist * ps_multiplier)
                    try:
                        tr.trim(starttime=s_pick - pre_pick,
                                endtime=s_pick + (s_pick - p_modelled) *
                                winlen)
                    except ValueError:
                        continue
                else:
                    # In this case we only have a P pick
                    p_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'P']
                    p_pick = min(p_pick)
                    s_modelled = p_pick + (hypo_dist * ps_multiplier)
                    Logger.info('P_pick=%s' % str(p_pick))
                    Logger.info('hypo_dist: %s' % str(hypo_dist))
                    Logger.info('S modelled=%s' % str(s_modelled))
                    try:
                        tr.trim(starttime=s_modelled - pre_pick,
                                endtime=s_modelled + (s_modelled - p_pick) *
                                winlen)
                        Logger.debug(tr)
                    except ValueError:
                        continue
                # Work out the window length based on p-s time or distance
            elif 'S' in [picktypes[i] for i in sta_picks]:
                # If the window is fixed we still need to find the start time,
                # which can be based either on the S-pick (this elif), or
                # on the hypocentral distance and the P-pick

                # Take the minimum S-pick time if more than one S-pick is
                # available
                s_pick = [picktimes[i] for i in sta_picks
                          if picktypes[i] == 'S']
                s_pick = min(s_pick)
                try:
                    tr.trim(starttime=s_pick - pre_pick,
                            endtime=s_pick + winlen)
                except ValueError:
                    continue
            else:
                # In this case, there is no S-pick and the window length is
                # fixed we need to calculate an expected S_pick based on the
                # hypocentral distance, this will be quite hand-wavey as we
                # are not using any kind of velocity model.
                p_pick = [picktimes[i] for i in sta_picks
                          if picktypes[i] == 'P']
                Logger.debug(picktimes)
                p_pick = min(p_pick)
                s_modelled = p_pick + hypo_dist * ps_multiplier
                try:
                    tr.trim(starttime=s_modelled - pre_pick,
                            endtime=s_modelled + winlen)
                except ValueError:
                    continue
            if len(tr.data) <= 10:
                Logger.warning('No data found for: ' + tr.stats.station)
                continue
            # Get the amplitude
            try:
                amplitude, period, delay = _max_p2t(tr.data, tr.stats.delta)
            except ValueError:
                Logger.error('No amplitude picked for tr %s' % str(tr))
                continue
            # Calculate the normalized noise amplitude
            noise_amplitude = np.sqrt(np.mean(np.square(tr.data)))
            if amplitude == 0.0:
                continue
            if amplitude / noise_amplitude < min_snr:
                Logger.info(
                    'Signal to noise ratio of %s is below threshold.' %
                    (amplitude / noise_amplitude))
                continue
            if plot:
                plt.plot(np.arange(len(tr.data)), tr.data, 'k')
                plt.scatter(tr.stats.sampling_rate * delay, amplitude / 2)
                plt.scatter(tr.stats.sampling_rate * (delay + period),
                            -amplitude / 2)
                plt.show()
            Logger.info('Amplitude picked: ' + str(amplitude))
            Logger.info('Signal-to-noise ratio is: %s' %
                        (amplitude / noise_amplitude))
            # Note, amplitude should be in meters at the moment!
            # Remove the pre-filter response
            if pre_filt:
                # Generate poles and zeros for the filter we used earlier: this
                # is how the filter is designed in the convenience methods of
                # filtering in obspy.
                z, p, k = iirfilter(
                    corners, [lowcut / (0.5 * tr.stats.sampling_rate),
                              highcut / (0.5 * tr.stats.sampling_rate)],
                    btype='band', ftype='butter', output='zpk')
                filt_paz = {'poles': list(p), 'zeros': list(z), 'gain': k,
                            'sensitivity': 1.0}
                amplitude /= (paz_2_amplitude_value_of_freq_resp(
                    filt_paz, 1 / period) * filt_paz['sensitivity'])
            if PAZ:
                amplitude /= 1000
            if seedresp:  # Seedresp method returns mm
                amplitude *= 1000000
            # Write out the half amplitude, approximately the peak amplitude as
            # used directly in magnitude calculations
            amplitude *= 0.5
            # Append an amplitude reading to the event
            _waveform_id = WaveformStreamID(
                station_code=tr.stats.station, channel_code=tr.stats.channel,
                network_code=tr.stats.network)
            pick_ind = len(event.picks)
            event.picks.append(Pick(
                waveform_id=_waveform_id, phase_hint='IAML',
                polarity='undecidable', time=tr.stats.starttime + delay,
                evaluation_mode='automatic'))
            if not velocity:
                event.amplitudes.append(Amplitude(
                    generic_amplitude=amplitude / 1e9, period=period,
                    pick_id=event.picks[pick_ind].resource_id,
                    waveform_id=event.picks[pick_ind].waveform_id, unit='m',
                    magnitude_hint='ML', type='AML', category='point'))
            else:
                event.amplitudes.append(Amplitude(
                    generic_amplitude=amplitude / 1e9, period=period,
                    pick_id=event.picks[pick_ind].resource_id,
                    waveform_id=event.picks[pick_ind].waveform_id, unit='m/s',
                    magnitude_hint='ML', type='AML', category='point'))
    return event


def amp_pick_sfile(*args, **kwargs):
    raise ImportError(
        "Sfile support is depreciated, read in using obspy.io.nordic")


def SVD_moments(U, s, V, stachans, event_list, n_SVs=4):
    """Depreciated."""
    Logger.warning('Depreciated, use svd_moments instead')
    return svd_moments(u=U, s=s, v=V, stachans=stachans,
                       event_list=event_list, n_svs=n_SVs)


def svd_moments(u, s, v, stachans, event_list, n_svs=2):
    """
    Calculate relative moments/amplitudes using singular-value decomposition.

    Convert basis vectors calculated by singular value \
    decomposition (see the SVD functions in clustering) into relative \
    moments.

    For more information see the paper by \
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
    :param event_list: List of events for which you have data, such that \
        event_list[i] corresponds to stachans[i], U[i] etc. and \
        event_list[i][j] corresponds to event j in U[i].  These are a series \
        of indexes that map the basis vectors to their relative events and \
        channels - if you have every channel for every event generating these \
        is trivial (see example).
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
