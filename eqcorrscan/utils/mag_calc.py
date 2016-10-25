"""Functions to calculate local magnitudes automatically, and to calcualte \
relative moments for near-repeating earthquakes using singular-value \
decomposition techniques.

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
import os
import glob
import matplotlib.pyplot as plt
import datetime as dt
import itertools
import sys
import shutil
import copy
import random
import pickle

from scipy.signal import iirfilter
from collections import Counter
from obspy.signal.invsim import simulate_seismometer as seis_sim
from obspy.signal.invsim import evalresp, paz_2_amplitude_value_of_freq_resp
from obspy import UTCDateTime, read
from obspy.core.event import Amplitude, Pick, WaveformStreamID

from eqcorrscan.utils import sfile_util


def dist_calc(loc1, loc2):
    """
    Function to calculate the distance in km between two points.

    Uses the flat Earth approximation. Better things are available for this,
    like `gdal <http://www.gdal.org/>`_.

    :type loc1: tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: tuple
    :param loc2: Tuple of lat, lon, depth (in decimal degrees and km)

    :returns: Distance between points in km.
    :rtype: float
    """
    R = 6371.009  # Radius of the Earth in km
    dlat = np.radians(abs(loc1[0] - loc2[0]))
    dlong = np.radians(abs(loc1[1] - loc2[1]))
    ddepth = abs(loc1[2] - loc2[2])
    mean_lat = np.radians((loc1[0] + loc2[0]) / 2)
    dist = R * np.sqrt(dlat ** 2 + (np.cos(mean_lat) * dlong) ** 2)
    dist = np.sqrt(dist ** 2 + ddepth ** 2)
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

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.mag_calc import calc_max_curv
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime('2012-03-26T00:00:00')
    >>> t2 = t1 + (3 * 86400)
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
    >>> magnitudes = [event.magnitudes[0].mag for event in catalog]
    >>> calc_max_curv(magnitudes, plotvar=False)
    3.1000000000000001
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
            warnings.warn('Not computing completeness at %s, above max_mag' %
                          str(m_c))
            break
        complete_mags = []
        complete_freq = []
        for i, mag in enumerate(mag_steps):
            if mag >= m_c <= max_mag:
                complete_mags.append(mag)
                complete_freq.append(np.log10(cdf[i]))
        if len(complete_mags) < 4:
            warnings.warn('Not computing completeness above ' + str(m_c) +
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
        b_vals = ax1.scatter(zip(*b_values)[0], zip(*b_values)[1], c='k')
        resid = ax1.scatter(zip(*b_values)[0],
                            [100 - b for b in zip(*b_values)[2]], c='r')
        ax1.set_ylabel('b-value and residual')
        plt.xlabel('Completeness magnitude')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Number of events used in fit')
        n_ev = ax2.scatter(zip(*b_values)[0], zip(*b_values)[3], c='g')
        fig.legend((b_vals, resid, n_ev),
                   ('b-values', 'residuals', 'number of events'),
                   'lower right')
        ax1.set_title('Possible completeness values')
        plt.show()
    return b_values


def _sim_WA(trace, PAZ, seedresp, water_level):
    """
    Remove the instrument response from a trace and simulate a Wood-Anderson.

    Returns a de-meaned, de-trended, Wood Anderson simulated trace in
    it's place.

    Works in-place on data and will destroy your original data, copy the \
    trace before giving it to this function!

    :type trace: obspy.core.trace.Trace
    :param trace: A standard obspy trace, generally should be given without
                    pre-filtering, if given with pre-filtering for use with
                    amplitude determiniation for magnitudes you will need to
                    worry about how you cope with the response of this filter
                    yourself.
    :type PAZ: dict
    :param PAZ: Dictionary containing lists of poles and zeros, the gain and
                the sensitivity. If unset will expect seedresp.
    :type seedresp: dict
    :param seedresp: Seed response information - if unset will expect PAZ.
    :type water_level: int
    :param water_level: Water level for the simulation.

    :returns: Trace of Wood-Anderson simulated data
    :rtype: :class:`obspy.core.trace.Trace`
    """
    # Note Wood anderson sensitivity is 2080 as per Uhrhammer & Collins 1990
    PAZ_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
              'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    # De-trend data
    trace.detrend('simple')
    # Simulate Wood Anderson
    if PAZ:
        trace.data = seis_sim(trace.data, trace.stats.sampling_rate,
                              paz_remove=PAZ, paz_simulate=PAZ_WA,
                              water_level=water_level, remove_sensitivity=True)
    elif seedresp:
        trace.data = seis_sim(trace.data, trace.stats.sampling_rate,
                              paz_remove=None, paz_simulate=PAZ_WA,
                              water_level=water_level, seedresp=seedresp)
    else:
        UserWarning('No response given to remove, will just simulate WA')
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
    debug_plot = False
    turning_points = []  # A list of tuples of (amplitude, sample)
    for i in range(1, len(data) - 1):
        if (data[i] < data[i - 1] and data[i] < data[i + 1]) or\
           (data[i] > data[i - 1] and data[i] > data[i + 1]):
            turning_points.append((data[i], i))
    if len(turning_points) >= 1:
        amplitudes = np.empty([len(turning_points) - 1],)
        half_periods = np.empty([len(turning_points) - 1],)
    else:
        plt.plot(data)
        plt.show()
        print('Turning points has length: ' + str(len(turning_points)) +
              ' data have length: ' + str(len(data)))
        return 0.0, 0.0, 0.0
    for i in range(1, len(turning_points)):
        half_periods[i - 1] = (delta * (turning_points[i][1] -
                                        turning_points[i - 1][1]))
        amplitudes[i - 1] = np.abs(turning_points[i][0] -
                                   turning_points[i - 1][0])
    amplitude = np.max(amplitudes)
    period = 2 * half_periods[np.argmax(amplitudes)]
    if debug_plot:
        plt.plot(data, 'k')
        plt.plot([turning_points[np.argmax(amplitudes)][1],
                  turning_points[np.argmax(amplitudes) - 1][1]],
                 [turning_points[np.argmax(amplitudes)][0],
                  turning_points[np.argmax(amplitudes) - 1][0]], 'r')
        plt.show()
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
        print('Reading response from: ' + respfile)
        if respfile.split(os.path.sep)[-1][0:4] == 'RESP':
            # Read from a resp file
            seedresp = {'filename': respfile, 'date': UTCDateTime(time),
                        'units': 'DIS', 'network': network, 'station': station,
                        'channel': channel, 'location': '*'}
            try:
                # Attempt to evaluate the response for this information, if not
                # then this is not the correct response info!
                freq_resp, freqs = evalresp(delta, 100, seedresp['filename'],
                                            seedresp['date'],
                                            units=seedresp['units'], freq=True,
                                            network=seedresp['network'],
                                            station=seedresp['station'],
                                            channel=seedresp['channel'])
            except:
                print('Issues with RESP file')
                seedresp = []
                continue
        elif respfile[-3:] == 'GSE':
            PAZ, pazdate, pazstation, pazchannel, pazsensor =\
                _GSE2_PAZ_read(respfile)
            # check that the date is good!
            if pazdate >= time and pazchannel != channel and\
               pazstation != station:
                print('Issue with GSE file')
                print('date: ' + str(pazdate) + ' channel: ' + pazchannel +
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


def amp_pick_event(event, st, respdir, chans=['Z'], var_wintype=True,
                   winlen=0.9, pre_pick=0.2, pre_filt=True, lowcut=1.0,
                   highcut=20.0, corners=4, min_snr=1.0, plot=False,
                   remove_old=False):
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
        body waves.  When P-S times are not available we hard-wire a P-S
        at 0.34 x hypocentral distance.

    :type event: obspy.core.event.event.Event
    :param event: Event to pick
    :type st: obspy.core.stream.Stream
    :param st: Stream associated with event
    :type respdir: str
    :param respdir: Path to the response information directory
    :type chans: list
    :param chans: List of the channels to pick on, defaults to ['Z'] - should \
        just be the orientations, e.g. Z,1,2,N,E
    :type var_wintype: bool
    :param var_wintype: If True, the winlen will be \
        multiplied by the P-S time if both P and S picks are \
        available, otherwise it will be multiplied by the \
        hypocentral distance*0.34 - derived using a p-s ratio of \
        1.68 and S-velocity of 1.5km/s to give a large window, \
        defaults to True
    :type winlen: float
    :param winlen: Length of window, see above parameter, if var_wintype is \
        False then this will be in seconds, otherwise it is the \
        multiplier to the p-s time, defaults to 0.5.
    :type pre_pick: float
    :param pre_pick: Time before the s-pick to start the cut window, defaults \
        to 0.2
    :type pre_filt: bool
    :param pre_filt: To apply a pre-filter or not, defaults to True
    :type lowcut: float
    :param lowcut: Lowcut in Hz for the pre-filter, defaults to 1.0
    :type highcut: float
    :param highcut: Highcut in Hz for the pre-filter, defaults to 20.0
    :type corners: int
    :param corners: Number of corners to use in the pre-filter
    :type min_snr: float
    :param min_snr: Minimum signal-to-noise ratio to allow a pick - see note \
        below on signal-to-noise ratio calculation.
    :type plot: bool
    :param plot: Turn plotting on or off.
    :type remove_old: bool
    :param remove_old:
        If True, will remove old amplitude picks from event and overwrite
        with new picks. Defaults to False.

    :returns: Picked event
    :rtype: :class:`obspy.core.event.Event`

    .. Note:: Signal-to-noise ratio is calculated using the filtered data by \
        dividing the maximum amplitude in the signal window (pick window) \
        by the normalized noise amplitude (taken from the whole window \
        supplied).
    """
    # Hardwire a p-s multiplier of hypocentral distance based on p-s ratio of
    # 1.68 and an S-velocity 0f 1.5km/s, deliberately chosen to be quite slow
    ps_multiplier = 0.34
    # Convert these picks into a lists
    stations = []  # List of stations
    channels = []  # List of channels
    picktimes = []  # List of pick times
    picktypes = []  # List of pick types
    distances = []  # List of hypocentral distances
    picks_out = []
    if remove_old and event.amplitudes:
        for amp in event.amplitudes:
            # Find the pick and remove it too
            pick = [p for p in event.picks if p.resource_id == amp.pick_id][0]
            event.picks.remove(pick)
            event.amplitudes.remove(amp)
    for pick in event.picks:
        if pick.phase_hint in ['P', 'S']:
            picks_out.append(pick)  # Need to be able to remove this if there
            # isn't data for a station!
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code)
            picktimes.append(pick.time)
            picktypes.append(pick.phase_hint)
            arrival = [arrival for arrival in event.origins[0].arrivals
                       if arrival.pick_id == pick.resource_id][0]
            distances.append(arrival.distance)
    st.merge()  # merge the data, just in case!
    # For each station cut the window
    uniq_stas = list(set(stations))
    del(arrival)
    for sta in uniq_stas:
        for chan in chans:
            print('Working on ' + sta + ' ' + chan)
            tr = st.select(station=sta, channel='*' + chan)
            if not tr:
                # Remove picks from file
                # picks_out=[picks_out[i] for i in range(len(picks))\
                # if picks_out[i].station+picks_out[i].channel != \
                # sta+chan]
                warnings.warn('There is no station and channel match in the '
                              'wavefile!')
                continue
            else:
                tr = tr[0]
            # Apply the pre-filter
            if pre_filt:
                try:
                    tr.split().detrend('simple').merge(fill_value=0)
                except:
                    print('Some issue splitting this one')
                    dummy = tr.split()
                    dummy.detrend('simple')
                    tr = dummy.merge(fill_value=0)
                try:
                    tr.filter('bandpass', freqmin=lowcut, freqmax=highcut,
                              corners=corners)
                except NotImplementedError:
                    print('For some reason trace is not continuous:')
                    print(tr)
                    continue
            # Find the response information
            resp_info = _find_resp(tr.stats.station, tr.stats.channel,
                                   tr.stats.network, tr.stats.starttime,
                                   tr.stats.delta, respdir)
            PAZ = []
            seedresp = []
            if resp_info and 'gain' in resp_info:
                PAZ = resp_info
            elif resp_info:
                seedresp = resp_info
            # Simulate a Wood Anderson Seismograph
            if PAZ and len(tr.data) > 10:
                # Set ten data points to be the minimum to pass
                tr = _sim_WA(tr, PAZ, None, 10)
            elif seedresp and len(tr.data) > 10:
                tr = _sim_WA(tr, None, seedresp, 10)
            elif len(tr.data) > 10:
                warnings.warn('No PAZ for ' + tr.stats.station + ' ' +
                              tr.stats.channel + ' at time: ' +
                              str(tr.stats.starttime))
                continue
            noise = tr.copy()  # Copy the data to use for noise calculation
            sta_picks = [i for i in range(len(stations))
                         if stations[i] == sta]
            pick_id = event.picks[sta_picks[0]].resource_id
            arrival = [arrival for arrival in event.origins[0].arrivals
                       if arrival.pick_id == pick_id][0]
            hypo_dist = arrival.distance
            if var_wintype and hypo_dist:
                if 'S' in [picktypes[i] for i in sta_picks] and\
                   'P' in [picktypes[i] for i in sta_picks]:
                    # If there is an S-pick we can use this :D
                    S_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'S']
                    S_pick = min(S_pick)
                    P_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'P']
                    P_pick = min(P_pick)
                    try:
                        tr.trim(starttime=S_pick - pre_pick,
                                endtime=S_pick + (S_pick - P_pick) * winlen)
                    except ValueError:
                        continue
                elif 'S' in [picktypes[i] for i in sta_picks]:
                    S_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'S']
                    S_pick = min(S_pick)
                    P_modelled = S_pick - hypo_dist * ps_multiplier
                    try:
                        tr.trim(starttime=S_pick - pre_pick,
                                endtime=S_pick + (S_pick - P_modelled) *
                                winlen)
                    except ValueError:
                        continue
                else:
                    # In this case we only have a P pick
                    P_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'P']
                    P_pick = min(P_pick)
                    S_modelled = P_pick + hypo_dist * ps_multiplier
                    try:
                        tr.trim(starttime=S_modelled - pre_pick,
                                endtime=S_modelled + (S_modelled - P_pick) *
                                winlen)
                    except ValueError:
                        continue
                # Work out the window length based on p-s time or distance
            elif 'S' in [picktypes[i] for i in sta_picks]:
                # If the window is fixed we still need to find the start time,
                # which can be based either on the S-pick (this elif), or
                # on the hypocentral distance and the P-pick

                # Take the minimum S-pick time if more than one S-pick is
                # available
                S_pick = [picktimes[i] for i in sta_picks
                          if picktypes[i] == 'S']
                S_pick = min(S_pick)
                try:
                    tr.trim(starttime=S_pick - pre_pick,
                            endtime=S_pick + winlen)
                except ValueError:
                    continue
            else:
                # In this case, there is no S-pick and the window length is
                # fixed we need to calculate an expected S_pick based on the
                # hypocentral distance, this will be quite hand-wavey as we
                # are not using any kind of velocity model.
                P_pick = [picktimes[i] for i in sta_picks
                          if picktypes[i] == 'P']
                P_pick = min(P_pick)
                hypo_dist = [distances[i] for i in sta_picks
                             if picktypes[i] == 'P'][0]
                S_modelled = P_pick + hypo_dist * ps_multiplier
                try:
                    tr.trim(starttime=S_modelled - pre_pick,
                            endtime=S_modelled + winlen)
                except ValueError:
                    continue
            if len(tr.data) <= 10:
                # Should remove the P and S picks if len(tr.data)==0
                warnings.warn('No data found for: ' + tr.stats.station)
                # print 'No data in miniseed file for '+tr.stats.station+\
                # ' removing picks'
                # picks_out=[picks_out[i] for i in range(len(picks_out))\
                # if i not in sta_picks]
                continue
            # Get the amplitude
            amplitude, period, delay = _max_p2t(tr.data, tr.stats.delta)
            # Calculate the normalized noise amplitude
            noise_amplitude = np.sqrt(np.mean(np.square(noise.data)))
            if amplitude == 0.0:
                continue
            if amplitude / noise_amplitude < min_snr:
                print('Signal to noise ratio of %s is below threshold.' %
                      (amplitude / noise_amplitude))
                continue
            if plot:
                plt.plot(np.arange(len(tr.data)), tr.data, 'k')
                plt.scatter(tr.stats.sampling_rate * delay, amplitude / 2)
                plt.scatter(tr.stats.sampling_rate * (delay + period),
                            -amplitude / 2)
                plt.show()
            print('Amplitude picked: ' + str(amplitude))
            print('Signal-to-noise ratio is: %s' %
                  (amplitude / noise_amplitude))
            # Note, amplitude should be in meters at the moment!
            # Remove the pre-filter response
            if pre_filt:
                # Generate poles and zeros for the filter we used earlier: this
                # is how the filter is designed in the convenience methods of
                # filtering in obspy.
                z, p, k = iirfilter(corners, [lowcut / (0.5 *
                                                        tr.stats.
                                                        sampling_rate),
                                              highcut / (0.5 *
                                                         tr.stats.
                                                         sampling_rate)],
                                    btype='band', ftype='butter', output='zpk')
                filt_paz = {'poles': list(p),
                            'zeros': list(z),
                            'gain': k,
                            'sensitivity': 1.0}
                amplitude /= (paz_2_amplitude_value_of_freq_resp(filt_paz,
                                                                 1 / period) *
                              filt_paz['sensitivity'])
            # Convert amplitude to mm
            if PAZ:  # Divide by Gain to get to nm (returns pm? 10^-12)
                # amplitude *=PAZ['gain']
                amplitude /= 1000
            if seedresp:  # Seedresp method returns mm
                amplitude *= 1000000
            # Write out the half amplitude, approximately the peak amplitude as
            # used directly in magnitude calculations
            # Page 343 of Seisan manual:
            #   Amplitude (Zero-Peak) in units of nm, nm/s, nm/s^2 or counts
            amplitude *= 0.5
            # Append an amplitude reading to the event
            _waveform_id = WaveformStreamID(station_code=tr.stats.station,
                                            channel_code=tr.stats.channel,
                                            network_code=tr.stats.network)
            pick_ind = len(event.picks)
            event.picks.append(Pick(waveform_id=_waveform_id,
                                    phase_hint='IAML',
                                    polarity='undecidable',
                                    time=tr.stats.starttime + delay,
                                    evaluation_mode='automatic'))
            event.amplitudes.append(Amplitude(generic_amplitude=amplitude /
                                              1e9, period=period,
                                              pick_id=event.
                                              picks[pick_ind].resource_id,
                                              waveform_id=event.
                                              picks[pick_ind].waveform_id,
                                              unit='m',
                                              magnitude_hint='ML',
                                              type='AML',
                                              category='point'))
    return event


def amp_pick_sfile(sfile, datapath, respdir, chans=['Z'], var_wintype=True,
                   winlen=0.9, pre_pick=0.2, pre_filt=True, lowcut=1.0,
                   highcut=20.0, corners=4, min_snr=1.0, plot=False,
                   remove_old=False):
    """
    Function to pick amplitudes for local magnitudes from NORDIC s-files.

    Reads information from a SEISAN s-file, load the data and the \
    picks, cut the data for the channels given around the S-window, simulate \
    a Wood Anderson seismometer, then pick the maximum peak-to-trough \
    amplitude.

    Output will be put into a mag_calc.out file which will be in full S-file \
    format and can be copied to a REA database.

    See docs for :func:`eqcorrscan.utils.mag_calc.amp_pick_event` for methods
    used here for stabilisation.

    :type sfile: str
    :param sfile: Path to NORDIC format s-file
    :type datapath: str
    :param datapath: Path to the waveform files - usually the path to the WAV \
        directory
    :type respdir: str
    :param respdir: Path to the response information directory
    :type chans: list
    :param chans: List of the channels to pick on, defaults to ['Z'] - should \
        just be the orientations, e.g. Z,1,2,N,E
    :type var_wintype: bool
    :param var_wintype: If True, the winlen will be \
        multiplied by the P-S time if both P and S picks are \
        available, otherwise it will be multiplied by the \
        hypocentral distance*0.34 - derived using a p-s ratio of \
        1.68 and S-velocity of 1.5km/s to give a large window, \
        defaults to True
    :type winlen: float
    :param winlen: Length of window, see above parameter, if var_wintype is \
        False then this will be in seconds, otherwise it is the \
        multiplier to the p-s time, defaults to 0.5.
    :type pre_pick: float
    :param pre_pick: Time before the s-pick to start the cut window, defaults \
        to 0.2
    :type pre_filt: bool
    :param pre_filt: To apply a pre-filter or not, defaults to True
    :type lowcut: float
    :param lowcut: Lowcut in Hz for the pre-filter, defaults to 1.0
    :type highcut: float
    :param highcut: Highcut in Hz for the pre-filter, defaults to 20.0
    :type corners: int
    :param corners: Number of corners to use in the pre-filter
    :type min_snr: float
    :param min_snr: Minimum signal-to-noise ratio to allow a pick - see note \
        in amp_pick_event on signal-to-noise ratio calculation.
    :type plot: bool
    :param plot: Turn plotting on or off.
    :type remove_old: bool
    :param remove_old:
        If True, will remove old amplitude picks from event and overwrite with
        new picks. Defaults to False.

    :returns: Picked event
    :rtype: :class:`obspy.core.event.event.Event`
    """
    # First we need to work out what stations have what picks
    event = sfile_util.readpicks(sfile)
    # Read in waveforms
    try:
        stream = read(os.path.join(datapath,
                                   sfile_util.readwavename(sfile)[0]))
    except IOError:
        stream = read(os.path.join(datapath,
                                   str(event.origins[0].time.year),
                                   str(event.origins[0].time.month).zfill(2),
                                   sfile_util.readwavename(sfile)[0]))
    if len(sfile_util.readwavename(sfile)) > 1:
        for wavfile in sfile_util.readwavename(sfile):
            try:
                stream += read(os.path.join(datapath, wavfile))
            except IOError:
                stream += read(os.path.join(datapath,
                                            str(event.origins[0].time.year),
                                            str(event.origins[0].time.month).
                                            zfill(2),
                                            wavfile))
    stream.merge()  # merge the data, just in case!
    event_picked = amp_pick_event(event=event, st=stream, respdir=respdir,
                                  chans=chans, var_wintype=var_wintype,
                                  winlen=winlen, pre_pick=pre_pick,
                                  pre_filt=pre_filt, lowcut=lowcut,
                                  highcut=highcut, corners=corners,
                                  min_snr=min_snr, plot=plot,
                                  remove_old=remove_old)
    new_sfile = sfile_util.eventtosfile(event=event_picked, userID=str('EQCO'),
                                        evtype=str('L'), outdir=str('.'),
                                        wavefiles=sfile_util.
                                        readwavename(sfile))
    shutil.move(new_sfile, 'mag_calc.out')
    return event_picked


def SVD_moments(U, s, V, stachans, event_list, n_SVs=4):
    """
    Calculate relative moments/amplitudes using singular-value decomposition.

    Convert basis vectors calculated by singular value \
    decomposition (see the SVD functions in clustering) into relative \
    moments.

    For more information see the paper by \
    `Rubinstein & Ellsworth (2010).
    <http://www.bssaonline.org/content/100/5A/1952.short>`_

    :type U: list
    :param U:
        List of the :class:`numpy.ndarray` input basis vectors from the SVD,
        one array for each channel used.
    :type s: list
    :param s:
        List of the :class:`numpy.ndarray` of singular values, one array for
        each channel.
    :type V: list
    :param V:
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
    :type n_SVs: int
    :param n_SVs: Number of singular values to use, defaults to 4.

    :returns: M, array of relative moments
    :rtype: :class:`numpy.ndarray`
    :returns: events_out, list of events that relate to M (in order), \
        does not include the magnitude information in the events, see note.
    :rtype: :class:`obspy.core.event.event.Event`

    .. note:: M is an array of relative moments (or amplitudes), these cannot
        be directly compared to true moments without calibration.

    .. rubric:: Example

    >>> from eqcorrscan.utils.mag_calc import SVD_moments
    >>> from obspy import read
    >>> import glob
    >>> import os
    >>> from eqcorrscan.utils.clustering import SVD
    >>> import numpy as np
    >>> # Do the set-up
    >>> testing_path = 'eqcorrscan/tests/test_data/similar_events'
    >>> stream_files = glob.glob(os.path.join(testing_path, '*'))
    >>> stream_list = [read(stream_file) for stream_file in stream_files]
    >>> event_list = []
    >>> for i, stream in enumerate(stream_list):
    ...     st_list = []
    ...     for tr in stream:
    ...         if (tr.stats.station, tr.stats.channel) not in
    ...           [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]:
    ...             stream.remove(tr)
    ...             continue
    ...         tr.detrend('simple')
    ...         tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
    ...         tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
    ...         st_list.append(i)
    ...     event_list.append(st_list) # doctest: +SKIP
    >>> event_list = np.asarray(event_list).T.tolist()
    >>> SVec, SVal, U, stachans = SVD(stream_list=stream_list) # doctest: +SKIP
    ['GCSZ.EHZ', 'WV04.SHZ', 'WHAT2.SH1']
    >>> M, events_out = SVD_moments(U=U, s=SVal, V=SVec, stachans=stachans,
    ...                             event_list=event_list) # doctest: +SKIP

    """
    # Define maximum number of events, will be the width of K
    K_width = max([max(ev_list) for ev_list in event_list]) + 1
    # Sometimes the randomisation generates a singular matrix - rather than
    # attempting to regulerize this matrix I propose undertaking the
    # randomisation step a further time
    if len(stachans) == 1:
        print('Only provided data from one station-channel - '
              'will not try to invert')
        return U[0][:, 0], event_list[0]
    for i, stachan in enumerate(stachans):
        k = []  # Small kernel matrix for one station - channel
        # Copy the relevant vectors so as not to destroy them
        U_working = copy.deepcopy(U[i])
        V_working = copy.deepcopy(V[i])
        s_working = copy.deepcopy(s[i])
        ev_list = event_list[i]
        if len(ev_list) > len(V_working):
            print('V is : ' + str(len(V_working)))
            f_dump = open('mag_calc_V_working.pkl', 'wb')
            pickle.dump(V_working, f_dump)
            f_dump.close()
            raise IOError('More events than represented in V')
        # Set all non-important singular values to zero
        s_working[n_SVs:len(s_working)] = 0
        s_working = np.diag(s_working)
        # Convert to numpy matrices
        U_working = np.matrix(U_working)
        V_working = np.matrix(V_working)
        s_working = np.matrix(s_working)

        SVD_weights = U_working[:, 0]
        # If all the weights are negative take the abs
        if np.all(SVD_weights < 0):
            warnings.warn('All weights are negative - flipping them')
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
            min_weight = min(_weights)
            max_weight = max(_weights)
            min_index = _indeces[np.argmin(_weights)]
            max_index = _indeces[np.argmax(_weights)]
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
    print("\nCreated Kernel matrix: ")
    del row
    print('\n'.join([''.join([str(round(float(item), 3)).ljust(6)
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

    return M, events_out


def pick_db(indir, outdir, calpath, startdate, enddate, wavepath=None):
    """
    Wrapper to loop through a SEISAN database and make a lot of magnitude \
    picks.

    :type indir: str
    :param indir: Path to the seisan REA directory (not including yyyy/mm)
    :type outdir: str
    :param outdir: Path to output seisan REA directory (not including yyyy/mm)
    :type calpath: str
    :param calpath: Path to the directory containing the response files
    :type startdate: datetime.datetime
    :param startdate: Date to start looking for S-files
    :type enddate: datetime.datetime
    :param enddate: Date to stop looking for S-files
    :type wavepath: str
    :param wavepath: Path to the seisan WAV directory (not including yyyy/mm)
    """
    kdays = ((enddate + dt.timedelta(1)) - startdate).days
    for i in range(kdays):
        day = startdate + dt.timedelta(i)
        print('Working on ' + str(day))
        sfiles = glob.glob(os.path.join(indir, str(day.year),
                                        str(day.month).zfill(2),
                                        str(day.day).zfill(2) + '-*L.S' +
                                        str(day.year) +
                                        str(day.month).zfill(2)))
        datetimes = [dt.datetime.strptime(os.path.split(sfiles[i])[-1],
                                          '%d-%H%M-%SL.S%Y%m')
                     for i in range(len(sfiles))]
        sfiles = [sfiles[i] for i in range(len(sfiles))
                  if startdate < datetimes[i] < enddate]
        if not wavepath:
            wavedir = os.path.join(os.path.split(indir), 'WAV',
                                   os.path.split(indir)[-1],
                                   str(day.year), str(day.month).zfill(2))
        else:
            wavedir = os.path.join(wavepath, str(day.year),
                                   str(day.month).zfill(2))
        sfiles.sort()
        for sfile in sfiles:
            # Make the picks!
            print('\tWorking on Sfile: ' + sfile)
            event = amp_pick_sfile(sfile, wavedir, calpath)
            del event
            # Copy the mag_calc.out file to the correct place
            shutil.copyfile('mag_calc.out',
                            os.path.join(outdir, str(day.year),
                                         str(day.month).zfill(2),
                                         os.path.split(sfile)[-1]))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
