#!/usr/bin/python
"""
Functions to simulate Wood Anderson traces, pick maximum peak-to-peak \
amplitudes write these amplitudes and periods to SEISAN s-files and to \
calculate magnitudes from this and the informaiton within SEISAN s-files.

Written as part of the EQcorrscan package by Calum Chamberlain - first \
written to impliment magnitudes for the 2015 Wanaka aftershock sequence, \
written up by Warren-Smith [2014/15].

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
import warnings


def dist_calc(loc1, loc2):
    """
    Function to calcualte the distance in km between two points, uses the \
    flat Earth approximation.

    :type loc1: tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: tuple
    :param loc2: Tuple of lat, lon, depth (in decimal degrees and km)

    :returns: float, Distance between points.
    """
    R = 6371.009  # Radius of the Earth in km
    dlat = np.radians(abs(loc1[0] - loc2[0]))
    dlong = np.radians(abs(loc1[1] - loc2[1]))
    ddepth = abs(loc1[2] - loc2[2])
    mean_lat = np.radians((loc1[0] + loc2[0]) / 2)
    dist = R * np.sqrt(dlat ** 2 + (np.cos(mean_lat) * dlong) ** 2)
    dist = np.sqrt(dist ** 2 + ddepth ** 2)
    return dist


def _sim_WA(trace, PAZ, seedresp, water_level):
    """
    Function to remove the insturment response from a trace and return a \
    de-meaned, de-trended, Wood Anderson simulated trace in it's place.

    Works in-place on data and will destroy your original data, copy the \
    trace before giving it to this function!

    :type trace: obspy.Trace
    :param trace: A standard obspy trace, generally should be given without
                    pre-filtering, if given with pre-filtering for use with
                    amplitude determiniation for magnitudes you will need to
                    worry about how you cope with the response of this filter
                    yourself.
    :type PAZ: dict
    :param PAZ: Dictionary containing lists of poles and zeros, the gain and
                the sensitivity.
    :type water_level: int
    :param water_level: Water level for the simulation.

    :returns: obspy.Trace
    """
    # Note Wood anderson sensitivity is 2080 as per Uhrhammer & Collins 1990
    PAZ_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
              'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    from obspy.signal import seisSim
    # De-trend data
    trace.detrend('simple')
    # Simulate Wood Anderson
    if PAZ:
        trace.data = seisSim(trace.data, trace.stats.sampling_rate,
                             paz_remove=PAZ, paz_simulate=PAZ_WA,
                             water_level=water_level, remove_sensitivity=True)
    elif seedresp:
        trace.data = seisSim(trace.data, trace.stats.sampling_rate,
                             paz_remove=None, paz_simulate=PAZ_WA,
                             water_level=water_level, seedresp=seedresp)
    else:
        UserWarning('No response given to remove, will just simulate WA')
        trace.data = seisSim(trace.data, trace.stats.samplng_rate,
                             paz_remove=None, water_level=water_level)
    return trace


def _max_p2t(data, delta):
    """
    Function to find the maximum peak-to-trough amplitude and period of this \
    amplitude.  Originally designed to be used to calculate magnitudes (by \
    taking half of the peak-to-trough amplitude as the peak amplitude).

    :type data: ndarray
    :param data: waveform trace to find the peak-to-trough in.
    :type delta: float
    :param delta: Sampling interval in seconds

    :returns: tuple of (amplitude, period, time) with amplitude in the same \
        scale as given in the input data, and period in seconds, and time in \
        seconds from the start of the data window.
    """
    import matplotlib.pyplot as plt
    debug_plot = False
    turning_points = []  # A list of tuples of (amplitude, sample)
    for i in range(1, len(data) - 1):
        if (data[i] < data[i-1] and data[i] < data[i+1]) or\
           (data[i] > data[i-1] and data[i] > data[i+1]):
            turning_points.append((data[i], i))
    if len(turning_points) >= 1:
        amplitudes = np.empty([len(turning_points)-1],)
        half_periods = np.empty([len(turning_points)-1],)
    else:
        plt.plot(data)
        plt.show()
        print 'Turning points has length: '+str(len(turning_points)) +\
            ' data have length: '+str(len(data))
        return (0.0, 0.0, 0.0)
    for i in range(1, len(turning_points)):
        half_periods[i-1] = (delta * (turning_points[i][1] -
                                      turning_points[i-1][1]))
        amplitudes[i-1] = np.abs(turning_points[i][0]-turning_points[i-1][0])
    amplitude = np.max(amplitudes)
    period = 2 * half_periods[np.argmax(amplitudes)]
    if debug_plot:
        plt.plot(data, 'k')
        plt.plot([turning_points[np.argmax(amplitudes)][1],
                  turning_points[np.argmax(amplitudes)-1][1]],
                 [turning_points[np.argmax(amplitudes)][0],
                  turning_points[np.argmax(amplitudes)-1][0]], 'r')
        plt.show()
    return (amplitude, period, delta*turning_points[np.argmax(amplitudes)][1])


def _GSE2_PAZ_read(GSEfile):
    """
    Function to read the instrument response information from a GSE Poles and \
    Zeros file as generated by the SEISAN program RESP.

    Format must be CAL2, not coded for any other format at the moment, \
    contact the authors to add others in.

    :type GSEfile: string
    :param GSEfile: Path to GSE file

    :returns: Dict of poles, zeros, gain and sensitivity
    """
    import datetime as dt
    f = open(GSEfile)
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
        poles.append(complex(float(pole.split()[0]), float(pole.split()[1])))
    zeros = []
    for i in range(kzeros):
        zero = f.readline()
        zeros.append(complex(float(zero.split()[0]), float(zero.split()[1])))
    # Have Poles and Zeros, but need Gain and Sensitivity
    # Gain should be in the DIG2 line:
    for line in f:
        if line[0:4] == 'DIG2':
            sensitivity = float(line.split()[2])  # measured in counts/muVolt
    f.close()
    PAZ = {'poles': poles, 'zeros': zeros, 'gain': gain,
           'sensitivity': sensitivity}
    return PAZ, date, station, channel, sensor


def _find_resp(station, channel, network, time, delta, directory):
    """
    Helper function to find the response information for a given station and \
    channel at a given time and return a dictionary of poles and zeros, gain \
    and sensitivity.

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

    :returns: dict, response information
    """
    import glob
    from obspy.signal.invsim import evalresp
    from obspy import UTCDateTime
    possible_respfiles = glob.glob(directory+'/RESP.'+network+'.'+station +
                                   '.*.' + channel)  # GeoNet RESP naming
    possible_respfiles += glob.glob(directory+'/RESP.'+network+'.'+channel +
                                    '.' + station)  # RDseed RESP naming
    possible_respfiles += glob.glob(directory+'/RESP.'+station+'.'+network)
    # WIZARD resp naming
    # GSE format, station needs to be 5 charectars padded with _, channel is 4
    # characters padded with _
    possible_respfiles += glob.glob(directory+'/'+station.ljust(5, '_') +
                                    channel[0:len(channel)-1].ljust(3, '_') +
                                    channel[-1]+'.*_GSE')
    PAZ = []
    seedresp = []
    for respfile in possible_respfiles:
        print 'Reading response from: '+respfile
        if respfile.split('/')[-1][0:4] == 'RESP':
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
                print 'Issues with RESP file'
                seedresp = []
                continue
        elif respfile[-3:] == 'GSE':
            PAZ, pazdate, pazstation, pazchannel, pazsensor =\
                _GSE2_PAZ_read(respfile)
            # check that the date is good!
            if pazdate >= time and pazchannel != channel and\
               pazstation != station:
                print 'Issue with GSE file'
                print 'date: '+str(pazdate)+' channel: '+pazchannel +\
                    ' station: '+pazstation
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
    import itertools
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


def Amp_pick_sfile(sfile, datapath, respdir, chans=['Z'], var_wintype=True,
                   winlen=0.9, pre_pick=0.2, pre_filt=True, lowcut=1.0,
                   highcut=20.0, corners=4):
    """
    Function to read information from a SEISAN s-file, load the data and the \
    picks, cut the data for the channels given around the S-window, simulate \
    a Wood Anderson seismometer, then pick the maximum peak-to-trough \
    amplitude.

    Output will be put into a mag_calc.out file which will be in full S-file \
    format and can be copied to a REA database.

    :type sfile: string
    :type datapath: string
    :param datapath: Path to the waveform files - usually the path to the WAV \
        directory
    :type respdir: string
    :param respdir: Path to the response information directory
    :type chans: List of strings
    :param chans: List of the channels to pick on, defaults to ['Z'] - should \
        just be the orientations, e.g. Z,1,2,N,E
    :type var_wintype: bool
    :param var_wintype: If True, the winlen will be \
        multiplied by the P-S time if both P and S picks are \
        available, otherwise it will be multiplied by the \
        hypocentral distance*0.34 - dervided using a p-s ratio of \
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
    """
    # Hardwire a p-s multiplier of hypocentral distance based on p-s ratio of
    # 1.68 and an S-velocity 0f 1.5km/s, deliberately chosen to be quite slow
    ps_multiplier = 0.34
    from eqcorrscan.utils import Sfile_util
    from obspy import read
    from scipy.signal import iirfilter
    from obspy.signal.invsim import paz2AmpValueOfFreqResp
    import warnings
    # First we need to work out what stations have what picks
    event = Sfile_util.readpicks(sfile)[0]
    # Convert these picks into a lists
    stations = []  # List of stations
    channels = []  # List of channels
    picktimes = []  # List of pick times
    picktypes = []  # List of pick types
    distances = []  # List of hypocentral distances
    picks_out = []
    for pick in event.picks:
        if pick.phase_hint in ['P', 'S']:
            picks_out.append(pick)  # Need to be able to remove this if there
            # isn't data for a station!
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code)
            picktimes.append(pick.time)
            picktypes.append(pick.phase_hint)
            arrival = [arrival for arrival in event.origins[0].arrivals
                       if arrival.pick_id == pick.resource_id]
            distances.append(arrival.distance)
    # Read in waveforms
    stream = read(datapath+'/'+Sfile_util.readwavename(sfile)[0])
    if len(Sfile_util.readwavename(sfile)) > 1:
        for wavfile in Sfile_util.readwavename(sfile):
            stream += read(datapath+'/'+wavfile)
    stream.merge()  # merge the data, just in case!
    # For each station cut the window
    uniq_stas = list(set(stations))
    del arrival
    for sta in uniq_stas:
        for chan in chans:
            print 'Working on '+sta+' '+chan
            tr = stream.select(station=sta, channel='*'+chan)
            if not tr:
                # Remove picks from file
                # picks_out=[picks_out[i] for i in xrange(len(picks))\
                # if picks_out[i].station+picks_out[i].channel != \
                # sta+chan]
                warnings.warn('There is no station and channel match in the ' +
                              'wavefile!')
                break
            else:
                tr = tr[0]
            # Apply the pre-filter
            if pre_filt:
                try:
                    tr.detrend('simple')
                except:
                    dummy = tr.split()
                    dummy.detrend('simple')
                    tr = dummy.merge()[0]
                tr.filter('bandpass', freqmin=lowcut, freqmax=highcut,
                          corners=corners)
            sta_picks = [i for i in xrange(len(stations))
                         if stations[i] == sta]
            pick_id = event.picks[sta_picks[0]].resource_id
            arrival = [arrival for arrival in event.origins[0].arrivals
                       if arrival.pick_id == pick_id]
            hypo_dist = arrival.distance
            CAZ = arrival.azimuth
            if var_wintype:
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
                        tr.trim(starttime=S_pick-pre_pick,
                                endtime=S_pick+(S_pick-P_pick)*winlen)
                    except:
                        break
                elif 'S' in [picktypes[i] for i in sta_picks]:
                    S_pick = [picktimes[i] for i in sta_picks
                              if picktypes[i] == 'S']
                    S_pick = min(S_pick)
                    P_modelled = S_pick - hypo_dist * ps_multiplier
                    try:
                        tr.trim(starttime=S_pick-pre_pick,
                                endtime=S_pick + (S_pick - P_modelled) *
                                winlen)
                    except:
                        break
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
                    except:
                        break
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
                except:
                    break
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
                except:
                    break
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
                warnings.warn('No PAZ for '+tr.stats.station+' ' +
                              tr.stats.channel+' at time: ' +
                              str(tr.stats.starttime))
                continue
            if len(tr.data) <= 10:
                # Should remove the P and S picks if len(tr.data)==0
                warnings.warn('No data found for: '+tr.stats.station)
                # print 'No data in miniseed file for '+tr.stats.station+\
                # ' removing picks'
                # picks_out=[picks_out[i] for i in xrange(len(picks_out))\
                # if i not in sta_picks]
                break
            # Get the amplitude
            amplitude, period, delay = _max_p2t(tr.data, tr.stats.delta)
            if amplitude == 0.0:
                break
            print 'Amplitude picked: ' + str(amplitude)
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
                            'sensitivity':  1.0}
                amplitude /= (paz2AmpValueOfFreqResp(filt_paz, 1/period) *
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
            # Generate a PICK type object for this pick
            picks_out.append(Sfile_util.PICK(station=tr.stats.station,
                                             channel=tr.stats.channel,
                                             impulsivity=' ',
                                             phase='IAML',
                                             weight='', polarity=' ',
                                             time=tr.stats.starttime+delay,
                                             coda=999, amplitude=amplitude,
                                             peri=period, azimuth=float('NaN'),
                                             velocity=float('NaN'), AIN=999,
                                             SNR='',
                                             azimuthres=999,
                                             timeres=float('NaN'),
                                             finalweight=999,
                                             distance=hypo_dist,
                                             CAZ=CAZ))
    # Copy the header from the sfile to a new local S-file
    fin = open(sfile, 'r')
    fout = open('mag_calc.out', 'w')
    for line in fin:
        if not line[79] == '7':
            fout.write(line)
        else:
            fout.write(line)
            break
    fin.close()
    for pick in picks_out:
        fout.write(pick)
        # Note this uses the legacy pick class
    fout.close()
    # Write picks out to new s-file
    for pick in picks_out:
        print pick
    # Sfile_util.populateSfile('mag_calc.out', picks_out)
    return picks_out


def SVD_moments(U, s, V, stachans, event_list, n_SVs=4):
    """
    Function to convert basis vectors calculated by singular value \
    decomposition (see the SVD functions in clustering) into relative \
    magnitudes.

    :type U: List of np.ndarray
    :param U: List of the input basis vectors from the SVD, one array for \
        each channel used.
    :type s: List of nd.array
    :param s: List of the singular values, one array for each channel
    :type V: List of np.ndarry
    :param V: List of output basis vectors from SVD, one array per channel.
    :type stachans: List of string
    :param stachans: List of station.channel input
    :type event_list: List of list
    :param event_list: List of events for which you have data, such that \
        event_list[i] corresponds to stachans[i], U[i] etc. and \
        event_list[i][j] corresponds to event j in U[i]
    type n_SVs: int
    :param n_SVs: Number of singular values to use, defaults to 4.

    :returns: M, np.array of relative moments
    """
    import copy
    import random
    import pickle

    # Copying script from one obtained from John Townend.
    # Define maximum number of events, will be the width of K
    K_width = max([max(ev_list) for ev_list in event_list])+1
    # Sometimes the randomisation generates a singular matrix - rather than
    # attempting to regulerize this matrix I propose undertaking the
    # randomisation step a further time
    for i, stachan in enumerate(stachans):
        k = []  # Small kernel matrix for one station - channel
        # Copy the relevant vectors so as not to detroy them
        U_working = copy.deepcopy(U[i])
        V_working = copy.deepcopy(V[i])
        s_working = copy.deepcopy(s[i])
        ev_list = event_list[i]
        if len(ev_list) > len(V_working):
            print 'V is : '+str(len(V_working))
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
            _weights = zip(*list(pairs[pairsIndex]))[0]
            _indeces = zip(*list(pairs[pairsIndex]))[1]
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
                    normalised = max_weight/min_weight
                    result = float(normalised)
                else:
                    result = 0
                row.append(result)
            print row
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
    print "\nCreated Kernel matrix: "
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
    Kinv = np.array(np.linalg.inv(K.T*W*K) * K.T * W)

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
    import datetime as dt
    import glob
    import shutil

    kdays = ((enddate + dt.timedelta(1)) - startdate).days
    for i in xrange(kdays):
        day = startdate + dt.timedelta(i)
        print 'Working on ' + str(day)
        sfiles = glob.glob(indir + '/' + str(day.year) + '/' +
                           str(day.month).zfill(2) + '/' +
                           str(day.day).zfill(2) + '-*L.S' + str(day.year) +
                           str(day.month).zfill(2))
        datetimes = [dt.datetime.strptime(sfiles[i].split('/')[-1],
                                          '%d-%H%M-%SL.S%Y%m')
                     for i in xrange(len(sfiles))]
        sfiles = [sfiles[i] for i in xrange(len(sfiles))
                  if datetimes[i] > startdate and
                  datetimes[i] < enddate]
        if not wavepath:
            wavedir = "/".join(indir.split('/')[:-2])+'/WAV/' +\
                indir.split('/')[-1]+'/'+str(day.year)+'/' +\
                str(day.month).zfill(2)
        else:
        	wavedir = wavepath+'/'+str(day.year)+'/' +\
                    str(day.month).zfill(2)
        sfiles.sort()
        for sfile in sfiles:
            # Make the picks!
            print '				Working on Sfile: '+sfile
            picks = Amp_pick_sfile(sfile, wavedir, calpath)
            del picks
            # Copy the mag_calc.out file to the correct place
            shutil.copyfile('mag_calc.out', outdir+'/'+str(day.year)+'/' +
                            str(day.month).zfill(2)+'/'+sfile.split('/')[-1])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
