#!/usr/bin/python
"""
Functions to simulate Wood Anderson traces, pick maximum peak-to-peak amplitudes
write these amplitudes and periods to SEISAN s-files and to calculate magnitudes
from this and the informaiton within SEISAN s-files.

Written as part of the EQcorrscan package by Calum Chamberlain - first written
to impliment magnitudes for the 2015 Wanaka aftershock sequence, written up
by Warren-Smith [2014/15].

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

def dist_calc(loc1, loc2):
    """
    Function to calcualte the distance in km between two points, uses the flat
    Earth approximation

    :type loc1: Tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: Tuple
    :param loc2: Tuple of lat, lon, depth (in decimal degrees and km)
    """
    R=6371.009  # Radius of the Earth in km
    dlat=np.radians(abs(loc1[0]-loc2[0]))
    dlong=np.radians(abs(loc1[1]-loc2[1]))
    ddepth=abs(loc1[2]-loc2[2])
    mean_lat=np.radians((loc1[0]+loc2[0])/2)
    dist=R*np.sqrt(dlat**2+(np.cos(mean_lat)*dlong)**2)
    dist=np.sqrt(dist**2+ddepth**2)
    return dist

def _sim_WA(trace, PAZ, seedresp, water_level):
    """
    Function to remove the insturment response from a trace and return a
    de-meaned, de-trended, Wood Anderson simulated trace in it's place.

    Works in-place on data and will destroy your original data, copy the
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
    PAZ_WA={'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
    from obspy.signal import seisSim
    # De-trend data
    trace.detrend('simple')
    # Simulate Wood Anderson
    if PAZ:
        trace.data=seisSim(trace.data, trace.stats.sampling_rate, paz_remove=PAZ,\
                           paz_simulate=PAZ_WA, water_level=water_level,\
                           remove_sensitivity=True)
    elif seedresp:
        trace.data=seisSim(trace.data, trace.stats.sampling_rate, paz_remove=None,\
                           paz_simulate=PAZ_WA, water_level=water_level,\
                           seedresp=seedresp)
    else:
        UserWarning('No response given to remove, will just simulate WA')
        trace.data=seisSim(trace.data, trace.stats.samplng_rate, paz_remove=None,\
                            water_level=water_level)
    return trace

def _max_p2t(data, delta):
    """
    Function to find the maximum peak to trough amplitude and period of this
    amplitude.  Originally designed to be used to calculate magnitudes (by
    taking half of the peak-to-trough amplitude as the peak amplitude).

    :type data: ndarray
    :param data: waveform trace to find the peak-to-trough in.
    :type delta: float
    :param delta: Sampling interval in seconds

    :returns: tuple of (amplitude, period, time) with amplitude in the same scale as
            given in the input data, and period in seconds, and time in seconds
            from the start of the data window.
    """
    import matplotlib.pyplot as plt
    debug_plot=False
    turning_points=[] # A list of tuples of (amplitude, sample)
    for i in xrange(1,len(data)-1):
        if (data[i] < data[i-1] and data[i] < data[i+1]) or\
           (data[i] > data[i-1] and data[i] > data[i+1]):
            turning_points.append((data[i], i))
    if len(turning_points) >=1:
        amplitudes=np.empty([len(turning_points)-1],)
        half_periods= np.empty([len(turning_points)-1],)
    else:
        plt.plot(data)
        plt.show()
        print 'Turning points has length: '+str(len(turning_points))+\
                         ' data have length: '+str(len(data))
        return (0.0, 0.0, 0.0)
    for i in xrange(1,len(turning_points)):
        half_periods[i-1]=(delta*(turning_points[i][1]-turning_points[i-1][1]))
        amplitudes[i-1]=np.abs(turning_points[i][0]-turning_points[i-1][0])
    amplitude=np.max(amplitudes)
    period=2*half_periods[np.argmax(amplitudes)]
    if debug_plot:
        plt.plot(data,'k')
        plt.plot([turning_points[np.argmax(amplitudes)][1],\
                  turning_points[np.argmax(amplitudes)-1][1]],\
                 [turning_points[np.argmax(amplitudes)][0],\
                  turning_points[np.argmax(amplitudes)-1][0]], 'r')
        plt.show()
    return (amplitude, period, delta*turning_points[np.argmax(amplitudes)][1])

def _GSE2_PAZ_read(GSEfile):
    """
    Function to read the instrument response information from a GSE Poles and
    Zeros file as generated by the SEISAN program RESP.

    Format must be CAL2, not coded for any other format at the moment, contact
    the author to add others in.

    :type GSEfile: Str
    :param GSEfile: Path to GSE file

    :returns: Dict of poles, zeros, gain and sensitivity
    """
    import datetime as dt
    f=open(GSEfile)
    # First line should start with CAL2
    header=f.readline()
    if not header[0:4] == 'CAL2':
        raise IOError('Unknown format for GSE file, only coded for CAL2')
    station=header.split()[1]
    channel=header.split()[2]
    sensor=header.split()[3]
    date=dt.datetime.strptime(header.split()[7], '%Y/%m/%d')
    header=f.readline()
    if not header[0:4] == 'PAZ2':
        raise IOError('Unknown format for GSE file, only coded for PAZ2')
    gain=float(header.split()[3]) # Measured in nm/counts
    kpoles=int(header.split()[4])
    kzeros=int(header.split()[5])
    poles=[]
    for i in xrange(kpoles):
        pole=f.readline()
        poles.append(complex(float(pole.split()[0]),float(pole.split()[1])))
    zeros=[]
    for i in xrange(kzeros):
        zero=f.readline()
        zeros.append(complex(float(zero.split()[0]),float(zero.split()[1])))
    # Have Poles and Zeros, but need Gain and Sensitivity
    # Gain should be in the DIG2 line:
    for line in f:
        if line[0:4] == 'DIG2':
            sensitivity=float(line.split()[2]) # measured in counts/muVolt
    f.close()
    PAZ={'poles': poles, 'zeros': zeros, 'gain': gain,
         'sensitivity': sensitivity}
    return PAZ, date, station, channel, sensor

def _find_resp(station, channel, network, time, delta, directory):
    """
    Helper function to find the response information for a given station and
    channel at a given time and return a dictionary of poles and zeros, gain
    and sensitivity.

    :type station: String
    :param station: Station name (as in the response files)
    :type channel: String
    :param channel: Channel name (as in the response files)
    :type network: String
    :param network: Network to scan for, can be a wildcard
    :type time: datetime.datetime
    :param time: Date-time to look for repsonse information
    :type delta: float
    :param delta: Sample interval in seconds
    :type directory: String
    :param directory: Directory to scan for response information

    :returns: Dictionary
    """
    import glob
    from obspy.signal.invsim import evalresp
    from obspy import UTCDateTime
    possible_respfiles=glob.glob(directory+'/RESP.'+network+'.'+station+'.*.'+\
                                 channel) # GeoNet RESP naming
    possible_respfiles+=glob.glob(directory+'/RESP.'+network+'.'+channel+'.'+\
                                 station) # RDseed RESP naming
    possible_respfiles+=glob.glob(directory+'/RESP.'+station+'.'+network)
                                    # WIZARD resp naming
    # GSE format, station needs to be 5 charectars padded with _, channel is 4
    # characters padded with _
    possible_respfiles+=glob.glob(directory+'/'+station.ljust(5,'_')+\
                                  channel[0:len(channel)-1].ljust(3,'_')+\
                                  channel[-1]+'.*_GSE')
    PAZ=[]
    seedresp=[]
    for respfile in possible_respfiles:
        print 'Reading response from: '+respfile
        if respfile.split('/')[-1][0:4]=='RESP':
            # Read from a resp file
            seedresp={'filename': respfile, 'date': UTCDateTime(time),
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
                seedresp=[]
                continue
        elif respfile[-3:]=='GSE':
            PAZ, pazdate, pazstation, pazchannel, pazsensor =\
                    _GSE2_PAZ_read(respfile)
            # check that the date is good!
            if pazdate >= time and pazchannel != channel and\
               pazstation != station:
                print 'Issue with GSE file'
                print 'date: '+str(pazdate)+' channel: '+pazchannel+\
                        ' station: '+pazstation
                PAZ=[]
        else:
            continue
        # Check that PAZ are for the correct station, channel and date
        if PAZ or seedresp:
            break
    if PAZ:
        return PAZ
    elif seedresp:
        return seedresp

def Amp_pick_sfile(sfile, datapath, respdir, chans=['Z'], var_wintype=True, \
                   winlen=0.9, pre_pick=0.2, pre_filt=True, lowcut=1.0,\
                   highcut=20.0, corners=4):
    """
    Function to read information from a SEISAN s-file, load the data and the
    picks, cut the data for the channels given around the S-window, simulate
    a Wood Anderson seismometer, then pick the maximum peak-to-trough
    amplitude.

    Output will be put into a mag_calc.out file which will be in full S-file
    format and can be copied to a REA database.

    :type sfile: String
    :type datapath: String
    :param datapath: Path to the waveform files - usually the path to the WAV directory
    :type respdir: String
    :param respdir: Path to the response information directory
    :type chans: List of strings
    :param chans: List of the channels to pick on, defaults to ['Z'] - should
                just be the orientations, e.g. Z,1,2,N,E
    :type var_wintype: Bool
    :param var_wintype: If True, the winlen will be
                    multiplied by the P-S time if both P and S picks are
                    available, otherwise it will be multiplied by the hypocentral
                    distance*0.34 - dervided using a p-s ratio of 1.68 and
                    S-velocity of 1.5km/s to give a large window, defaults to True
    :type winlen: Float
    :param winlen: Length of window, see above parameter, if var_wintype is False
                    Then this will be in seconds, otherwise it is the multiplier
                    to the p-s time, defaults to 0.5
    :type pre_pick: Float
    :param pre_pick: Time before the s-pick to start the cut window, defaults
                    to 0.2
    :type pre_filt: Bool
    :param pre_filt: To apply a pre-filter or not, defaults to True
    :type lowcut: Float
    :param lowcut: Lowcut in Hz for the pre-filter, defaults to 1.0
    :type highcut: Float
    :param highcut: Highcut in Hz for the pre-filter, defaults to 20.0
    :type corners: Int
    :param corners: Number of corners to use in the pre-filter
    """
    # Hardwire a p-s multiplier of hypocentral distance based on p-s ratio of
    # 1.68 and an S-velocity 0f 1.5km/s, deliberately chosen to be quite slow
    ps_multiplier=0.34
    try:
        from utils import Sfile_util
    except:
        import Sfile_util
    from obspy import read
    from scipy.signal import iirfilter
    from obspy.signal.invsim import paz2AmpValueOfFreqResp
    import warnings
    # First we need to work out what stations have what picks
    picks=Sfile_util.readpicks(sfile)
    # Convert these picks into a lists
    stations=[] # List of stations
    channels=[] # List of channels
    picktimes=[] # List of pick times
    picktypes=[] # List of pick types
    distances=[] # List of hypocentral distances
    picks_out=[]
    for pick in picks:
        if pick.phase in ['P','S']:
            picks_out.append(pick) # Need to be able to remove this if there
                                   # isn't data for a station!
            stations.append(pick.station)
            channels.append(pick.channel)
            picktimes.append(pick.time)
            picktypes.append(pick.phase)
            distances.append(pick.distance)
    # Read in waveforms
    stream=read(datapath+'/'+Sfile_util.readwavename(sfile)[0])
    if len(Sfile_util.readwavename(sfile)) > 1:
        for wavfile in Sfile_util.readwavename(sfile):
            stream+=read(datapath+'/'+wavfile)
    stream.merge() # merge the data, just in case!
    # For each station cut the window
    uniq_stas=list(set(stations))
    for sta in uniq_stas:
        for chan in chans:
            print 'Working on '+sta+' '+chan
            tr=stream.select(station=sta, channel='*'+chan)
            if not tr:
                # Remove picks from file
                # picks_out=[picks_out[i] for i in xrange(len(picks))\
                           # if picks_out[i].station+picks_out[i].channel != \
                           # sta+chan]
            	warnings.warn('There is no station and channel match in the wavefile!')
                break
            else:
                tr=tr[0]
            # Apply the pre-filter
            if pre_filt:
                try:
                    tr.detrend('simple')
                except:
                    dummy=tr.split()
                    dummy.detrend('simple')
                    tr=dummy.merge()[0]
                tr.filter('bandpass',freqmin=lowcut, freqmax=highcut,\
                             corners=corners)
            sta_picks=[i for i in xrange(len(stations)) \
                           if stations[i]==sta]
            hypo_dist=picks[sta_picks[0]].distance
            CAZ=picks[sta_picks[0]].CAZ
            if var_wintype:
                if 'S' in [picktypes[i] for i in sta_picks] and\
                   'P' in [picktypes[i] for i in sta_picks]:
                    # If there is an S-pick we can use this :D
                    S_pick=[picktimes[i] for i in sta_picks \
                            if picktypes[i]=='S']
                    S_pick=min(S_pick)
                    P_pick=[picktimes[i] for i in sta_picks \
                            if picktypes[i]=='P']
                    P_pick=min(P_pick)
                    try:
                    	tr.trim(starttime=S_pick-pre_pick, \
                               endtime=S_pick+(S_pick-P_pick)*winlen)
                    except:
                    	break
                elif 'S' in [picktypes[i] for i in sta_picks]:
                    S_pick=[picktimes[i] for i in sta_picks \
                            if picktypes[i]=='S']
                    S_pick=min(S_pick)
                    P_modelled=S_pick-hypo_dist*ps_multiplier
                    try:
                    	tr.trim(starttime=S_pick-pre_pick,\
                            endtime=S_pick+(S_pick-P_modelled)*winlen)
                    except:
                    	break
                else:
                    # In this case we only have a P pick
                    P_pick=[picktimes[i] for i in sta_picks \
                            if picktypes[i]=='P']
                    P_pick=min(P_pick)
                    S_modelled=P_pick+hypo_dist*ps_multiplier
                    try:
                    	tr.trim(starttime=S_modelled-pre_pick,\
                        	    endtime=S_modelled+(S_modelled-P_pick)*winlen)
                    except:
                    	break
                # Work out the window length based on p-s time or distance
            elif 'S' in [picktypes[i] for i in sta_picks]:
                # If the window is fixed we still need to find the start time,
                # which can be based either on the S-pick (this elif), or
                # on the hypocentral distance and the P-pick

                # Take the minimum S-pick time if more than one S-pick is available
                S_pick=[picktimes[i] for i in sta_picks \
                           if picktypes[i]=='S']
                S_pick=min(S_pick)
                try:
                	tr.trim(starttime=S_pick-pre_pick, endtime=S_pick+winlen)
                except:
                	break
            else:
                # In this case, there is no S-pick and the window length is fixed
                # We need to calculate an expected S_pick based on the hypocentral
                # distance, this will be quite hand-wavey as we are not using
                # any kind of velocity model.
                P_pick=[picktimes[i] for i in sta_picks \
                           if picktypes[i]=='P']
                P_pick=min(P_pick)
                hypo_dist=[distances[i] for i in sta_picks\
                           if picktypes[i]=='P'][0]
                S_modelled=P_pick+hypo_dist*ps_multiplier
                try:
                	tr.trim(starttime=S_modelled-pre_pick,\
                    	       endtime=S_modelled+winlen)
                except:
        	       break
            # Find the response information
            resp_info=_find_resp(tr.stats.station, tr.stats.channel,\
                           tr.stats.network, tr.stats.starttime, tr.stats.delta,\
                                 respdir)
            PAZ=[]
            seedresp=[]
            if resp_info and 'gain' in resp_info:
                PAZ=resp_info
            elif resp_info:
                seedresp=resp_info
            # Simulate a Wood Anderson Seismograph
            if PAZ and len(tr.data) > 10: # Set ten data points to be the minimum to pass
                tr=_sim_WA(tr, PAZ, None, 10)
            elif seedresp and len(tr.data) > 10:
                tr=_sim_WA(tr, None, seedresp, 10)
            elif len(tr.data) > 10:
                warnings.warn('No PAZ for '+tr.stats.station+' '+\
                                 tr.stats.channel+' at time: '+\
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
            amplitude, period, delay= _max_p2t(tr.data, tr.stats.delta)
            if amplitude==0.0:
                break
            print 'Amplitude picked: '+str(amplitude)
            # Note, amplitude should be in meters at the moment!
            # Remove the pre-filter response
            if pre_filt:
                # Generate poles and zeros for the filter we used earlier - this
                # is how the filter is designed in the convenience methods of
                # filtering in obspy.
                z, p, k=iirfilter(corners, [lowcut/(0.5*tr.stats.sampling_rate),\
                                            highcut/(0.5*tr.stats.sampling_rate)],\
                                  btype='band', ftype='butter', output='zpk')
                filt_paz={'poles': list(p),
                          'zeros': list(z),
                          'gain': k,
                          'sensitivity':  1.0}
                amplitude /= (paz2AmpValueOfFreqResp(filt_paz, 1/period) * \
                              filt_paz['sensitivity'])
            # Convert amplitude to mm
            if PAZ: # Divide by Gain to get to nm (returns pm? 10^-12)
                # amplitude *=PAZ['gain']
                amplitude /= 1000
            if seedresp: # Seedresp method returns mm
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
                                         velocity=float('NaN'), AIN=999, SNR='',
                                         azimuthres=999, timeres=float('NaN'),
                                         finalweight=999, distance=hypo_dist,
                                         CAZ=CAZ))
    # Copy the header from the sfile to a new local S-file
    fin=open(sfile,'r')
    fout=open('mag_calc.out','w')
    for line in fin:
        if not line[79]=='7':
            fout.write(line)
        else:
            fout.write(line)
            break
    fin.close()
    fout.close()
    # Write picks out to new s-file
    for pick in picks_out:
        print pick
    Sfile_util.populateSfile('mag_calc.out',picks_out)
    return picks

if __name__ == '__main__':
    """
    Code to loop through a database and make a shit-tonne of magnitude picks-boi!

    Coded to run through SEISAN databases.
    """
    import sys
    if len(sys.argv) < 6 or len(sys.argv) > 7:
        msg='Insufficient arguments, needs the database to calculate over, and'+\
            ' the ouptut database, paths to REA dir (not year/mm dirs) for both'+\
            ' please, and the path to the CAL directory, and the start and'+\
            ' stop dates as yyyymmddhhmmss'
        raise IOError(msg)
    indir=str(sys.argv[1])
    outdir=str(sys.argv[2])
    calpath=str(sys.argv[3])
    startdate=str(sys.argv[4])
    enddate=str(sys.argv[5])
    if len(sys.argv) == 7:
    	wavepath=sys.argv[6]
    elif len(sys.argv) == 6:
    	wavepath=None
    import glob, shutil
    import datetime as dt
    try:
        startdate=dt.datetime.strptime(startdate.ljust(14,'0'), '%Y%m%d%H%M%S')
    except:
        raise IOError('start date is not yyyymmddhhmmss form')
    try:
        stopdate=dt.datetime.strptime(enddate.ljust(14,'0'), '%Y%m%d%H%M%S')
    except:
        raise IOError('end date is not yyyymmddhhmmss form')
    kdays=((stopdate+dt.timedelta(1))-startdate).days
    for i in xrange(kdays):
        day=startdate+dt.timedelta(i)
        print 'Working on '+str(day)
        sfiles=glob.glob(indir+'/'+str(day.year)+'/'+str(day.month).zfill(2)+\
                         '/'+str(day.day).zfill(2)+'-*L.S'+str(day.year)+\
                         str(day.month).zfill(2))
        datetimes=[dt.datetime.strptime(sfiles[i].split('/')[-1], '%d-%H%M-%SL.S%Y%m')\
                   for i in xrange(len(sfiles))]
        sfiles=[sfiles[i] for i in xrange(len(sfiles)) if datetimes[i] > startdate and
                    datetimes[i] < stopdate]
        if not wavepath:
        	wavedir="/".join(indir.split('/')[:-2])+'/WAV/'+\
            	    indir.split('/')[-1]+'/'+str(day.year)+'/'+\
                	str(day.month).zfill(2)
        else:
        	wavedir=wavepath+'/'+str(day.year)+'/'+\
                	str(day.month).zfill(2)
        sfiles.sort()
        for sfile in sfiles:
            # Make the picks!
            print '				Working on Sfile: '+sfile
            picks=Amp_pick_sfile(sfile, wavedir, calpath)
            # Copy the mag_calc.out file to the correct place
            shutil.copyfile('mag_calc.out', outdir+'/'+str(day.year)+'/'+\
                            str(day.month).zfill(2)+'/'+sfile.split('/')[-1])
