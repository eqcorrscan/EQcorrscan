#!/usr/bin/env python
"""
Script to validate the synthetic detector method against test earthquakes
"""
import warnings

def synth_from_sfile(sfile, samp_rate, length=10.0, PS_ratio=1.68):
    """
    Function to generate a synthetic template for a given s-file

    :type path: str
    :param path: Path to the sfile
    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz for template
    :type length: float
    :param length: Length of templates in seconds, defaults to 10
    :type PS_ratio: float
    :param PS_ratio: S-P ratio to use if only one pick found for station,\
            defaults to 1.68

    :returns: :class: obspy.Stream
    """
    from utils import Sfile_util
    from utils import synth_seis
    from obspy import Stream, Trace, UTCDateTime
    # Get the picks and the origin time
    picks=Sfile_util.readpicks(sfile)
    ori_time=Sfile_util.readheader(sfile).time.datetime
    # We only want P and S phases
    picks=[p for p in picks if p.phase in ['P','S']]
    # We want a list of the stations that we have picks for
    stations=list(set([p.station for p in picks]))
    # Loop through the stations
    synths=Stream()
    for station in stations:
        # Find the relevant picks
        sta_picks=[p for p in picks if p.station==station]
        if len(sta_picks) == 1:
            msg='Only '+sta_picks[0].phase+' phase picked for station '+\
                station+' will use an S-P ratio of '+str(PS_ratio)
            warnings.warn(msg)
            # Calculate the pick travel time
            tt = (sta_picks[0].time.datetime-ori_time).total_seconds()
            if sta_picks[0].phase == 'P':
                SP_time=(tt*PS_ratio)-tt
                P_pick=sta_picks[0].time.datetime
                S_pick=(UTCDateTime(P_pick)+SP_time).datetime
            else:
                SP_time=tt-(tt/PS_ratio)
                P_pick=(sta_picks[0].time-SP_time).datetime
                S_pick=sta_picks[0].time.datetime
        else:
            if len([p for p in sta_picks if p.phase=='P']) > 1:
                warnings.warn('Multiple P picks found for station '+station+\
                                ', will use earliest')
                P_pick=min([p.time for p in sta_picks if p.phase=='P'])
                channel=sta_picks[0].channel
            else:
                P_pick=[p.time for p in sta_picks if p.phase=='P'][0]
                channel=sta_picks[0].channel
            if len([p for p in sta_picks if p.phase=='S']) > 1:
                warnings.warn('Multiple S picks found for station '+station+\
                                ', will use earliest')
                S_pick=min([p.time for p in sta_picks if p.phase=='S'])
            else:
                S_pick=[p.time for p in sta_picks if p.phase=='S'][0]
            if P_pick > S_pick:
                raise ValueError('P pick is after S pick')
            SP_time=(S_pick.datetime-P_pick.datetime).total_seconds()
        # Loop through the picks available
        for p in sta_picks:
            tr=Trace(synth_seis.seis_sim(int(SP_time*samp_rate),\
                        flength=length*samp_rate, phaseout=p.phase))
            tr.stats.sampling_rate=samp_rate
            tr.stats.station=station
            tr.stats.channel=p.channel
            # Sythetics start 10 samples before P
            if p.phase in ['all', 'P']:
                tr.stats.starttime=UTCDateTime(P_pick)-(10.0/samp_rate)
            else:
                tr.stats.starttime=UTCDateTime(S_pick)-(10.0/samp_rate)
            synths+=tr
    return synths

def match_synth(sfile, cont_base, freqmin=2.0, freqmax=10.0, samp_rate=100.0,\
                threshold=8.0, threshold_type='MAD', trig_int=6.0, plotvar=True,\
                save_template=True):
    """
    Function to generate a basic synthetic from a real event, given by an s-file
    and cross-correlate this with the day of continuous data including the event

    :type sfile: str
    :param sfile: Path to the s-file for the event
    :type cont_base: str
    :param cont_base: Path to the continuous data, should be in Yyyyy/Rjjj.01\
                directories
    :type freqmin: float
    :param freqmin: Low-cut for bandpass in Hz, defualts to 2.0
    :type freqmax: float
    :param freqmax: High-cut for bandpass in Hz, defaults to 10.0
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz, defaults to 100.0
    :type threshold: float
    :param threshold: Threshold for detection in cccsum, defaults to 8.0
    :type threshold_type: str
    :param threshold_type: Type to threshold, either MAD or ABS, defaults to MAD
    :type trig_int: float
    :param trig_int: Trigger interval in seconds, defaults to 6.0
    :type plotvar: bool
    :param plotvar: To plot or not, defaults to true

    :returns: detections
    """
    # import matplotlib.pyplot as plt
    from core import match_filter, template_gen
    from utils import Sfile_util, pre_processing
    import glob
    from obspy import read, Stream, UTCDateTime
    from obspy.signal.cross_correlation import xcorr
    from joblib import Parallel, delayed
    from multiprocessing import cpu_count
    import numpy as np
    # Generate the synthetic
    synth_template=synth_from_sfile(sfile, samp_rate, length=1.0,\
                                    PS_ratio=1.68)
    synth_template.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
    for tr in synth_template:
        tr.data=(tr.data*1000).astype(np.int32)
    if save_template:
        synth_template.write('Synthetic_'+sfile.split('/')[-1],
                            format='MSEED', encoding='STEIM2')
    # Find the date from the sfile
    event_date=Sfile_util.readheader(sfile).time.datetime
    day=UTCDateTime(event_date.date())
    # Work out which stations we have template info for
    stachans=[(tr.stats.station, tr.stats.channel) for tr in synth_template]
    # Read in the day of data
    for stachan in stachans:
        wavfile=glob.glob(cont_base+event_date.strftime('/Y%Y/R%j.01/')+\
                            stachan[0]+'.*.'+stachan[1][0]+'?'+stachan[1][-1]+\
                            '.'+event_date.strftime('%Y.%j'))
        if len(wavfile) != 0:
            for wavf in wavfile:
                if not 'st' in locals():
                    st=read(wavf)
                else:
                    st+=read(wavf)
    st=st.merge(fill_value='interpolate')
    cores=cpu_count()
    if len(st) < cores:
        jobs=len(st)
    else:
        jobs=cores
    st=Parallel(n_jobs=jobs)(delayed(pre_processing.dayproc)(tr, freqmin,\
                                                             freqmax, 3,\
                                                             samp_rate, 0,\
                                                             day)
                            for tr in st)
    st=Stream(st)
    # Make the real template
    picks=Sfile_util.readpicks(sfile)
    real_template=template_gen._template_gen(picks, st, 1.0, 'all',\
                                            prepick=10/samp_rate)
    for tr in real_template:
        tr.data=tr.data.astype(np.int32)
    if save_template:
        real_template.write('Real_'+sfile.split('/')[-1], format='MSEED',\
                            encoding='STEIM2')
    # Shift the synthetic to better align with the real one
    for tr in real_template:
        synth_tr=synth_template.select(station=tr.stats.station,\
                                        channel=tr.stats.channel)[0]
        shift, corr = xcorr(tr.data, synth_tr.data, 20)
        print tr.stats.station+'.'+tr.stats.channel+' shift='+str(shift)+' corr='+str(corr)
        if corr < 0:
            synth_tr.data*=-1
        # Apply a pad
        pad=np.zeros(abs(shift))
        if shift < 0:
            synth_tr.data=np.append(synth_tr.data, pad)[abs(shift):]
        elif shift > 0:
            synth_tr.data=np.append(pad, synth_tr.data)[0:-shift]
        # plt.plot(tr.data, 'k')
        # plt.plot(synth_tr.data*(max(tr.data)/max(synth_tr.data)), 'r')
        # print max(tr.data)
        # print max(synth_tr.data)
        # plt.title(tr.stats.station+'.'+tr.stats.channel+' Correlates at '+str(corr))
        # plt.show()
    # Now we have processed data and a template, we can try and detect!
    detections=match_filter.match_filter(['Synthetic_'+sfile.split('/')[-1],
                                        'Real_'+sfile.split('/')[-1]],\
                                        [synth_template, real_template],\
                                        st, threshold, \
                                        threshold_type, trig_int,\
                                        plotvar, 'synth_temp')
    f=open('Synthetic_test.csv', 'w')
    f.write('template, detect-time, cccsum, threshold, number of channels\n')
    for detection in detections:
        # output detections to file
        f.write(detection.template_name+', '+str(detection.detect_time)+\
                ', '+str(detection.detect_val)+', '+str(detection.threshold)+\
                ', '+str(detection.no_chans)+'\n')
        print 'template: '+detection.template_name+' detection at: '\
            +str(detection.detect_time)+' with a cccsum of: '+\
            str(detection.detect_val)
    if detections:
        f.write('\n')
    f.close()


if __name__ =='__main__':
    import sys
    if not len(sys.argv)==10:
        raise IOError('Insufficient arguments '+str(len(sys.argv)))
    sfile=str(sys.argv[1])
    cont_base=str(sys.argv[2])
    freqmin=float(sys.argv[3])
    freqmax=float(sys.argv[4])
    samp_rate=float(sys.argv[5])
    threshold=float(sys.argv[6])
    threshold_type=str(sys.argv[7])
    trig_int=float(sys.argv[8])
    plotvar=sys.argv[9]
    match_synth(sfile, cont_base, freqmin, freqmax, samp_rate, threshold,\
                threshold_type, trig_int, plotvar)
