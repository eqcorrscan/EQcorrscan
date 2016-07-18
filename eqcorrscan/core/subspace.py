#!/usr/bin/python
r"""This module contains functions relevant to executing subspace detection \
for earthquake catalogs. The overarching function calls _channel_loop and \
_template_loop inner functions from match_filter in the same way as \
the normal matched filtering workflow.

We recommend that you read Harris' detailed report on subspace detection \
theory which can be found here: https://e-reports-ext.llnl.gov/pdf/335299.pdf

:copyright:
    Calum Chamberlain, Chet Hopp.

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


def det_statistic(detector, data):
    """
    Base function to calculate the subspace detection statistic.

    Calculates for a given subspace detector and data stream. \
    The statistic is calculated by \
    projecting the data onto the N dimensional subspace defined by the given \
    detector following the equation: :math:'\\gamma = y^TUU^Ty' where y is \
    the data stream, U is the subspace detector and :math:'\\gamma' is the \
    detection statistic from 0 to 1.
    """
    day_stats = []
    for i in range(len(data) - len(detector[0]) + 1):
        y = data[i:i + len(detector[0])]
        y = y / np.sqrt((y ** 2).sum(-1))
        day_stats.append(y.dot(detector.T).dot(detector).dot(y))
    day_stats = np.asarray(day_stats)
    if np.all(np.isnan(day_stats)):
        day_stats = np.zeros(len(day_stats))
    return day_stats


def subspace_detect(detector_names, detector_list, st, threshold,
                    threshold_type, trig_int, plotvar, plotdir='.', cores=1,
                    tempdir=False, debug=0, plot_format='jpg',
                    output_cat=False, extract_detections=False):
    """
    Overseer function to handle subspace detection.

    Modelled after match_filter.match_filter().
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.ioff()
    import copy
    from eqcorrscan.core import match_filter
    from eqcorrscan.core.match_filter import DETECTION
    from eqcorrscan.utils import plotting, findpeaks
    from obspy import Trace, Catalog, UTCDateTime, Stream
    from obspy.core.event import Event, Pick, CreationInfo, ResourceIdentifier
    from obspy.core.event import Comment, WaveformStreamID
    import time

    # Copy the stream here because we will muck about with it
    stream = st.copy()
    detectors = copy.deepcopy(detector_list)
    # Collect detector and data channel information
    detector_stachan = list((tr.stats.station, tr.stats.channel)
                            for detector in detectors
                            for det_st in detector
                            for tr in det_st)
    data_stachan = list((tr.stats.station, tr.stats.channel) for tr in stream)
    detector_stachan = list(set(detector_stachan))
    data_stachan = list(set(data_stachan))
    if debug >= 3:
        print('I have detector info for these stations:')
        print(detector_stachan)
        print('I have daylong data for these stations:')
        print(data_stachan)
    # Perform a check that the daylong vectors are daylong
    for tr in stream:
        if not tr.stats.sampling_rate * 86400 == tr.stats.npts:
            msg = ' '.join(['Data are not daylong for', tr.stats.station,
                            tr.stats.channel])
            raise ValueError(msg)
    # Check that detector_list is list of lists of obspy.Stream
    for i, detector in enumerate(detectors):
        for j, det_st in enumerate(detector):
            if isinstance(det_st, Stream):
                continue
            else:
                msg = ('Detector %s: stream %02d is not actually a stream! \
                       Check that detector_list is being given a list of \
                       lists of obspy.Stream objects.' %
                       (detector_names[i], j))
    # Perform check that all detector lengths are internally consistent
    for i, det in enumerate(detectors):
        for j, det_st in enumerate(det):
            if len(set([tr.stats.npts for tr in det_st])) > 1:
                msg = 'Detector %s, stream %02d contains traces of differing \
                      length!! THIS WILL CAUSE ISSUES' % (detector_names[i], j)
                raise ValueError(msg)
    # Call the _template_loop function to do all the correlation work
    outtic = time.clock()
    # Edit here from previous, stable, but slow match_filter
    # Would be worth testing without an if statement, but with every station in
    # the possible template stations having data, but for those without real
    # data make the data NaN to return NaN ccc_sum
    # Note: this works
    if debug >= 2:
        print('Ensuring all detector channels have matches in daylong data')
    for stachan in detector_stachan:
        if not stream.select(station=stachan[0], channel=stachan[1]):
            # Remove detector traces rather than adding NaN data
            if debug >= 3:
                print('Removing %s from detectors' % str(stachan))
            for detector in detectors:
                for det_st in detector:
                    if det_st.select(station=stachan[0], channel=stachan[1]):
                        for tr in det_st.select(station=stachan[0],
                                                channel=stachan[1]):
                            det_st.remove(tr)
    # Remove un-needed channels
    for tr in stream:
        if not (tr.stats.station, tr.stats.channel) in detector_stachan:
            if debug >= 3:
                print('Removed sta %s, chan %s from data' % (tr.stats.station,
                                                             tr.stats.channel))
            stream.remove(tr)
    # Also pad out detectors to have all channels
    for detector, detector_name in zip(detectors, detector_names):
        for i, det_st in enumerate(detector):
            if debug >= 3:
                print('Working out padding for detector %s, svec %02d'
                      % (detector_name, i))
            if len(det_st) == 0:
                msg = ('No channels matching in continuous data for ' +
                       'detector %s: stream %02d' % (detector_name, i))
                warnings.warn(msg)
                detectors.remove(detector)
                detector_names.remove(detector_name)
                continue
            for stachan in detector_stachan:
                if not det_st.select(station=stachan[0], channel=stachan[1]):
                    nulltrace = Trace()
                    nulltrace.stats.station = stachan[0]
                    nulltrace.stats.channel = stachan[1]
                    nulltrace.stats.sampling_rate = det_st[0].stats.sampling_rate
                    nulltrace.stats.starttime = det_st[0].stats.starttime
                    nulltrace.data = np.array([np.NaN] * len(det_st[0].data),
                                              dtype=np.float32)
                    det_st += nulltrace
    if debug >= 2:
        print('Starting the correlation run for this day')
    [cccsums, no_chans, chans] = match_filter._channel_loop(detectors, stream,
                                                            cores,
                                                            do_subspace=True,
                                                            debug=debug)
    if len(cccsums[0]) == 0:
        raise ValueError('Correlation has not run, zero length cccsum')
    outtoc = time.clock()
    print(' '.join(['Looping over detectors and streams took:',
                    str(outtoc - outtic), 's']))
    if debug >= 2:
        print(' '.join(['The shape of the returned cccsums is:',
                        str(np.shape(cccsums))]))
        print(' '.join(['This is from', str(len(detectors)), 'detectors']))
        print(' '.join(['Correlated with', str(len(stream)),
                        'channels of data']))
    detections = []
    if output_cat:
        det_cat = Catalog()
    # XXX TODO: Need to investigate thresholding for correct detection stats
    for i, cccsum in enumerate(cccsums):
        detector = detectors[i]
        # Unless I can prove otherwise, we will average cccsums to obtain our detection statistic from 0 to 1
        cccsum /= no_chans[i]
        if threshold_type == 'MAD':
            rawthresh = threshold * np.median(np.abs(cccsum))
        elif threshold_type == 'absolute':
            rawthresh = threshold
        elif threshold == 'av_chan_corr':
            rawthresh = threshold * (cccsum / len(detector[0]))
        else:
            print('You have not selected the correct threshold type, I will' +
                  'use MAD as I like it')
            rawthresh = threshold * np.mean(np.abs(cccsum))
        # Findpeaks returns a list of tuples in the form [(cccsum, sample)]
        print(' '.join(['Threshold is set at:', str(rawthresh)]))
        print(' '.join(['Max of data is:', str(max(cccsum))]))
        print(' '.join(['Mean of data is:', str(np.mean(cccsum))]))
        # Set up a trace object for the cccsum as this is easier to plot and
        # maintains timing
        if plotvar:
            stream_plot = copy.deepcopy(stream[0])
            # Downsample for plotting
            stream_plot.decimate(int(stream[0].stats.sampling_rate / 10))
            cccsum_plot = Trace(cccsum)
            cccsum_plot.stats.sampling_rate = stream[0].stats.sampling_rate
            # Resample here to maintain shape better
            cccsum_hist = cccsum_plot.copy()
            cccsum_hist = cccsum_hist.decimate(int(stream[0].stats.
                                                   sampling_rate / 10)).data
            cccsum_plot = plotting.chunk_data(cccsum_plot, 10,
                                              'Maxabs').data
            # Enforce same length
            stream_plot.data = stream_plot.data[0:len(cccsum_plot)]
            cccsum_plot = cccsum_plot[0:len(stream_plot.data)]
            cccsum_hist = cccsum_hist[0:len(stream_plot.data)]
            plotting.triple_plot(cccsum_plot, cccsum_hist,
                                 stream_plot, rawthresh, True,
                                 plotdir + '/cccsum_plot_' +
                                 detector_names[i] + '_' +
                                 stream[0].stats.starttime.
                                 datetime.strftime('%Y-%m-%d') +
                                 '.' + plot_format)
            if debug >= 4:
                print(' '.join(['Saved the cccsum to:', detector_names[i],
                                stream[0].stats.starttime.datetime.
                                strftime('%Y%j')]))
                np.save(detector_names[i] +
                        stream[0].stats.starttime.datetime.strftime('%Y%j'),
                        cccsum)
        tic = time.clock()
        if debug >= 4:
            np.save('cccsum_' + str(i) + '.npy', cccsum)
        if debug >= 3 and max(cccsum) > rawthresh:
            peaks = findpeaks.find_peaks2_short(cccsum, rawthresh,
                                                trig_int * stream[0].stats.
                                                sampling_rate, debug,
                                                stream[0].stats.starttime,
                                                stream[0].stats.sampling_rate)
        elif max(cccsum) > rawthresh:
            peaks = findpeaks.find_peaks2_short(cccsum, rawthresh,
                                                trig_int * stream[0].stats.
                                                sampling_rate, debug)
        else:
            print('No peaks found above threshold')
            peaks = False
        toc = time.clock()
        if debug >= 1:
            print(' '.join(['Finding peaks took:', str(toc - tic), 's']))
        if peaks:
            for peak in peaks:
                detecttime = stream[0].stats.starttime +\
                    peak[1] / stream[0].stats.sampling_rate
                rid = ResourceIdentifier(id=detector_names[i] + '_' +
                                         str(detecttime),
                                         prefix='smi:local')
                ev = Event(resource_id=rid)
                cr_i = CreationInfo(author='EQcorrscan',
                                    creation_time=UTCDateTime())
                ev.creation_info = cr_i
                # All detection info in Comments for lack of a better idea
                thresh_str = 'threshold=' + str(rawthresh)
                ccc_str = 'detect_val=' + str(peak[0])
                used_chans = 'channels used: ' +\
                             ' '.join([str(pair) for pair in chans[i]])
                ev.comments.append(Comment(text=thresh_str))
                ev.comments.append(Comment(text=ccc_str))
                ev.comments.append(Comment(text=used_chans))
                for det_st in detector:
                    for tr in det_st:
                        if (tr.stats.station, tr.stats.channel) not in chans[i]:
                            continue
                        else:
                            pick_tm = detecttime + (tr.stats.starttime -
                                                    detecttime)
                            wv_id = WaveformStreamID(network_code=tr.stats.network,
                                                     station_code=tr.stats.station,
                                                     channel_code=tr.stats.channel)
                            ev.picks.append(Pick(time=pick_tm, waveform_id=wv_id))
                detections.append(DETECTION(detector_names[i],
                                            detecttime,
                                            no_chans[i], peak[0], rawthresh,
                                            'corr', chans[i], event=ev))
                if output_cat:
                    det_cat.append(ev)
        if extract_detections:
            detection_streams = extract_from_stream(stream, detections)
    del stream, detectors
    if output_cat and not extract_detections:
        return detections, det_cat
    elif not extract_detections:
        return detections
    elif extract_detections and not output_cat:
        return detections, detection_streams
    else:
        return detections, det_cat, detection_streams


if __name__ == '__main__':
    import doctest
    doctest.testmod()