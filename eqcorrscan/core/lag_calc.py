#!/usr/bin/python
"""
Functions to generate lag-times for events detected by correlation.
Part of the EQcorrscan module to integrate seisan nordic files into a full
cross-channel correlation for detection routine.
EQcorrscan is a python module designed to run match filter routines for
seismology, within it are routines for integration to seisan and obspy.
With obspy integration (which is necessary) all main waveform formats can be
read in and output.
This main section contains a script, LFE_search.py which demonstrates the usage
of the built in functions from template generation from picked waveforms
through detection by match filter of continuous data to the generation of lag
times to be used for relative locations.
The match-filter routine described here was used a previous Matlab code for the
Chamberlain et al. 2014 G-cubed publication.  The basis for the lag-time
generation section is outlined in Hardebeck & Shelly 2011, GRL.
Code written by Calum Chamberlain and Chet Hopp, VUW, 2016.
.. rubric:: Note
Pre-requisites:
    - gcc             - for the installation of the openCV correlation routine
    - python-cv2      - Python bindings for the openCV routines
    - python-joblib   - used for parallel processing
    - python-obspy    - used for lots of common seismological processing
                        - requires:
                            - numpy
                            - scipy
                            - matplotlib
    - NonLinLoc       - used outside of all codes for travel-time generation
Copyright 2015, 2016, the authors.
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
from match_filter import normxcorr2


def _channel_loop(detection, template, min_cc, i=0):
    """
    Utility function to take a stream of data for the detected event and write
    maximum correlation to absolute time as picks in an obspy.core.event.Event
    object
    :type detection: obspy.Stream
    :param detection: Stream of data for the slave event detected using \
    template
    :type template: obspy.Stream
    :param template: Stream of data as the template for the detection.
    :type i: int, optional
    :param i: Used to track which process has occured when running in parallel
    :returns: Event object containing net, sta, chan information
    """
    from obspy.core.event import Event, Pick, WaveformStreamID
    from obspy.core.event import ResourceIdentifier
    event = Event()
    s_stachans = {}
    used_s_sta = []
    for tr in template:
        temp_net = tr.stats.network
        temp_sta = tr.stats.station
        temp_chan = tr.stats.channel
        image = detection.select(station=temp_sta,
                                 channel=temp_chan)
        if image:  # Ideally this if statement would be removed.
            ccc = normxcorr2(tr.data, image[0].data)
            # shiftlen = len(ccc)*image[0].stats.sampling_rate
            # Convert the maximum cross-correlation time to an actual time
            if np.amax(ccc) > min_cc:
                picktime = image[0].stats.starttime + (np.argmax(ccc) *
                                                       image[0].stats.delta)
            else:
                picktime = image[0].stats.starttime
            # Perhaps weight each pick by the cc val or cc val^2?
            # weight = np.amax(ccc) ** 2
            if temp_chan[-1:] == 'Z':
                phase = 'P'
            # Only take the S-pick with the best correlation
            elif temp_chan[-1:] in ['E', 'N']:
                phase = 'S'
                if temp_sta not in s_stachans and np.amax(ccc) > min_cc:
                    s_stachans[temp_sta] = ((temp_chan, np.amax(ccc),
                                             picktime))
                elif temp_sta in s_stachans and np.amax(ccc) > min_cc:
                    if np.amax(ccc) > s_stachans[temp_sta][1]:
                        picktime = picktime
                    else:
                        picktime = s_stachans[temp_sta][2]
                        temp_chan = s_stachans[temp_sta][0]
                elif np.amax(ccc) < min_cc and temp_sta not in used_s_sta:
                    used_s_sta.append(temp_sta)
                else:
                    continue
            _waveform_id = WaveformStreamID(network_code=temp_net,
                                            station_code=temp_sta,
                                            channel_code=temp_chan)
            event.picks.append(Pick(waveform_id=_waveform_id,
                                    time=picktime,
                                    method_id=ResourceIdentifier('EQcorrscan'),
                                    phase_hint=phase))
    return (i, event)


def _day_loop(detection_streams, template, min_cc):
    """
    Function to loop through multiple detections for one template - ostensibly
    designed to run for the same day of data for I/O simplicity, but as you
    are passing stream objects it could run for all the detections ever, as
    long as you have the RAM!
    :type detection_streams: List of obspy.Stream
    :param detection_streams: List of all the detections for this template that
                              you want to compute the optimum pick for.
    :type template: obspy.Stream
    :param template: The original template used to detect the detections passed
    :returns: Catalog object containing Event objects for each detection
              created by this template.
    """
    from multiprocessing import Pool, cpu_count
    # Used to run detections in parallel
    from obspy.core.event import Catalog
    num_cores = cpu_count()
    if num_cores > len(detection_streams):
        num_cores = len(detection_streams)
    pool = Pool(processes=num_cores, maxtasksperchild=None)
    # Parallelize generation of events for each detection:
    # results is a list of (i, event class)
    results = [pool.apply_async(_channel_loop, args=(detection_streams[i],
                                                     template, min_cc, i))
               for i in xrange(len(detection_streams))]
    pool.close()
    events_list = [p.get() for p in results]
    events_list.sort(key=lambda tup: tup[0])  # Sort based on i.
    temp_catalog = Catalog()
    temp_catalog.events = [event_tup[1] for event_tup in events_list]
    # Generates moveout times from event.picks.
    # mintime = template.sort(['starttime'])[0].stats.starttime
    # lags = []
    # for event in events:
    #     delay = template.select(station=lag,
    #                             channel=lag[3])[0].stats.starttime - mintime
    #     lag[0] += delay
    return temp_catalog


def lag_calc(detections, detect_data, templates, shift_len=0.2, min_cc=0.4):
    r"""
    Overseer function to take a list of detection objects, cut the data for
    them to lengths of the same length of the template + shift_len on
    either side. This will then write out SEISAN s-file or QuakeML for the
    detections with pick times based on the lag-times found at the maximum
    correlation, providing that correlation is above the min_cc.
    :type detections: List of DETECTION
    :param detections: List of DETECTION objects
    :type detect_data: obspy.Stream
    :param detect_data: All the data needed to cut from - can be a gappy Stream
    :type templates: List of tuple of String, obspy.Stream
    :param templates: List of the templates used as tuples of template name, \
        template
    :type shift_len: float
    :param shift_len: Shift length allowed for the pick in seconds, will be
                    plus/minus this amount - default=0.2
    :type min_cc: float
    :param min_cc: Minimum cross-correlation value to be considered a pick,
                    default=0.4
    :returns: obspy.core.event.Catalog of events with picks.  No origin\
        imformation is included, these events can then be written out via\
        obspy.core.event functions, or to seisan Sfiles using Sfile_util\
        and located.
    """
    from eqcorrscan.utils import Sfile_util
    from obspy import Stream
    from obspy.core.event import Catalog

    # Establish plugin directory relative to this module

    # First work out the delays for each template
    delays = []  # List of tuples of (tempname, (sta, chan, delay))
    for template in templates:
        temp_delays = []
        for tr in template[1]:
            temp_delays.append((tr.stats.station, tr.stats.channel,
                                tr.stats.starttime - template[1].
                                sort(['starttime'])[0].stats.starttime))
        delays.append((template[0], temp_delays))
    # List of tuples of (template name, Stream()) for each detection
    detect_streams = []
    for detection in detections:
        # Stream to be saved for new detection
        detect_stream = []
        for tr in detect_data:
            tr_copy = tr.copy()
            # Right now, copying each trace hundreds of times...
            template = [t for t in templates if t[0] == detection.
                        template_name][0]
            template = template[1].select(station=tr.stats.station,
                                          channel=tr.stats.channel)
            if template:
                # Save template trace length in seconds
                template_len = len(template[0]) / \
                    template[0].stats.sampling_rate
            else:
                continue
                # If there is no template-data match then skip the rest
                # of the trace loop.
            # Grab the delays for the desired template: [(sta, chan, delay)]
            delay = [delay for delay in delays if delay[0] == detection.
                     template_name][0][1]
            # Now grab the delay for the desired trace for this template
            delay = [d for d in delay if d[0] == tr.stats.station and
                     d[1] == tr.stats.channel][0][2]
            detect_stream.append(tr_copy.trim(starttime=detection.detect_time -
                                              shift_len + delay,
                                              endtime=detection.detect_time +
                                              delay + shift_len +
                                              template_len))
            del tr_copy
        # Create tuple of (template name, data stream)
        detect_streams.append((detection.template_name, Stream(detect_stream)))
    # Segregate detections by template, then feed to day_loop
    initial_cat = Catalog()
    for template in templates:
        template_detections = [detect[1] for detect in detect_streams
                               if detect[0] == template[0]]
        if len(template_detections) > 0:  # Way to remove this?
            initial_cat += _day_loop(template_detections, template[1], min_cc)
    return initial_cat
