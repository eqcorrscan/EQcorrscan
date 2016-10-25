"""
Functions to generate pick-corrections for events detected by correlation.

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
import scipy
import warnings
import logging
import sys

from multiprocessing import Pool, cpu_count
from collections import Counter

from obspy import Stream
from obspy.core.event import Catalog
from obspy.core.event import Event, Pick, WaveformStreamID
from obspy.core.event import ResourceIdentifier, Comment

from eqcorrscan.utils.plotting import plot_repicked, detection_multiplot
from eqcorrscan.core.match_filter import normxcorr2

# Set up logging
log = logging.getLogger(__name__)
ch = logging.StreamHandler(stream=sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s -' +
                              ' %(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)


class LagCalcError(Exception):
    """
    Default error for issues within lag-calc.
    """
    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return self.value

    def __str__(self):
        return 'LagCalcError: ' + self.value


def _xcorr_interp(ccc, dt):
    """
    Intrpolate around the maximum correlation value for sub-sample precision.

    :param ccc: Cross-correlation array
    :type ccc: numpy.ndarray
    :param dt: sample interval
    :type dt: float

    :return: Position of interpolated maximum in seconds from start of ccc
    :rtype: float
    """
    if ccc.shape[0] == 1:
        cc = ccc[0]
    else:
        cc = ccc
    # Code borrowed from obspy.signal.cross_correlation.xcorr_pick_correction
    cc_curvature = np.concatenate((np.zeros(1), np.diff(cc, 2), np.zeros(1)))
    cc_t = np.arange(0, len(cc) * dt, dt)
    peak_index = cc.argmax()
    first_sample = peak_index
    # XXX this could be improved..
    while first_sample > 0 and cc_curvature[first_sample - 1] <= 0:
        first_sample -= 1
    last_sample = peak_index
    while last_sample < len(cc) - 1 and cc_curvature[last_sample + 1] <= 0:
        last_sample += 1
    num_samples = last_sample - first_sample + 1
    if num_samples < 3:
        msg = "Less than 3 samples selected for fit to cross " + \
              "correlation: %s" % num_samples
        raise IndexError(msg)
    if num_samples < 5:
        msg = "Less than 5 samples selected for fit to cross " + \
              "correlation: %s" % num_samples
        warnings.warn(msg)
    coeffs, residual = scipy.polyfit(
        cc_t[first_sample:last_sample + 1],
        cc[first_sample:last_sample + 1], deg=2, full=True)[:2]
    # check results of fit
    if coeffs[0] >= 0:
        msg = "Fitted parabola opens upwards!"
        warnings.warn(msg)
    if residual > 0.1:
        msg = "Residual in quadratic fit to cross correlation maximum " + \
              "larger than 0.1: %s" % residual
        warnings.warn(msg)
    # X coordinate of vertex of parabola gives time shift to correct
    # differential pick time. Y coordinate gives maximum correlation
    # coefficient.
    shift = -coeffs[1] / 2.0 / coeffs[0]
    coeff = (4 * coeffs[0] * coeffs[2] - coeffs[1] ** 2) / (4 * coeffs[0])
    return shift, coeff


def _channel_loop(detection, template, min_cc, detection_id, interpolate, i,
                  pre_lag_ccsum=None, detect_chans=0):
    """
    Inner loop for correlating and assigning picks.

    Utility function to take a stream of data for the detected event and write
    maximum correlation to absolute time as picks in an obspy.core.event.Event
    object.
    Only outputs picks for picks above min_cc.

    :type detection: obspy.core.stream.Stream
    :param detection:
        Stream of data for the slave event detected using template.
    :type template: obspy.core.stream.Stream
    :param template: Stream of data as the template for the detection.
    :type min_cc: float
    :param min_cc: Minimum cross-correlation value to allow a pick to be made.
    :type detection_id: str
    :param detection_id: Detection ID to associate the event with.
    :type interpolate: bool
    :param interpolate:
        Interpolate the correlation function to achieve sub-sample precision.
    :type i: int
    :param i:
        Used to track which process has occurred when running in parallel.
    :type pre_lag_ccsum: float
    :param pre_lag_ccsum:
        Cross-correlation sum before lag-calc, will check that the
        cross-correlation sum is increased by lag-calc (using all channels,
        ignoring min_cc)
    :type detect_chans: int
    :param detect_chans:
        Number of channels originally used in detections, must match the number
        used here to allow for cccsum checking.

    :returns:
        Event object containing network, station, channel and pick information.
    :rtype: :class:`obspy.core.event.Event`
    """
    event = Event()
    s_stachans = {}
    used_s_sta = []
    cccsum = 0
    checksum = 0
    used_chans = 0
    for tr in template:
        temp_net = tr.stats.network
        temp_sta = tr.stats.station
        temp_chan = tr.stats.channel
        image = detection.select(station=temp_sta,
                                 channel=temp_chan)
        if image:
            if interpolate:
                try:
                    ccc = normxcorr2(tr.data, image[0].data)
                    shift, cc_max = _xcorr_interp(ccc=ccc,
                                                  dt=image[0].stats.delta)
                except IndexError:
                    log.error('Could not interpolate ccc, not smooth')
                    ccc = normxcorr2(tr.data, image[0].data)
                    cc_max = np.amax(ccc)
                    shift = np.argmax(ccc) * image[0].stats.delta
                # Convert the maximum cross-correlation time to an actual time
                picktime = image[0].stats.starttime + shift
            else:
                # Convert the maximum cross-correlation time to an actual time
                ccc = normxcorr2(tr.data, image[0].data)
                cc_max = np.amax(ccc)
                picktime = image[0].stats.starttime + (np.argmax(ccc) *
                                                       image[0].stats.delta)
            log.debug('********DEBUG: Maximum cross-corr=%s' % cc_max)
            checksum += cc_max
            used_chans += 1
            if cc_max < min_cc:
                continue
            cccsum += cc_max
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
            else:
                phase = None
            _waveform_id = WaveformStreamID(network_code=temp_net,
                                            station_code=temp_sta,
                                            channel_code=temp_chan)
            event.picks.append(Pick(waveform_id=_waveform_id,
                                    time=picktime,
                                    method_id=ResourceIdentifier('EQcorrscan'),
                                    phase_hint=phase,
                                    creation_info='eqcorrscan.core.lag_calc',
                                    comments=[Comment(text='cc_max=%s'
                                                      % cc_max)]))
            event.resource_id = detection_id
    ccc_str = ("detect_val=%s" % cccsum)
    event.comments.append(Comment(text=ccc_str))
    if used_chans == detect_chans:
        if pre_lag_ccsum is not None and checksum - pre_lag_ccsum < -0.05:
            msg = ('lag-calc has decreased cccsum from %f to %f - '
                   'report this error' % (pre_lag_ccsum, checksum))
            raise LagCalcError(msg)
    else:
        warnings.warn('Cannot check is cccsum is better, used %i channels '
                      'for detection, but %i are used here'
                      % (detect_chans, used_chans))
    return i, event


def _day_loop(detection_streams, template, min_cc, detections, interpolate,
              cores, parallel):
    """
    Function to loop through multiple detections for one template.

    Designed to run for the same day of data for I/O simplicity, but as you
    are passing stream objects it could run for all the detections ever, as
    long as you have the RAM!

    :type detection_streams: list
    :param detection_streams:
        List of all the detections for this template that you want to compute
        the optimum pick for. Individual things in list should be of
        :class:`obspy.core.stream.Stream` type.
    :type template: obspy.core.stream.Stream
    :param template: The original template used to detect the detections passed
    :type min_cc: float
    :param min_cc: Minimum cross-correlation value to be allowed for a pick.
    :type detections: list
    :param detections:
        List of detections to associate events with an input detection.
    :type interpolate: bool
    :param interpolate:
        Interpolate the correlation function to achieve sub-sample precision.

    :returns:
        Catalog object containing Event objects for each detection created by
        this template.
    :rtype: :class:`obspy.core.event.Catalog`
    """
    if not cores:
        num_cores = cpu_count()
    else:
        num_cores = cores
    if num_cores > len(detection_streams):
        num_cores = len(detection_streams)
    if parallel:
        pool = Pool(processes=num_cores)
        log.debug('Made pool of %i workers' % num_cores)
        # Parallelize generation of events for each detection:
        # results will be a list of (i, event class)
        results = [pool.apply_async(_channel_loop, args=(detection_streams[i],
                                                         template, min_cc,
                                                         detections[i].id,
                                                         interpolate, i,
                                                         detections[i].
                                                         detect_val,
                                                         detections[i].
                                                         no_chans))
                   for i in range(len(detection_streams))]
        pool.close()
        events_list = [p.get() for p in results]
        pool.join()
        events_list.sort(key=lambda tup: tup[0])  # Sort based on index.
    else:
        events_list = []
        for i in range(len(detection_streams)):
            events_list.append(_channel_loop(detection_streams[i],
                                             template, min_cc,
                                             detections[i].id, interpolate, i,
                                             detections[i].detect_val,
                                             detections[i].no_chans))
    temp_catalog = Catalog()
    temp_catalog.events = [event_tup[1] for event_tup in events_list]
    return temp_catalog


def _prepare_data(detect_data, detections, zipped_templates, delays,
                  shift_len, plot):
    """
    Prepare data for lag_calc - reduce memory here.

    :type detect_data: obspy.core.stream.Stream
    :param detect_data: Stream to extract detection streams from.
    :type detections: list
    :param detections:
        List of :class:`eqcorrscan.core.match_filter.DETECTION` to get
        data for.
    :type zipped_templates: zip
    :param zipped_templates: Zipped list of (template_name, template)
    :type delays: list
    :param delays: List of lists of the delays for each template
    :type shift_len: float
    :param shift_len: Shift length in seconds allowed for picking.
    :type plot: bool
    :param plot:
        Whether to plot the data extracted or not, used for debugging.

    :returns: List of detect_streams to be worked on
    :rtype: list
    """
    detect_streams = []
    for detection in detections:
        # Stream to be saved for new detection
        detect_stream = []
        max_delay = 0
        for tr in detect_data:
            tr_copy = tr.copy()
            # Right now, copying each trace hundreds of times...
            template = [t for t in zipped_templates
                        if str(t[0]) == str(detection.template_name)]
            if len(template) > 0:
                template = template[0]
            else:
                warnings.warn('No template with name: %s' %
                              detection.template_name)
                for t in zipped_templates:
                    print(t)
                continue
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
            if delay > max_delay:
                max_delay = delay
            detect_stream.append(tr_copy.trim(starttime=detection.detect_time -
                                              shift_len + delay,
                                              endtime=detection.detect_time +
                                              delay + shift_len +
                                              template_len))
            del tr_copy
        for tr in detect_stream:
            if len(tr.data) == 0:
                msg = ('No data in %s.%s for detection at time %s' %
                       (tr.stats.station, tr.stats.channel,
                        detection.detect_time))
                log.debug(msg)
                warnings.warn(msg)
                detect_stream.remove(tr)
            if tr.stats.endtime - tr.stats.starttime < template_len:
                msg = ("Insufficient data for %s.%s will not use."
                       % (tr.stats.station, tr.stats.channel))
                log.debug(msg)
                warnings.warn(msg)
                detect_stream.remove(tr)
        # Check for duplicate traces
        stachans = [(tr.stats.station, tr.stats.channel)
                    for tr in detect_stream]
        c_stachans = Counter(stachans)
        for key in c_stachans.keys():
            if c_stachans[key] > 1:
                msg = ('Multiple channels for %s.%s, likely a data issue'
                       % (key[0], key[1]))
                raise LagCalcError(msg)
        if plot:
            background = detect_data.copy().trim(starttime=detection.
                                                 detect_time - (shift_len + 5),
                                                 endtime=detection.
                                                 detect_time + shift_len +
                                                 max_delay + 7)
            for tr in background:
                if len(tr.data) == 0:
                    background.remove(tr)
            detection_multiplot(stream=background,
                                template=Stream(detect_stream),
                                times=[detection.detect_time -
                                       shift_len],
                                title='Detection Extracted')
        if not len(detect_stream) == 0:
            # Create tuple of (template name, data stream)
            detect_streams.append((detection.template_name,
                                   Stream(detect_stream)))
    return detect_streams


def lag_calc(detections, detect_data, template_names, templates,
             shift_len=0.2, min_cc=0.4, cores=1, interpolate=False,
             plot=False, parallel=True, debug=0):
    """
    Main lag-calculation function for detections of specific events.

    Overseer function to take a list of detection objects, cut the data for
    them to lengths of the same length of the template + shift_len on
    either side. This will output a :class:`obspy.core.event.Catalog` of
    picked events. Pick times are based on the lag-times found at the maximum
    correlation, providing that correlation is above the min_cc.

    :type detections: list
    :param detections:
        List of :class:`eqcorrscan.core.match_filter.DETECTION` objects.
    :type detect_data: obspy.core.stream.Stream
    :param detect_data:
        All the data needed to cut from - can be a gappy Stream.
    :type template_names: list
    :param template_names:
        List of the template names, used to help identify families of events.
        Must be in the same order as templates.
    :type templates: list
    :param templates:
        List of the templates, templates should be
        :class:`obspy.core.stream.Stream` objects.
    :type shift_len: float
    :param shift_len:
        Shift length allowed for the pick in seconds, will be plus/minus this
        amount - default=0.2
    :type min_cc: float
    :param min_cc:
        Minimum cross-correlation value to be considered a pick, default=0.4.
    :type cores: int
    :param cores:
        Number of cores to use in parallel processing, defaults to one.
    :type interpolate: bool
    :param interpolate:
        Interpolate the correlation function to achieve sub-sample precision.
    :type plot: bool
    :param plot:
        To generate a plot for every detection or not, defaults to False
    :type parallel: bool
    :param parallel: Turn parallel processing on or off.
    :type debug: int
    :param debug: Debug output level, 0-5 with 5 being the most output.
    .

    :returns:
        Catalog of events with picks.  No origin information is included.
        These events can then be written out via
        :func:`obspy.core.event.Catalog.write`, or to Nordic Sfiles using
        :func:`eqcorrscan.utils.sfile_util.eventtosfile` and located
        externally.
    :rtype: obspy.core.event.Catalog

    .. note:: Picks output in catalog are generated relative to the template \
        start-time.  For example, if you generated your template with a \
        pre_pick time of 0.2 seconds, you should expect picks generated by \
        lag_calc to occur 0.2 seconds before the true phase-pick.  This \
        is because we do not currently store template meta-data alongside the \
        templates.

    .. warning:: Because of the above note, origin times will be consistently \
        shifted by the static pre_pick applied to the templates.

    .. note::
        Individual channel cross-correlations are stored as a
        :class:`obspy.core.event.Comment` for each pick, and the summed
        cross-correlation value resulting from these is stored as a
        :class:`obspy.core.event.Comment` in the main
        :class:`obspy.core.event.Event` object.

    .. note::
        The order of events is preserved (e.g. detections[n] == output[n]),
        providing picks have been made for that event.  If no picks have
        been made for an event, it will not be included in the output.
        However, as each detection has an ID associated with it, these can
        be mapped to the output resource_id for each Event in the output
        Catalog. e.g.

            detections[n].id == output[m].resource_id

        if the output[m] is for the same event as detections[n].
    """
    if debug > 0:
        log.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
    if debug > 2:
        log.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    if debug > 3:
        log.setLevel(0)
        ch.setLevel(0)
    log.addHandler(ch)
    if debug > 2 and plot:
        prep_plot = True
    else:
        prep_plot = False
    # First check that sample rates are equal for everything
    for tr in detect_data:
        if tr.stats.sampling_rate != detect_data[0].stats.sampling_rate:
            raise LagCalcError('Sampling rates are not equal')
    for template in templates:
        for tr in template:
            if tr.stats.sampling_rate != detect_data[0].stats.sampling_rate:
                raise LagCalcError('Sampling rates are not equal')
    # Work out the delays for each template
    delays = []  # List of tuples of (tempname, (sta, chan, delay))
    zipped_templates = list(zip(template_names, templates))
    detect_stachans = [(tr.stats.station, tr.stats.channel)
                       for tr in detect_data]
    for template in zipped_templates:
        temp_delays = []
        # Remove channels not present in continuous data
        _template = template[1].copy()
        for tr in _template:
            if (tr.stats.station, tr.stats.channel) not in detect_stachans:
                _template.remove(tr)
        for tr in _template:
            temp_delays.append((tr.stats.station, tr.stats.channel,
                                tr.stats.starttime - _template.
                                sort(['starttime'])[0].stats.starttime))
        delays.append((template[0], temp_delays))
        del _template
    # Segregate detections by template, then feed to day_loop
    initial_cat = Catalog()
    for template in zipped_templates:
        log.info('Running lag-calc for template %s' % template[0])
        template_detections = [detection for detection in detections
                               if detection.template_name == template[0]]
        log.info('There are %i detections' % len(template_detections))
        detect_streams = _prepare_data(detect_data=detect_data,
                                       detections=template_detections,
                                       zipped_templates=zipped_templates,
                                       delays=delays, shift_len=shift_len,
                                       plot=prep_plot)
        detect_streams = [detect_stream[1] for detect_stream in detect_streams]
        if len(template_detections) > 0:
            template_cat = _day_loop(detection_streams=detect_streams,
                                     template=template[1], min_cc=min_cc,
                                     detections=template_detections,
                                     interpolate=interpolate, cores=cores,
                                     parallel=parallel)
            initial_cat += template_cat
            if plot:
                for i, event in enumerate(template_cat):
                    if len(event.picks) == 0:
                        log.warning('Made no picks for event at time %s' %
                                    event)
                        continue
                    plot_stream = detect_streams[i].copy()
                    pick_stachans = [(pick.waveform_id.station_code,
                                      pick.waveform_id.channel_code)
                                     for pick in event.picks]
                    for tr in plot_stream:
                        if (tr.stats.station, tr.stats.channel) \
                                not in pick_stachans:
                            plot_stream.remove(tr)
                    template_plot = template[1].copy()
                    for tr in template_plot:
                        if (tr.stats.station, tr.stats.channel) \
                                not in pick_stachans:
                            template_plot.remove(tr)
                    plot_repicked(template=template_plot, picks=event.picks,
                                  det_stream=plot_stream)
    sys.stdout.flush()
    # Order the catalogue to match the input
    output_cat = Catalog()
    for det in detections:
        event = [e for e in initial_cat if e.resource_id == det.id]
        if len(event) == 1:
            output_cat.append(event[0])
        elif len(event) == 0:
            print('No picks made for detection:')
            print(det)
        else:
            raise NotImplementedError('Multiple events with same id,'
                                      ' should not happen')
    return output_cat


if __name__ == '__main__':
    import doctest
    doctest.testmod()
