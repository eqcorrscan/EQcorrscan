"""
Functions to generate pick-corrections for events detected by correlation.


"""
import numpy as np
from eqcorrscan.core.match_filter import normxcorr2
import scipy
import warnings


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
    return shift

def _channel_loop(detection, template, min_cc, interpolate=False, i=0,
                  debug=0):
    """
    Inner loop for correlating and assigning picks.

    Utility function to take a stream of data for the detected event and write
    maximum correlation to absolute time as picks in an obspy.core.event.Event
    object.
    Only outputs picks for picks above min_cc.

    :type detection: obspy.core.stream.Stream
    :param detection: Stream of data for the slave event detected using \
        template.
    :type template: obspy.core.stream.Stream
    :param template: Stream of data as the template for the detection.
    :type interpolate: bool
    :param interpolate: Interpolate the correlation function to achieve \
        sub-sample precision.
    :type i: int
    :param i: Used to track which process has occurred when running in \
        parallel.

    :returns: Event object containing net, sta, chan information
    :rtype: obspy.core.event.Event
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
        if image:
            ccc = normxcorr2(tr.data, image[0].data)
            # Convert the maximum cross-correlation time to an actual time
            if debug > 3:
                print('********DEBUG: Maximum cross-corr=%s' % np.amax(ccc))
            if np.amax(ccc) > min_cc:
                if interpolate:
                    try:
                        interp_max = _xcorr_interp(ccc=ccc,
                                                   dt=image[0].stats.delta)
                    except IndexError:
                        print('Could not interpolate ccc, not smooth')
                        interp_max = np.argmax(ccc) * image[0].stats.delta
                    picktime = image[0].stats.starttime + interp_max
                else:
                    picktime = image[0].stats.starttime + (np.argmax(ccc) *
                                                           image[0].stats.delta)
            else:
                continue
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
                                    phase_hint=phase))
    return (i, event)


def _day_loop(detection_streams, template, min_cc, interpolate=False,
              cores=False, debug=0):
    """
    Function to loop through multiple detections for one template.

    Designed to run for the same day of data for I/O simplicity, but as you
    are passing stream objects it could run for all the detections ever, as
    long as you have the RAM!

    :type detection_streams: list
    :param detection_streams: List of all the detections for this template that
        you want to compute the optimum pick for. Individual things in list
        should be of obspy.core.stream.Stream type.
    :type template: obspy.core.stream.Stream
    :param template: The original template used to detect the detections passed
    :type min_cc: float
    :param min_cc: Minimum cross-correlation value to be allowed for a pick.
    :type interpolate: bool
    :param interpolate: Interpolate the correlation function to achieve \
        sub-sample precision.

    :returns: Catalog object containing Event objects for each detection
              created by this template.
    :rtype: obspy.core.event.Catalog
    """
    from multiprocessing import Pool, cpu_count
    # Used to run detections in parallel
    from obspy.core.event import Catalog
    if not cores:
        num_cores = cpu_count()
    else:
        num_cores = cores
    if num_cores > len(detection_streams):
        num_cores = len(detection_streams)
    pool = Pool(processes=num_cores)
    # Parallelize generation of events for each detection:
    # results is a list of (i, event class)
    results = [pool.apply_async(_channel_loop, args=(detection_streams[i],
                                                     template, min_cc,
                                                     interpolate, i, debug))
               for i in range(len(detection_streams))]
    pool.close()
    events_list = [p.get() for p in results]
    pool.join()
    events_list.sort(key=lambda tup: tup[0])  # Sort based on i.
    temp_catalog = Catalog()
    temp_catalog.events = [event_tup[1] for event_tup in events_list]
    return temp_catalog


def _prepare_data(detect_data, detections, zipped_templates, delays,
                  shift_len):
    """Prepare data for lag_calc - reduce memory here.

    :type detect_data: obspy.core.Stream
    :type detections: list
    :type zipped_templates: zip
    :type delays: list
    :type shift_len: float

    :returns: List of detect_streams to be worked on
    :rtype: list
    """
    from obspy import Stream
    detect_streams = []
    for detection in detections:
        # Stream to be saved for new detection
        detect_stream = []
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
            detect_stream.append(tr_copy.trim(starttime=detection.detect_time -
                                              shift_len + delay,
                                              endtime=detection.detect_time +
                                              delay + shift_len +
                                              template_len))
            del tr_copy
        for tr in detect_stream:
            if len(tr.data) == 0:
                detect_stream.remove(tr)
        if not len(detect_stream) == 0:
            # Create tuple of (template name, data stream)
            detect_streams.append((detection.template_name,
                                   Stream(detect_stream)))
    return detect_streams


def lag_calc(detections, detect_data, template_names, templates,
             shift_len=0.2, min_cc=0.4, cores=1, interpolate=False, plot=False):
    """
    Main lag-calculation function for detections of specific events.

    Overseer function to take a list of detection objects, cut the data for
    them to lengths of the same length of the template + shift_len on
    either side. This will then write out SEISAN s-file or QuakeML for the
    detections with pick times based on the lag-times found at the maximum
    correlation, providing that correlation is above the min_cc.

    :type detections: list
    :param detections: List of DETECTION objects
    :type detect_data: obspy.core.stream.Stream
    :param detect_data: All the data needed to cut from - can be a gappy Stream
    :type template_names: list
    :param template_names: List of the template names, used to help identify \
        families of events. Must be in the same order as templates.
    :type templates: list
    :param templates: List of the templates, templates are of type: \
        obspy.core.stream.Stream.
    :type shift_len: float
    :param shift_len: Shift length allowed for the pick in seconds, will be
        plus/minus this amount - default=0.2
    :type min_cc: float
    :param min_cc: Minimum cross-correlation value to be considered a pick,
        default=0.4
    :type cores: int
    :param cores: Number of cores to use in parallel processing, defaults to \
        one.
    :type interpolate: bool
    :param interpolate: Interpolate the correlation function to achieve \
        sub-sample precision.
    :type plot: bool
    :param plot: To generate a plot for every detection or not, defaults to \
        False.

    :returns: Catalog of events with picks.  No origin information is \
        included, these events can then be written out via \
        obspy.core.event functions, or to seisan Sfiles using Sfile_util \
        and located.
    :rtype: obspy.core.event.Catalog

    .. rubric: Example

    >>> from eqcorrscan.core import lag_calc

    .. note:: Picks output in catalog are generated relative to the template \
        start-time.  For example, if you generated your template with a \
        pre_pick time of 0.2 seconds, you should expect picks generated by \
        lag_calc to occur 0.2 seconds before the true phase-pick.  This \
        is because we do not currently store template meta-data alongside the \
        templates.

    .. warning:: Because of the above note, origin times will be consistently \
        shifted by the static pre_pick applied to the templates.
    """
    from obspy.core.event import Catalog
    from eqcorrscan.utils.plotting import plot_repicked

    # First work out the delays for each template
    delays = []  # List of tuples of (tempname, (sta, chan, delay))
    zipped_templates = list(zip(template_names, templates))
    for template in zipped_templates:
        temp_delays = []
        for tr in template[1]:
            temp_delays.append((tr.stats.station, tr.stats.channel,
                                tr.stats.starttime - template[1].
                                sort(['starttime'])[0].stats.starttime))
        delays.append((template[0], temp_delays))
    # Segregate detections by template, then feed to day_loop
    initial_cat = Catalog()
    for template in zipped_templates:
        template_detections = [detection for detection in detections
                               if detection.template_name == template[0]]
        detect_streams = _prepare_data(detect_data=detect_data,
                                       detections=template_detections,
                                       zipped_templates=zipped_templates,
                                       delays=delays, shift_len=shift_len)
        detect_streams = [detect_stream[1] for detect_stream in detect_streams]
        if len(template_detections) > 0:
            template_cat = _day_loop(detection_streams=detect_streams,
                                     template=template[1], min_cc=min_cc,
                                     interpolate=interpolate, cores=cores)
            initial_cat += template_cat
            if plot:
                for i, event in enumerate(template_cat):
                    if len(event.picks) == 0:
                        print('Made no picks for event')
                        print(event)
                        continue
                    plot_repicked(template=template[1], picks=event.picks,
                                  det_stream=detect_streams[i])
    return initial_cat


if __name__ == '__main__':
    import doctest
    doctest.testmod()