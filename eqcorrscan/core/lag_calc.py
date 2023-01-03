"""
Functions to generate pick-corrections for events detected by correlation.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import scipy
import logging
import os

from collections import Counter, namedtuple

from obspy import Stream, Trace
from obspy.core.event import Catalog
from obspy.core.event import Event, Pick, WaveformStreamID
from obspy.core.event import ResourceIdentifier, Comment

from eqcorrscan.utils.correlate import get_stream_xcorr
from eqcorrscan.core.match_filter.family import Family
from eqcorrscan.core.match_filter.template import Template
from eqcorrscan.utils.plotting import plot_repicked
from eqcorrscan.utils.pre_processing import _stream_quick_select

show_interp_deprec_warning = True

Logger = logging.getLogger(__name__)


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


def _xcorr_interp(ccc, dt, resample_factor=10, use_new_resamp_method=False,
                  **kwargs):
    """
    Resample correlation-trace and check if there is a better CCC peak for
    sub-sample precision.

    :param ccc: Cross-correlation array
    :type ccc: numpy.ndarray
    :param dt: sample interval
    :type dt: float
    :param resample_factor:
        Factor for upsampling CC-values (only for use_new_resamp_method=True)
    :type resample_factor: int

    :return: Position of interpolated maximum in seconds from start of ccc
    :rtype: float
    """
    if ccc.shape[0] == 1:
        cc = ccc[0]
    else:
        cc = ccc

    # New method with resampling - make this the default in a future version
    if use_new_resamp_method:
        cc_resampled = scipy.signal.resample(cc, len(cc) * resample_factor + 1)
        dt_resampled = dt / resample_factor
        cc_t = np.arange(0, len(cc_resampled) * dt_resampled, dt_resampled)
        peak_index = cc_resampled.argmax()
        cc_peak = max(cc_resampled)

        shift = cc_t[peak_index]
        if (cc_peak < np.amax(cc) or cc_peak > 1.0 or
                not 0 < shift < len(ccc) * dt):
            # Sometimes the interpolation returns a worse result.
            Logger.warning("Interpolation did not give an accurate result, "
                           "returning maximum in data")
            return np.argmax(ccc) * dt, np.amax(ccc)
        return shift, cc_peak

    # Otherwise use old interpolation method, but warn with deprcation message
    # (but show it only once):
    global show_interp_deprec_warning
    if show_interp_deprec_warning:
        Logger.warning(
            'This method for interpolating cross-correlations is deprecated, '
            'use a more robust method with use_new_resamp_method=True')
        show_interp_deprec_warning = False
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
        Logger.warning(
            "Fewer than 3 samples selected for fit to cross correlation: "
            "{0}, returning maximum in data".format(num_samples))
        return np.argmax(cc) * dt, np.amax(cc)
    if num_samples < 5:
        Logger.debug(
            "Fewer than 5 samples selected for fit to cross correlation: "
            "{0}".format(num_samples))
    coeffs, residual = np.polyfit(
        cc_t[first_sample:last_sample + 1],
        cc[first_sample:last_sample + 1], deg=2, full=True)[:2]
    # check results of fit
    if coeffs[0] >= 0:
        Logger.info("Fitted parabola opens upwards!")
    if residual > 0.1:
        Logger.info(
            "Residual in quadratic fit to cross correlation maximum larger "
            "than 0.1: {0}".format(residual))
    # X coordinate of vertex of parabola gives time shift to correct
    # differential pick time. Y coordinate gives maximum correlation
    # coefficient.
    shift = -coeffs[1] / 2.0 / coeffs[0]
    coeff = (4 * coeffs[0] * coeffs[2] - coeffs[1] ** 2) / (4 * coeffs[0])
    if coeff < np.amax(ccc) or coeff > 1.0 or not 0 < shift < len(ccc) * dt:
        # Sometimes the interpolation returns a worse result.
        Logger.warning("Interpolation did not give an accurate result, "
                       "returning maximum in data")
        return np.argmax(ccc) * dt, np.amax(ccc)
    return shift, coeff


def _concatenate_and_correlate(streams, template, cores):
    """
    Concatenate a list of streams into one stream and correlate that with a
    template.

    All traces in a stream must have the same length.
    """
    UsedChannel = namedtuple("UsedChannel", "channel used")

    samp_rate = {tr.stats.sampling_rate for st in streams for tr in st}
    assert len(samp_rate) == 1, "Multiple sample rates found"
    samp_rate = samp_rate.pop()

    template_length = {tr.stats.npts for tr in template}
    assert len(template_length) == 1, "Multiple data lengths in template"
    template_length = template_length.pop()

    channel_length = {tr.stats.npts for st in streams for tr in st}
    if len(channel_length) > 1:
        Logger.debug("Multiple lengths of stream found, using the longest")
    channel_length = sorted(list(channel_length))[-1]
    # pre-define stream for efficiency
    chans = {tr.id for st in streams for tr in st}.intersection(
        {tr.id for tr in template})
    data = np.zeros((len(chans), channel_length * len(streams)),
                    dtype=np.float32)

    # concatenate detection streams together.
    used_chans = [[] for _ in range(len(streams))]
    concatenated_stream = Stream()
    for i, chan in enumerate(chans):
        start_index = 0
        for j, stream in enumerate(streams):
            # tr = stream.select(id=chan)
            tr = _stream_quick_select(stream, chan)
            if len(tr) == 0:
                # No data for this channel in this stream
                used_chans[j].append(UsedChannel(
                    channel=(chan.split('.')[1], chan.split('.')[-1]),
                    used=False))
                start_index += channel_length
                continue
            assert len(tr) == 1, "Multiple channels found for {0}".format(chan)
            data[i][start_index:start_index + tr[0].stats.npts] = tr[0].data
            start_index += channel_length
            used_chans[j].append(UsedChannel(
                channel=(chan.split('.')[1], chan.split('.')[-1]), used=True))
        net, sta, loc, chan = chan.split('.')
        concatenated_stream += Trace(
            data=data[i], header=dict(network=net, station=sta, channel=chan,
                                      location=loc, sampling_rate=samp_rate))
    # Remove unnecesary channels from template
    _template = Stream()
    for tr in template:
        if tr.id in chans:
            _template += tr
    # Do correlations
    xcorr_func = get_stream_xcorr(name_or_func="fftw")
    ccc, _, chan_order = xcorr_func(
        templates=[_template], stream=concatenated_stream, stack=False,
        cores=cores)
    # Re-order used_chans
    chan_order = chan_order[0]
    for i in range(len(used_chans)):
        _used_chans = used_chans[i]
        # Remove any channels that ended up not being used.
        _used_chans = [c for c in _used_chans if c.channel in chan_order]
        # Order the channels in the same way that they were correlated
        _used_chans.sort(key=lambda chan: chan_order.index(chan.channel))
        used_chans[i] = _used_chans  # Put back in.

    # Reshape ccc output
    ccc_out = np.zeros((len(streams), len(chans),
                        channel_length - template_length + 1),
                       dtype=np.float32)
    for i in range(len(streams)):
        for j, chan in enumerate(used_chans[i]):
            if not chan.used:
                continue
            index_start = i * channel_length
            index_end = index_start + channel_length - template_length + 1
            ccc_out[i][j] = ccc[0][j][index_start: index_end]
    return ccc_out, used_chans


def xcorr_pick_family(family, stream, shift_len=0.2, min_cc=0.4,
                      min_cc_from_mean_cc_factor=None,
                      all_vert=False, all_horiz=False, vertical_chans=['Z'],
                      horizontal_chans=['E', 'N', '1', '2', '3'],
                      cores=1, interpolate=False,
                      plot=False, plotdir=None, export_cc=False, cc_dir=None,
                      **kwargs):
    """
    Compute cross-correlation picks for detections in a family.

    :type family: `eqcorrscan.core.match_filter.family.Family`
    :param family: Family to calculate correlation picks for.
    :type stream: `obspy.core.stream.Stream`
    :param stream:
        Data stream containing data for all (or a subset of) detections in
        the Family
    :type shift_len: float
    :param shift_len:
        Shift length allowed for the pick in seconds, will be plus/minus this
        amount - default=0.2
    :type min_cc: float
    :param min_cc:
        Minimum cross-correlation value to be considered a pick, default=0.4.
    :type min_cc_from_mean_cc_factor: float
    :param min_cc_from_mean_cc_factor:
        If set to a value other than None, then the minimum cross-correlation
        value for a trace is set individually for each detection based on:
        min(detect_val / n_chans * min_cc_from_mean_cc_factor, min_cc).
    :type horizontal_chans: list
    :param horizontal_chans:
        List of channel endings for horizontal-channels, on which S-picks will
        be made.
    :type vertical_chans: list
    :param vertical_chans:
        List of channel endings for vertical-channels, on which P-picks will
        be made.
    :type cores: int
    :param cores:
        Number of cores to use in parallel processing, defaults to one.
    :type interpolate: bool
    :param interpolate:
        Interpolate the correlation function to achieve sub-sample precision.
    :type plot: bool
    :param plot:
        To generate a plot for every detection or not, defaults to False
    :type plotdir: str
    :param plotdir:
        Path to plotting folder, plots will be output here.
    :type export_cc: bool
    :param export_cc:
        To generate a binary file in NumPy for every detection or not,
        defaults to False
    :type cc_dir: str
    :param cc_dir:
        Path to saving folder, NumPy files will be output here.

    :return: Dictionary of picked events keyed by detection id.
    """
    picked_dict = {}
    delta = family.template.st[0].stats.delta
    detect_streams_dict = _prepare_data(
        family=family, detect_data=stream, shift_len=shift_len,
        all_vert=all_vert, all_horiz=all_horiz, vertical_chans=vertical_chans,
        horizontal_chans=horizontal_chans)
    detection_ids = list(detect_streams_dict.keys())
    detect_streams = [detect_streams_dict[detection_id]
                      for detection_id in detection_ids]
    if len(detect_streams) == 0:
        Logger.warning("No appropriate data found, check your family and "
                       "detections - make sure seed ids match")
        return picked_dict
    if len(detect_streams) != len(family):
        Logger.warning("Not all detections have matching data. "
                       "Proceeding anyway. HINT: Make sure SEED IDs match")
    # Correlation function needs a list of streams, we need to maintain order.
    ccc, chans = _concatenate_and_correlate(
        streams=detect_streams, template=family.template.st, cores=cores)
    for i, detection_id in enumerate(detection_ids):
        detection = [d for d in family.detections if d.id == detection_id][0]
        correlations = ccc[i]
        if export_cc:
            os.makedirs(cc_dir, exist_ok=True)
            fname = f"{detection_id}-cc.npy"
            np.save(os.path.join(cc_dir, f'{fname}'), correlations)
            Logger.info(f"Saved correlation statistic to {fname} (lag_calc)")
        picked_chans = chans[i]
        detect_stream = detect_streams_dict[detection_id]
        checksum, cccsum, used_chans = 0.0, 0.0, 0
        event = Event()
        if min_cc_from_mean_cc_factor is not None:
            cc_thresh = min(abs(detection.detect_val / detection.no_chans
                                * min_cc_from_mean_cc_factor),
                            min_cc)
            Logger.info('Setting minimum cc-threshold for detection %s to %s',
                        detection.id, str(cc_thresh))
        else:
            cc_thresh = min_cc
        for correlation, stachan in zip(correlations, picked_chans):
            if not stachan.used:
                continue
            tr = detect_stream.select(
                station=stachan.channel[0], channel=stachan.channel[1])[0]
            if interpolate:
                shift, cc_max = _xcorr_interp(correlation, dt=delta, **kwargs)
            else:
                cc_max = np.amax(correlation)
                shift = np.argmax(correlation) * delta
            if np.isnan(cc_max):  # pragma: no cover
                Logger.error(
                    'Problematic trace, no cross correlation possible')
                continue
            picktime = tr.stats.starttime + shift
            checksum += cc_max
            used_chans += 1
            if cc_max < cc_thresh:
                Logger.debug('Correlation of {0} is below threshold, not '
                             'using'.format(cc_max))
                continue
            cccsum += cc_max
            phase = None
            if stachan.channel[1][-1] in vertical_chans:
                phase = 'P'
            elif stachan.channel[1][-1] in horizontal_chans:
                phase = 'S'
            _waveform_id = WaveformStreamID(seed_string=tr.id)
            event.picks.append(Pick(
                waveform_id=_waveform_id, time=picktime,
                method_id=ResourceIdentifier('EQcorrscan'), phase_hint=phase,
                creation_info='eqcorrscan.core.lag_calc',
                evaluation_mode='automatic',
                comments=[Comment(text='cc_max={0}'.format(cc_max))]))
        event.resource_id = ResourceIdentifier(detection_id)
        event.comments.append(Comment(text="detect_val={0}".format(cccsum)))
        # Add template-name as comment to events
        event.comments.append(Comment(
            text="Detected using template: {0}".format(family.template.name)))
        if used_chans == detection.no_chans:  # pragma: no cover
            if detection.detect_val is not None and\
               checksum - detection.detect_val < -(0.3 * detection.detect_val):
                msg = ('lag-calc has decreased cccsum from %f to %f - '
                       % (detection.detect_val, checksum))
                Logger.error(msg)
                continue
        else:
            Logger.warning(
                'Cannot check if cccsum is better, used {0} channels for '
                'detection, but {1} are used here'.format(
                    detection.no_chans, used_chans))
        picked_dict.update({detection_id: event})
    if plot:  # pragma: no cover
        for i, event in enumerate(picked_dict.values()):
            if len(event.picks) == 0:
                continue
            plot_stream = detect_streams[i].copy()
            template_plot = family.template.st.copy()
            pick_stachans = [(pick.waveform_id.station_code,
                              pick.waveform_id.channel_code)
                             for pick in event.picks]
            for tr in plot_stream:
                if (tr.stats.station, tr.stats.channel) \
                        not in pick_stachans:
                    plot_stream.remove(tr)
            for tr in template_plot:
                if (tr.stats.station, tr.stats.channel) \
                        not in pick_stachans:
                    template_plot.remove(tr)
            if plotdir is not None:
                if not os.path.isdir(plotdir):
                    os.makedirs(plotdir)
                savefile = "{plotdir}/{rid}.png".format(
                    plotdir=plotdir, rid=event.resource_id.id)
                plot_repicked(template=template_plot, picks=event.picks,
                              det_stream=plot_stream, show=False, save=True,
                              savefile=savefile)
            else:
                plot_repicked(template=template_plot, picks=event.picks,
                              det_stream=plot_stream, show=True)
    return picked_dict


def _prepare_data(family, detect_data, shift_len, all_vert=False,
                  all_horiz=False, vertical_chans=['Z'],
                  horizontal_chans=['E', 'N', '1', '2']):
    """
    Prepare data for lag_calc - reduce memory here.

    :type family: `eqcorrscan.core.match_filter.family.Family`
    :param family:
        The Family containing the template and detections.
    :type detect_data: obspy.core.stream.Stream
    :param detect_data: Stream to extract detection streams from.
    :type shift_len: float
    :param shift_len: Shift length in seconds allowed for picking.

    :returns: Dictionary of detect_streams keyed by detection id
              to be worked on
    :rtype: dict
    """
    lengths = {tr.stats.npts * tr.stats.delta for tr in family.template.st}
    assert len(lengths) == 1, "Template contains channels of different length"
    length = lengths.pop() + (2 * shift_len)
    # Enforce length be an integer number of samples
    length_samples = length * family.template.samp_rate
    if length_samples != int(length_samples):
        length = round(length_samples) / family.template.samp_rate
        Logger.info("Setting length to {0}s to give an integer number of "
                    "samples".format(length))
    prepick = shift_len + family.template.prepick
    detect_streams_dict = family.extract_streams(
        stream=detect_data, length=length, prepick=prepick,
        all_vert=all_vert, all_horiz=all_horiz, vertical_chans=vertical_chans,
        horizontal_chans=horizontal_chans)
    for key, detect_stream in detect_streams_dict.items():
        # Split to remove trailing or leading masks
        for i in range(len(detect_stream) - 1, -1, -1):
            trace = detect_stream[i]
            if np.ma.is_masked(trace.data):
                detect_stream.remove(trace)
                Logger.warning("Masked array found for {0}, not supported, "
                               "removing.".format(trace.id))
        stachans = [(tr.stats.station, tr.stats.channel)
                    for tr in detect_stream]
        c_stachans = Counter(stachans)
        for key in c_stachans.keys():
            if c_stachans[key] > 1:
                raise LagCalcError(
                    'Multiple channels for {0}.{1}, likely a data '
                    'issue'.format(key[0], key[1]))
    detect_streams_dict_out = {}
    for key, detect_stream in detect_streams_dict.items():
        if len(detect_stream) == 0:
            Logger.info("No data for {0}".format(key))
            continue
        detect_streams_dict_out.update({key: detect_stream})
    return detect_streams_dict_out


def lag_calc(detections, detect_data, template_names, templates,
             shift_len=0.2, min_cc=0.4, min_cc_from_mean_cc_factor=None,
             all_vert=False, all_horiz=False,
             horizontal_chans=['E', 'N', '1', '2'],
             vertical_chans=['Z'], cores=1, interpolate=False,
             plot=False, plotdir=None, export_cc=False, cc_dir=None, **kwargs):
    """
    Cross-correlation derived picking of seismic events.

    Overseer function to take a list of detection objects, cut the data for
    them to lengths of the same length of the template + shift_len on
    either side. This will output a :class:`obspy.core.event.Catalog` of
    picked events. Pick times are based on the lag-times found at the maximum
    correlation, providing that correlation is above the min_cc.

    :type detections: list
    :param detections:
        List of :class:`eqcorrscan.core.match_filter.Detection` objects.
    :type detect_data: obspy.core.stream.Stream
    :param detect_data:
        All the data needed to cut from - can be a gappy Stream.
    :type template_names: list
    :param template_names:
        List of the template names, used to help identify families of events.
        Must be in the same order as templates.
    :type templates: list
    :param templates:
        List of the templates, templates must be a list of
         :class:`obspy.core.stream.Stream` objects.
    :type shift_len: float
    :param shift_len:
        Shift length allowed for the pick in seconds, will be plus/minus this
        amount - default=0.2
    :type min_cc: float
    :param min_cc:
        Minimum cross-correlation value to be considered a pick, default=0.4.
    :type min_cc_from_mean_cc_factor: float
    :param min_cc_from_mean_cc_factor:
        If set to a value other than None, then the minimum cross-correlation
        value for a trace is set individually for each detection based on:
        min(detect_val / n_chans * min_cc_from_mean_cc_factor, min_cc).
    :type horizontal_chans: list
    :param horizontal_chans:
        List of channel endings for horizontal-channels, on which S-picks will
        be made.
    :type vertical_chans: list
    :param vertical_chans:
        List of channel endings for vertical-channels, on which P-picks will
        be made.
    :type cores: int
    :param cores:
        Number of cores to use in parallel processing, defaults to one.
    :type interpolate: bool
    :param interpolate:
        Interpolate the correlation function to achieve sub-sample precision.
    :type plot: bool
    :param plot:
        To generate a plot for every detection or not, defaults to False
    :param plotdir:
        Path to plotting folder, plots will be output here.
    :type export_cc: bool
    :param export_cc:
        To generate a binary file in NumPy for every detection or not,
        defaults to False
    :type cc_dir: str
    :param cc_dir:
        Path to saving folder, NumPy files will be output here.

    :returns:
        Catalog of events with picks.  No origin information is included.
        These events can then be written out via
        :func:`obspy.core.event.Catalog.write`, or to Nordic Sfiles using
        :func:`eqcorrscan.utils.sfile_util.eventtosfile` and located
        externally.
    :rtype: obspy.core.event.Catalog

    .. note::
        Picks output in catalog are generated relative to the template
        start-time.  For example, if you generated your template with a
        pre_pick time of 0.2 seconds, you should expect picks generated by
        lag_calc to occur 0.2 seconds before the true phase-pick.  This
        is because we do not currently store template meta-data alongside the
        templates.

    .. warning::
        Because of the above note, origin times will be consistently
        shifted by the static pre_pick applied to the templates.

    .. warning::
        This routine requires only one template per channel (e.g. you should
        not use templates with a P and S template on a single channel).  If
        this does occur an error will be raised.

    .. note::
        S-picks will be made on horizontal channels, and P picks made on
        vertical channels - the default is that horizontal channels end in
        one of: 'E', 'N', '1' or '2', and that vertical channels end in 'Z'.
        The options vertical_chans and horizontal_chans can be changed to suit
        your dataset.

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

    .. note::
        The correlation data that are saved to the binary files can be useful
        to select an appropriate threshold for your data.
    """
    # First check that sample rates are equal for everything
    for tr in detect_data:
        if tr.stats.sampling_rate != detect_data[0].stats.sampling_rate:
            raise LagCalcError('Sampling rates are not equal')
    for template in templates:
        for tr in template:
            if tr.stats.sampling_rate != detect_data[0].stats.sampling_rate:
                raise LagCalcError('Sampling rates are not equal')
    initial_cat = dict()  # Dictionary keyed by detection id
    for template, template_name in zip(templates, template_names):
        Logger.info('Running lag-calc for template %s' % template[0])
        template_detections = [detection for detection in detections
                               if detection.template_name == template_name]
        for detection in template_detections:
            detection.event = detection.event or detection._calculate_event(
                template_st=template)
        family = Family(
            detections=template_detections,
            template=Template(
                name=template_name, st=template,
                samp_rate=template[0].stats.sampling_rate, prepick=0.0))
        # Make a sparse template
        if len(template_detections) > 0:
            template_dict = xcorr_pick_family(
                family=family, stream=detect_data, min_cc=min_cc,
                min_cc_from_mean_cc_factor=min_cc_from_mean_cc_factor,
                all_vert=all_vert, all_horiz=all_horiz,
                horizontal_chans=horizontal_chans,
                vertical_chans=vertical_chans, interpolate=interpolate,
                cores=cores, shift_len=shift_len, plot=plot, plotdir=plotdir,
                export_cc=export_cc, cc_dir=cc_dir, **kwargs)
            initial_cat.update(template_dict)
    # Order the catalogue to match the input
    output_cat = Catalog()
    for det in detections:
        event = initial_cat.get(det.id, None)
        if event:
            output_cat.append(event)
        else:
            Logger.error('No picks made for detection: \n{0}'.format(det))
    return output_cat


if __name__ == '__main__':
    import doctest
    doctest.testmod()
