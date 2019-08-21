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

from multiprocessing import Pool, cpu_count
from collections import Counter

from obspy import Stream
from obspy.core.event import Catalog
from obspy.core.event import Event, Pick, WaveformStreamID
from obspy.core.event import ResourceIdentifier, Comment

from eqcorrscan.core.match_filter.family import Family
from eqcorrscan.utils.plotting import plot_repicked, detection_multiplot


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


def _xcorr_interp(ccc, dt):
    """
    Interpolate around the maximum correlation value for sub-sample precision.

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
        raise IndexError(
            "Fewer than 3 samples selected for fit to cross correlation: "
            "{0}".format(num_samples))
    if num_samples < 5:
        Logger.warning(
            "Fewer than 5 samples selected for fit to cross correlation: "
            "{0}".format(num_samples))
    coeffs, residual = scipy.polyfit(
        cc_t[first_sample:last_sample + 1],
        cc[first_sample:last_sample + 1], deg=2, full=True)[:2]
    # check results of fit
    if coeffs[0] >= 0:
        Logger.warning("Fitted parabola opens upwards!")
    if residual > 0.1:
        Logger.warning(
            "Residual in quadratic fit to cross correlation maximum larger "
            "than 0.1: {0}".format(residual))
    # X coordinate of vertex of parabola gives time shift to correct
    # differential pick time. Y coordinate gives maximum correlation
    # coefficient.
    shift = -coeffs[1] / 2.0 / coeffs[0]
    coeff = (4 * coeffs[0] * coeffs[2] - coeffs[1] ** 2) / (4 * coeffs[0])
    return shift, coeff


def xcorr_pick_family(family, stream, shift_len=0.2, min_cc=0.4,
                      horizontal_chans=['E', 'N', '1', '2'],
                      vertical_chans=['Z'], cores=1, interpolate=False,
                      plot=False):
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

    :return: Catalog of events.
    """
    picked_cat = Catalog()
    detect_streams = _prepare_data(
        detect_data=stream, detections=family.detections,
        template=family.template.st, shift_len=shift_len)
    detect_streams = [detect_stream[1] for detect_stream in detect_streams]

    # Add template-name as comment to events
    for event in picked_cat:
        event.comments.append(Comment(
            text="Detected using template: {0}".format(family.template.name)))
    if plot:
        for i, event in enumerate(picked_cat):
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
            plot_repicked(template=template_plot, picks=event.picks,
                          det_stream=plot_stream)
    return picked_cat


def _prepare_data(family, detect_data, shift_len):
    """
    Prepare data for lag_calc - reduce memory here.

    :type family: `eqcorrscan.core.match_filter.family.Family`
    :param family:
        The Family containing the template and detections.
    :type detect_data: obspy.core.stream.Stream
    :param detect_data: Stream to extract detection streams from.
    :type shift_len: float
    :param shift_len: Shift length in seconds allowed for picking.

    :returns: List of detect_streams to be worked on
    :rtype: list
    """
    detect_streams = []
    for detection in family.detections:
        # Stream to be saved for new detection
        # TODO: write or find an extract_stream_for_detection function and use that
        #  - Needs to maintain a minimum length requirement
        #  - Can return masked data, but this needs to throw that out
        #  - Cope with multiple versions of the same channel.
        #  Write as method on Detection and Family?
    return detect_streams


def lag_calc(detections, detect_data, template_names, templates,
             shift_len=0.2, min_cc=0.4, horizontal_chans=['E', 'N', '1', '2'],
             vertical_chans=['Z'], cores=1, interpolate=False,
             plot=False, parallel=True):
    """
    Main lag-calculation function for detections of specific events.

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
    :type parallel: bool
    :param parallel: Turn parallel processing on or off.

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
    """
    # First check that sample rates are equal for everything
    for tr in detect_data:
        if tr.stats.sampling_rate != detect_data[0].stats.sampling_rate:
            raise LagCalcError('Sampling rates are not equal')
    for template in templates:
        for tr in template:
            if tr.stats.sampling_rate != detect_data[0].stats.sampling_rate:
                raise LagCalcError('Sampling rates are not equal')
    initial_cat = Catalog()
    for template in templates:
        Logger.info('Running lag-calc for template %s' % template[0])
        template_detections = [detection for detection in detections
                               if detection.template_name == template[0]]
        for detection in template_detections:
            if detection.event is not None:
                continue
            # TODO: Estimate event for detection.
        family = Family(template=template, detections=template_detections)
        if len(template_detections) > 0:
            template_cat = xcorr_pick_family(
                family=family, stream=detect_data,
                min_cc=min_cc, horizontal_chans=horizontal_chans,
                vertical_chans=vertical_chans, interpolate=interpolate,
                cores=cores)
            initial_cat += template_cat
    # Order the catalogue to match the input
    output_cat = Catalog()
    for det in detections:
        event = [e for e in initial_cat if str(e.resource_id) == str(det.id)]
        if len(event) == 1:
            output_cat.append(event[0])
        elif len(event) == 0:
            Logger.error('No picks made for detection: \n{0}'.format(det))
        else:
            raise NotImplementedError('Multiple events with same id,'
                                      ' should not happen')
    return output_cat


if __name__ == '__main__':
    import doctest
    doctest.testmod()
