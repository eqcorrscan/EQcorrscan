#!/usr/bin/python
"""
Functions to generate template waveforms and information to go with them
for the application of cross-correlation of seismic data for the detection of
repeating events.

.. note::
    All functions use obspy filters, which are implemented such that
    if both highcut and lowcut are set a bandpass filter will be used,
    but of highcut is not set (None) then a highpass filter will be used and
    if only the highcut is set then a lowpass filter will be used.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import numpy as np
import logging
import os
import warnings

from obspy import Stream, read, Trace, UTCDateTime, read_events
from obspy.core.event import Catalog
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.seishub import Client as SeisHubClient

from eqcorrscan.utils.sac_util import sactoevent
from eqcorrscan.utils import pre_processing
from eqcorrscan.core import EQcorrscanDeprecationWarning


Logger = logging.getLogger(__name__)


class TemplateGenError(Exception):
    """
    Default error for template generation errors.
    """
    def __init__(self, value):
        """
        Raise error.
        """
        self.value = value

    def __repr__(self):
        return self.value

    def __str__(self):
        return 'TemplateGenError: ' + self.value


def template_gen(method, lowcut, highcut, samp_rate, filt_order,
                 length, prepick, swin="all", process_len=86400,
                 all_horiz=False, delayed=True, plot=False, plotdir=None,
                 return_event=False, min_snr=None, parallel=False,
                 num_cores=False, save_progress=False, skip_short_chans=False,
                 check_full_seed=False, **kwargs):
    """
    Generate processed and cut waveforms for use as templates.

    :type method: str
    :param method:
        Template generation method, must be one of ('from_client',
        'from_seishub', 'from_sac', 'from_meta_file'). - Each method requires
        associated arguments, see note below.
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will not apply a lowcut.
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will not apply a highcut.
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz.
    :type filt_order: int
    :param filt_order: Filter level (number of corners).
    :type length: float
    :param length: Length of template waveform in seconds.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type process_len: int
    :param process_len: Length of data in seconds to download and process.
    :type all_horiz: bool
    :param all_horiz:
        To use both horizontal channels even if there is only a pick on one of
        them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Plot templates or not.
    :type plotdir: str
    :param plotdir:
        The path to save plots to. If `plotdir=None` (default) then the figure
        will be shown on screen.
    :type return_event: bool
    :param return_event: Whether to return the event and process length or not.
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.
    :type parallel: bool
    :param parallel: Whether to process data in parallel or not.
    :type num_cores: int
    :param num_cores:
        Number of cores to try and use, if False and parallel=True, will use
        either all your cores, or as many traces as in the data (whichever is
        smaller).
    :type save_progress: bool
    :param save_progress:
        Whether to save the resulting templates at every data step or not.
        Useful for long-running processes.
    :type skip_short_chans: bool
    :param skip_short_chans:
        Whether to ignore channels that have insufficient length data or not.
        Useful when the quality of data is not known, e.g. when downloading
        old, possibly triggered data from a datacentre
    :type check_full_seed: bool
    :param check_full_seed:
        If True, will check the trace header against the full SEED id,
        including Network, Station, Location and Channel. If False (default),
        will check only against Station and Channel. This behavior was
        originally necessary to cope with some software (i.e. SEISAN) not
        storing picks with full SEED info.

    :returns: List of :class:`obspy.core.stream.Stream` Templates
    :rtype: list

    .. note::
        By convention templates are generated with P-phases on the
        vertical channel and S-phases on the horizontal channels, normal
        seismograph naming conventions are assumed, where Z denotes vertical
        and N, E, R, T, 1 and 2 denote horizontal channels, either oriented
        or not.  To this end we will **only** use Z channels if they have a
        P-pick, and will use one or other horizontal channels **only** if
        there is an S-pick on it.

    .. warning::
        If there is no phase_hint included in picks, and swin=all, all channels
        with picks will be used.

    .. note::
        If swin=all, then all picks will be used, not just phase-picks (e.g. it
        will use amplitude picks). If you do not want this then we suggest
        that you remove any picks you do not want to use in your templates
        before using the event.

    .. note::
        *Method specific arguments:*

        - `from_client` requires:
            :param str client_id:
                string passable by obspy to generate Client, or any object
                with a `get_waveforms` method, including a Client instance.
            :param `obspy.core.event.Catalog` catalog:
                Catalog of events to generate template for
            :param float data_pad: Pad length for data-downloads in seconds
        - `from_seishub` requires:
            :param str url: url to seishub database
            :param `obspy.core.event.Catalog` catalog:
                Catalog of events to generate template for
            :param float data_pad: Pad length for data-downloads in seconds
        - `from_sac` requires:
            :param list sac_files:
                osbpy.core.stream.Stream of sac waveforms, or list of paths to
                sac waveforms.
            .. note::
                See `eqcorrscan.utils.sac_util.sactoevent` for details on
                how pick information is collected.
        - `from_meta_file` requires:
            :param str meta_file:
                Path to obspy-readable event file, or an obspy Catalog
            :param `obspy.core.stream.Stream` st:
                Stream containing waveform data for template. Note that this
                should be the same length of stream as you will use for the
                continuous detection, e.g. if you detect in day-long files,
                give this a day-long file!
            :param bool process:
                Whether to process the data or not, defaults to True.

    .. note::
        process_len should be set to the same length as used when computing
        detections using match_filter.match_filter, e.g. if you read
        in day-long data for match_filter, process_len should be 86400.

        .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from eqcorrscan.core.template_gen import template_gen
    >>> client = Client('NCEDC')
    >>> catalog = client.get_events(eventid='72572665', includearrivals=True)
    >>> # We are only taking two picks for this example to speed up the
    >>> # example, note that you don't have to!
    >>> catalog[0].picks = catalog[0].picks[0:2]
    >>> templates = template_gen(
    ...    method='from_client', catalog=catalog, client_id='NCEDC',
    ...    lowcut=2.0, highcut=9.0, samp_rate=20.0, filt_order=4, length=3.0,
    ...    prepick=0.15, swin='all', process_len=300, all_horiz=True)
    >>> templates[0].plot(equal_scale=False, size=(800,600)) # doctest: +SKIP

    .. figure:: ../../plots/template_gen.from_client.png

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.core.template_gen import template_gen
    >>> # Get the path to the test data
    >>> import eqcorrscan
    >>> import os
    >>> TEST_PATH = os.path.dirname(eqcorrscan.__file__) + '/tests/test_data'
    >>> st = read(TEST_PATH + '/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> quakeml = TEST_PATH + '/20130901T041115.xml'
    >>> templates = template_gen(
    ...    method='from_meta_file', meta_file=quakeml, st=st, lowcut=2.0,
    ...    highcut=9.0, samp_rate=20.0, filt_order=3, length=2, prepick=0.1,
    ...    swin='S', all_horiz=True)
    >>> print(len(templates[0]))
    10
    >>> templates = template_gen(
    ...    method='from_meta_file', meta_file=quakeml, st=st, lowcut=2.0,
    ...    highcut=9.0, samp_rate=20.0, filt_order=3, length=2, prepick=0.1,
    ...    swin='S_all', all_horiz=True)
    >>> print(len(templates[0]))
    15

    .. rubric:: Example

    >>> from eqcorrscan.core.template_gen import template_gen
    >>> import glob
    >>> # Get all the SAC-files associated with one event.
    >>> sac_files = glob.glob(TEST_PATH + '/SAC/2014p611252/*')
    >>> templates = template_gen(
    ...    method='from_sac', sac_files=sac_files, lowcut=2.0, highcut=10.0,
    ...    samp_rate=25.0, filt_order=4, length=2.0, swin='all', prepick=0.1,
    ...    all_horiz=True)
    >>> print(templates[0][0].stats.sampling_rate)
    25.0
    >>> print(len(templates[0]))
    15
    """
    if not check_full_seed:
        warnings.warn(
            "Deprecation warning: check_full_seed will default to"
            "True in a future release. Check the docs page here "
            "for how this will affect you: "
            "https://eqcorrscan.readthedocs.io/en/latest/faq.html")
    client_map = {'from_client': 'fdsn', 'from_seishub': 'seishub'}
    assert method in ('from_client', 'from_seishub', 'from_meta_file',
                      'from_sac')
    if not isinstance(swin, list):
        swin = [swin]
    process = True
    if method in ['from_client', 'from_seishub']:
        catalog = kwargs.get('catalog', Catalog())
        data_pad = kwargs.get('data_pad', 90)
        # Group catalog into days and only download the data once per day
        sub_catalogs = _group_events(
            catalog=catalog, process_len=process_len, template_length=length,
            data_pad=data_pad)
        if method == 'from_client':
            client_id = kwargs.get('client_id', None)
            if hasattr(client_id, 'get_waveforms'):
                client = client_id
            elif isinstance(client_id, str):
                client = FDSNClient(client_id)
            else:
                raise NotImplementedError(
                    "client_id must be an FDSN client string, or a Client "
                    "with a get_waveforms method"
                )
            available_stations = []
        else:
            client = SeisHubClient(kwargs.get('url', None), timeout=10)
            available_stations = client.waveform.get_station_ids()
    elif method == 'from_meta_file':
        if isinstance(kwargs.get('meta_file'), Catalog):
            catalog = kwargs.get('meta_file')
        elif kwargs.get('meta_file'):
            catalog = read_events(kwargs.get('meta_file'))
        else:
            catalog = kwargs.get('catalog')
        sub_catalogs = [catalog]
        st = kwargs.get('st', Stream())
        process = kwargs.get('process', True)
    elif method == 'from_sac':
        sac_files = kwargs.get('sac_files')
        if isinstance(sac_files, list):
            if isinstance(sac_files[0], (Stream, Trace)):
                # This is a list of streams...
                st = Stream(sac_files[0])
                for sac_file in sac_files[1:]:
                    st += sac_file
            else:
                sac_files = [read(sac_file)[0] for sac_file in sac_files]
                st = Stream(sac_files)
        else:
            st = sac_files
        # Make an event object...
        catalog = Catalog([sactoevent(st)])
        sub_catalogs = [catalog]

    temp_list = []
    process_lengths = []
    catalog_out = Catalog()

    if "P_all" in swin or "S_all" in swin or all_horiz:
        all_channels = True
    else:
        all_channels = False
    for sub_catalog in sub_catalogs:
        if method in ['from_seishub', 'from_client']:
            Logger.info("Downloading data")
            st = _download_from_client(
                client=client, client_type=client_map[method],
                catalog=sub_catalog, data_pad=data_pad,
                process_len=process_len, available_stations=available_stations,
                all_channels=all_channels)
        Logger.info('Pre-processing data')
        st.merge()
        if len(st) == 0:
            Logger.info("No data")
            continue
        if process:
            data_len = max([len(tr.data) / tr.stats.sampling_rate
                            for tr in st])
            if 80000 < data_len < 90000:
                daylong = True
                starttime = min([tr.stats.starttime for tr in st])
                min_delta = min([tr.stats.delta for tr in st])
                # Cope with the common starttime less than 1 sample before the
                #  start of day.
                if (starttime + min_delta).date > starttime.date:
                    starttime = (starttime + min_delta)
                # Check if this is stupid:
                if abs(starttime - UTCDateTime(starttime.date)) > 600:
                    daylong = False
                starttime = starttime.date
            else:
                daylong = False
            # Check if the required amount of data have been downloaded - skip
            # channels if arg set.
            for tr in st:
                if np.ma.is_masked(tr.data):
                    _len = np.ma.count(tr.data) * tr.stats.delta
                else:
                    _len = tr.stats.npts * tr.stats.delta
                if _len < process_len * .8:
                    Logger.info(
                        "Data for {0} are too short, skipping".format(
                            tr.id))
                    if skip_short_chans:
                        continue
                # Trim to enforce process-len
                tr.data = tr.data[0:int(process_len * tr.stats.sampling_rate)]
            if len(st) == 0:
                Logger.info("No data")
                continue
            if daylong:
                st = pre_processing.dayproc(
                    st=st, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, samp_rate=samp_rate,
                    parallel=parallel, starttime=UTCDateTime(starttime),
                    num_cores=num_cores)
            else:
                st = pre_processing.shortproc(
                    st=st, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, parallel=parallel,
                    samp_rate=samp_rate, num_cores=num_cores)
        data_start = min([tr.stats.starttime for tr in st])
        data_end = max([tr.stats.endtime for tr in st])

        for event in sub_catalog:
            stations, channels, st_stachans = ([], [], [])
            if len(event.picks) == 0:
                Logger.warning(
                    'No picks for event {0}'.format(event.resource_id))
                continue
            use_event = True
            # Check that the event is within the data
            for pick in event.picks:
                if not data_start < pick.time < data_end:
                    Logger.warning(
                        "Pick outside of data span: Pick time {0} Start "
                        "time {1} End time: {2}".format(
                            str(pick.time), str(data_start), str(data_end)))
                    use_event = False
            if not use_event:
                Logger.error('Event is not within data time-span')
                continue
            # Read in pick info
            Logger.debug("I have found the following picks")
            for pick in event.picks:
                if not pick.waveform_id:
                    Logger.warning(
                        'Pick not associated with waveforms, will not use:'
                        ' {0}'.format(pick))
                    continue
                Logger.debug(pick)
                stations.append(pick.waveform_id.station_code)
                channels.append(pick.waveform_id.channel_code)
            # Check to see if all picks have a corresponding waveform
            for tr in st:
                st_stachans.append('.'.join([tr.stats.station,
                                             tr.stats.channel]))
            # Cut and extract the templates
            template = _template_gen(
                event.picks, st, length, swin, prepick=prepick, plot=plot,
                all_horiz=all_horiz, delayed=delayed, min_snr=min_snr,
                plotdir=plotdir, check_full_seed=check_full_seed)
            process_lengths.append(len(st[0].data) / samp_rate)
            temp_list.append(template)
            catalog_out += event
        if save_progress:
            if not os.path.isdir("eqcorrscan_temporary_templates"):
                os.makedirs("eqcorrscan_temporary_templates")
            for template in temp_list:
                template.write(
                    "eqcorrscan_temporary_templates{0}{1}.ms".format(
                        os.path.sep, template[0].stats.starttime.strftime(
                            "%Y-%m-%dT%H%M%S")),
                    format="MSEED")
        del st
    if return_event:
        return temp_list, catalog_out, process_lengths
    return temp_list


def extract_from_stack(stack, template, length, pre_pick, pre_pad,
                       Z_include=False, pre_processed=True, samp_rate=None,
                       lowcut=None, highcut=None, filt_order=3):
    """
    Extract a multiplexed template from a stack of detections.

    Function to extract a new template from a stack of previous detections.
    Requires the stack, the template used to make the detections for the \
    stack, and we need to know if the stack has been pre-processed.

    :type stack: obspy.core.stream.Stream
    :param stack: Waveform stack from detections.  Can be of any length and \
        can have delays already included, or not.
    :type template: obspy.core.stream.Stream
    :param template: Template used to make the detections in the stack. Will \
        use the delays of this for the new template.
    :type length: float
    :param length: Length of new template in seconds
    :type pre_pick: float
    :param pre_pick: Extract additional data before the detection, seconds
    :type pre_pad: float
    :param pre_pad: Pad used in seconds when extracting the data, e.g. the \
        time before the detection extracted.  If using \
        clustering.extract_detections this half the length of the extracted \
        waveform.
    :type Z_include: bool
    :param Z_include: If True will include any Z-channels even if there is \
        no template for this channel, as long as there is a template for this \
        station at a different channel.  If this is False and Z channels are \
        included in the template Z channels will be included in the \
        new_template anyway.
    :type pre_processed: bool
    :param pre_processed: Have the data been pre-processed, if True (default) \
        then we will only cut the data here.
    :type samp_rate: float
    :param samp_rate: If pre_processed=False then this is required, desired \
        sampling rate in Hz, defaults to False.
    :type lowcut: float
    :param lowcut: If pre_processed=False then this is required, lowcut in \
        Hz, defaults to False.
    :type highcut: float
    :param highcut: If pre_processed=False then this is required, highcut in \
        Hz, defaults to False
    :type filt_order: int
    :param filt_order:
        If pre_processed=False then this is required, filter order, defaults
        to False

    :returns: Newly cut template.
    :rtype: :class:`obspy.core.stream.Stream`
    """
    new_template = stack.copy()
    # Copy the data before we trim it to keep the stack safe
    # Get the earliest time in the template as this is when the detection is
    # taken.
    mintime = min([tr.stats.starttime for tr in template])
    # Generate a list of tuples of (station, channel, delay) with delay in
    # seconds
    delays = [(tr.stats.station, tr.stats.channel[-1],
               tr.stats.starttime - mintime) for tr in template]

    #  Process the data if necessary
    if not pre_processed:
        new_template = pre_processing.shortproc(
            st=new_template, lowcut=lowcut, highcut=highcut,
            filt_order=filt_order, samp_rate=samp_rate)
    # Loop through the stack and trim!
    out = Stream()
    for tr in new_template:
        # Find the matching delay
        delay = [d[2] for d in delays if d[0] == tr.stats.station and
                 d[1] == tr.stats.channel[-1]]
        if Z_include and len(delay) == 0:
            delay = [d[2] for d in delays if d[0] == tr.stats.station]
        if len(delay) == 0:
            Logger.error("No matching template channel found for stack channel"
                         " {0}.{1}".format(tr.stats.station, tr.stats.channel))
        else:
            for d in delay:
                out += tr.copy().trim(
                    starttime=tr.stats.starttime + d + pre_pad - pre_pick,
                    endtime=tr.stats.starttime + d + pre_pad + length -
                    pre_pick)
    return out


def _download_from_client(client, client_type, catalog, data_pad, process_len,
                          available_stations=[], all_channels=False):
    """
    Internal function to handle downloading from either seishub or fdsn client
    """
    st = Stream()
    catalog = Catalog(sorted(catalog, key=lambda e: e.origins[0].time))
    all_waveform_info = []
    for event in catalog:
        for pick in event.picks:
            if not pick.waveform_id:
                Logger.warning(
                    "Pick not associated with waveforms, will not use:"
                    " {0}".format(pick))
                continue
            if all_channels:
                channel_code = pick.waveform_id.channel_code[0:2] + "?"
            else:
                channel_code = pick.waveform_id.channel_code
            if pick.waveform_id.station_code is None:
                Logger.error("No station code for pick, skipping")
                continue
            all_waveform_info.append((
                pick.waveform_id.network_code or "*",
                pick.waveform_id.station_code,
                channel_code, pick.waveform_id.location_code or "*"))
    starttime = UTCDateTime(
        catalog[0].origins[0].time - data_pad)
    endtime = starttime + process_len
    # Check that endtime is after the last event
    if not endtime > catalog[-1].origins[0].time + data_pad:
        raise TemplateGenError(
            'Events do not fit in processing window')
    all_waveform_info = sorted(list(set(all_waveform_info)))
    dropped_pick_stations = 0
    for waveform_info in all_waveform_info:
        net, sta, chan, loc = waveform_info
        if client_type == 'seishub' and sta not in available_stations:
            Logger.error("Station not found in SeisHub DB")
            dropped_pick_stations += 1
            continue
        Logger.info('Downloading for start-time: {0} end-time: {1}'.format(
            starttime, endtime))
        Logger.debug('.'.join([net, sta, loc, chan]))
        query_params = dict(
            network=net, station=sta, location=loc, channel=chan,
            starttime=starttime, endtime=endtime)
        try:
            st += client.get_waveforms(**query_params)
        except Exception as e:
            Logger.error(e)
            Logger.error('Found no data for this station: {0}'.format(
                query_params))
            dropped_pick_stations += 1
    if not st and dropped_pick_stations == len(event.picks):
        raise Exception('No data available, is the server down?')
    st.merge()
    # clients download chunks, we need to check that the data are
    # the desired length
    final_channels = []
    for tr in st:
        tr.trim(starttime, endtime)
        if len(tr.data) == (process_len * tr.stats.sampling_rate) + 1:
            tr.data = tr.data[1:len(tr.data)]
        if tr.stats.endtime - tr.stats.starttime < 0.8 * process_len:
            Logger.warning(
                "Data for {0}.{1} is {2} hours long, which is less than 80 "
                "percent of the desired length, will not use".format(
                    tr.stats.station, tr.stats.channel,
                    (tr.stats.endtime - tr.stats.starttime) / 3600))
        elif not pre_processing._check_daylong(tr):
            Logger.warning(
                "Data are mostly zeros, removing trace: {0}".format(tr.id))
        else:
            final_channels.append(tr)
    st.traces = final_channels
    return st


def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))


def _template_gen(picks, st, length, swin='all', prepick=0.05,
                  all_horiz=False, delayed=True, plot=False, min_snr=None,
                  plotdir=None, check_full_seed=False):
    """
    Master function to generate a multiplexed template for a single event.

    Function to generate a cut template as :class:`obspy.core.stream.Stream`
    from a given set of picks and data.  Should be given pre-processed
    data (downsampled and filtered).

    :type picks: list
    :param picks: Picks to extract data around, where each pick in the \
        list is an obspy.core.event.origin.Pick object.
    :type st: obspy.core.stream.Stream
    :param st: Stream to extract templates from
    :type length: float
    :param length: Length of template in seconds
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type prepick: float
    :param prepick:
        Length in seconds to extract before the pick time default is 0.05
        seconds.
    :type all_horiz: bool
    :param all_horiz:
        To use both horizontal channels even if there is only a pick on one
        of them.  Defaults to False.
    :type delayed: bool
    :param delayed:
        If True, each channel will begin relative to it's own pick-time, if
        set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot:
        To plot the template or not, default is False. Plots are saved as
        `template-starttime_template.png` and `template-starttime_noise.png`,
        where `template-starttime` is the start-time of the template
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.
    :type plotdir: str
    :param plotdir:
        The path to save plots to. If `plotdir=None` (default) then the figure
        will be shown on screen.
    :type check_full_seed: bool
    :param check_full_seed:
        If True, will check the trace header against the full SEED id,
        including Network, Station, Location and Channel. If False (default),
        will check only against Station and Channel. This behavior was
        originally necessary to cope with some software (i.e. SEISAN) not
        storing picks with full SEED info.

    :returns: Newly cut template.
    :rtype: :class:`obspy.core.stream.Stream`

    .. note::
        By convention templates are generated with P-phases on the
        vertical channel and S-phases on the horizontal channels, normal
        seismograph naming conventions are assumed, where Z denotes vertical
        and N, E, R, T, 1 and 2 denote horizontal channels, either oriented
        or not.  To this end we will **only** use Z channels if they have a
        P-pick, and will use one or other horizontal channels **only** if
        there is an S-pick on it.

    .. note::
        swin argument: Setting to `P` will return only data for channels
        with P picks, starting at the pick time (minus the prepick).
        Setting to `S` will return only data for channels with
        S picks, starting at the S-pick time (minus the prepick)
        (except if `all_horiz=True` when all horizontal channels will
        be returned if there is an S pick on one of them). Setting to `all`
        will return channels with either a P or S pick (including both
        horizontals if `all_horiz=True`) - with this option vertical channels
        will start at the P-pick (minus the prepick) and horizontal channels
        will start at the S-pick time (minus the prepick).
        `P_all` will return cut traces starting at the P-pick time for all
        channels. `S_all` will return cut traces starting at the S-pick
        time for all channels.

    .. warning::
        If there is no phase_hint included in picks, and swin=all, all
        channels with picks will be used.
    """
    from eqcorrscan.utils.plotting import pretty_template_plot as tplot
    from eqcorrscan.utils.plotting import noise_plot

    # the users picks intact.
    if not isinstance(swin, list):
        swin = [swin]
    for _swin in swin:
        assert _swin in ['P', 'all', 'S', 'P_all', 'S_all']
    picks_copy = []
    for pick in picks:
        if not pick.waveform_id:
            Logger.warning(
                "Pick not associated with waveform, will not use it: "
                "{0}".format(pick))
            continue
        if not pick.waveform_id.station_code or not \
                pick.waveform_id.channel_code:
            Logger.warning(
                "Pick not associated with a channel, will not use it:"
                " {0}".format(pick))
            continue
        picks_copy.append(pick)
    if len(picks_copy) == 0:
        return Stream()
    st_copy = Stream()
    for tr in st:
        # Check that the data can be represented by float16, and check they
        # are not all zeros
        if np.all(tr.data.astype(np.float16) == 0):
            Logger.error("Trace is all zeros at float16 level, either gain or "
                         "check. Not using in template: {0}".format(tr))
            continue
        st_copy += tr
    st = st_copy
    if len(st) == 0:
        return st
    # Get the earliest pick-time and use that if we are not using delayed.
    picks_copy.sort(key=lambda p: p.time)
    first_pick = picks_copy[0]
    if plot:
        stplot = st.slice(first_pick.time - 20,
                          first_pick.time + length + 90).copy()
        noise = stplot.copy()
    # Work out starttimes
    starttimes = []
    for _swin in swin:
        for tr in st:
            if check_full_seed:
                starttime = {'network': tr.stats.network,
                             'station': tr.stats.station,
                             'location': tr.stats.location,
                             'channel': tr.stats.channel,
                             'picks': []}
                station_picks = [pick for pick in picks_copy
                                 if pick.waveform_id.network_code ==
                                 tr.stats.network and
                                 pick.waveform_id.station_code ==
                                 tr.stats.station and
                                 pick.waveform_id.location_code ==
                                 tr.stats.location]
            else:
                Logger.warning(
                    'Not checking full SEED id compatibility between' +
                    ' picks and waveforms. Checking full net.sta.loc.chan ' +
                    'compatibility will be default behavior in future release')
                starttime = {'station': tr.stats.station,
                             'channel': tr.stats.channel, 'picks': []}
                station_picks = [pick for pick in picks_copy
                                 if pick.waveform_id.station_code ==
                                 tr.stats.station]
            if _swin == 'P_all':
                p_pick = [pick for pick in station_picks
                          if pick.phase_hint.upper()[0] == 'P']
                if len(p_pick) == 0:
                    continue
                starttime.update({'picks': p_pick})
            elif _swin == 'S_all':
                s_pick = [pick for pick in station_picks
                          if pick.phase_hint.upper()[0] == 'S']
                if len(s_pick) == 0:
                    continue
                starttime.update({'picks': s_pick})
            elif _swin == 'all':
                if all_horiz and tr.stats.channel[-1] in ['1', '2', '3',
                                                          'N', 'E']:
                    # Get all picks on horizontal channels
                    channel_pick = [
                        pick for pick in station_picks
                        if pick.waveform_id.channel_code[-1] in
                        ['1', '2', '3', 'N', 'E']]
                else:
                    channel_pick = [
                        pick for pick in station_picks
                        if pick.waveform_id.channel_code == tr.stats.channel]
                if len(channel_pick) == 0:
                    continue
                starttime.update({'picks': channel_pick})
            elif _swin == 'P':
                p_pick = [pick for pick in station_picks
                          if pick.phase_hint.upper()[0] == 'P' and
                          pick.waveform_id.channel_code == tr.stats.channel]
                if len(p_pick) == 0:
                    continue
                starttime.update({'picks': p_pick})
            elif _swin == 'S':
                if tr.stats.channel[-1] in ['Z', 'U']:
                    continue
                s_pick = [pick for pick in station_picks
                          if pick.phase_hint.upper()[0] == 'S']
                if not all_horiz:
                    s_pick = [pick for pick in s_pick
                              if pick.waveform_id.channel_code ==
                              tr.stats.channel]
                starttime.update({'picks': s_pick})
                if len(starttime['picks']) == 0:
                    continue
            if not delayed:
                starttime.update({'picks': [first_pick]})
            starttimes.append(starttime)
    # Cut the data
    st1 = Stream()
    for _starttime in starttimes:
        Logger.info(f"Working on channel {_starttime['station']}."
                    f"{_starttime['channel']}")
        if check_full_seed:
            tr = st.select(
                network=_starttime['network'],
                station=_starttime['station'],
                location=_starttime['location'],
                channel=_starttime['channel'])[0]
        else:
            tr = st.select(
                station=_starttime['station'],
                channel=_starttime['channel'])[0]
        Logger.info(f"Found Trace {tr}")
        used_tr = False
        for pick in _starttime['picks']:
            if not pick.phase_hint:
                Logger.warning(
                    "Pick for {0}.{1} has no phase hint given, you should not "
                    "use this template for cross-correlation"
                    " re-picking!".format(
                        pick.waveform_id.station_code,
                        pick.waveform_id.channel_code))
            starttime = pick.time - prepick
            Logger.debug("Cutting {0}".format(tr.id))
            noise_amp = _rms(
                tr.slice(starttime=starttime - 100, endtime=starttime).data)
            tr_cut = tr.slice(
                starttime=starttime, endtime=starttime + length,
                nearest_sample=False).copy()
            if plot:
                noise.select(
                    station=_starttime['station'],
                    channel=_starttime['channel']).trim(
                        noise[0].stats.starttime, starttime)
            if len(tr_cut.data) == 0:
                Logger.warning(
                    "No data provided for {0}.{1} starting at {2}".format(
                        tr.stats.station, tr.stats.channel, starttime))
                continue
            # Ensure that the template is the correct length
            if len(tr_cut.data) == (tr_cut.stats.sampling_rate *
                                    length) + 1:
                tr_cut.data = tr_cut.data[0:-1]
            Logger.debug(
                'Cut starttime = %s\nCut endtime %s' %
                (str(tr_cut.stats.starttime), str(tr_cut.stats.endtime)))
            if min_snr is not None and \
               max(tr_cut.data) / noise_amp < min_snr:
                Logger.warning(
                    "Signal-to-noise ratio {0} below threshold for {1}.{2}, "
                    "not using".format(
                        max(tr_cut.data) / noise_amp, tr_cut.stats.station,
                        tr_cut.stats.channel))
                continue
            st1 += tr_cut
            used_tr = True
        if not used_tr:
            Logger.warning('No pick for {0}'.format(tr.id))
    if plot and len(st1) > 0:
        plot_kwargs = dict(show=True)
        if plotdir is not None:
            if not os.path.isdir(plotdir):
                os.makedirs(plotdir)
            plot_kwargs.update(dict(show=False, save=True))
        tplot(st1, background=stplot, picks=picks_copy,
              title='Template for ' + str(st1[0].stats.starttime),
              savefile="{0}/{1}_template.png".format(
                  plotdir, st1[0].stats.starttime.strftime(
                      "%Y-%m-%dT%H%M%S")),
              **plot_kwargs)
        noise_plot(signal=st1, noise=noise,
                   savefile="{0}/{1}_noise.png".format(
                       plotdir, st1[0].stats.starttime.strftime(
                           "%Y-%m-%dT%H%M%S")),
                   **plot_kwargs)
        del stplot
    return st1


def _group_events(catalog, process_len, template_length, data_pad):
    """
    Internal function to group events into sub-catalogs based on process_len.

    :param catalog: Catalog to groups into sub-catalogs
    :type catalog: obspy.core.event.Catalog
    :param process_len: Length in seconds that data will be processed in
    :type process_len: int

    :return: List of catalogs
    :rtype: list
    """
    # case for catalog only containing one event
    if len(catalog) == 1:
        return [catalog]
    sub_catalogs = []
    # Sort catalog by date
    catalog.events = sorted(
        catalog.events,
        key=lambda e: (e.preferred_origin() or e.origins[0]).time)
    sub_catalog = Catalog([catalog[0]])
    for event in catalog[1:]:
        origin_time = (event.preferred_origin() or event.origins[0]).time
        last_pick = sorted(event.picks, key=lambda p: p.time)[-1]
        max_diff = (
            process_len - (last_pick.time - origin_time) - template_length)
        max_diff -= 2 * data_pad
        if origin_time - sub_catalog[0].origins[0].time < max_diff:
            sub_catalog.append(event)
        else:
            sub_catalogs.append(sub_catalog)
            sub_catalog = Catalog([event])
    sub_catalogs.append(sub_catalog)
    return sub_catalogs


if __name__ == "__main__":
    import doctest
    doctest.testmod()
