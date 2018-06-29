#!/usr/bin/python
"""
Functions to generate template waveforms and information to go with them
for the application of cross-correlation of seismic data for the detection of
repeating events.

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
    will use amplitude picks).  If you do not want this then we suggest that
    you remove any picks you do not want to use in your templates before using
    the event.

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
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import warnings
import copy
import os

from obspy import Stream, read, Trace, UTCDateTime, read_events
from obspy.core.event import Catalog
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.seishub import Client as SeisHubClient

from eqcorrscan.utils.debug_log import debug_print
from eqcorrscan.utils.sac_util import sactoevent
from eqcorrscan.utils import pre_processing
from eqcorrscan.core import EQcorrscanDeprecationWarning


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
                 length, prepick, swin, process_len=86400,
                 all_horiz=False, delayed=True, plot=False, debug=0,
                 return_event=False, min_snr=None, parallel=False,
                 num_cores=False, save_progress=False, **kwargs):
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
    :param length: Extract length in seconds.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type process_len: int
    :param process_len: Length of data in seconds to download and process.
    :type all_horiz: bool
    :param all_horiz: To use both horizontal channels even if there is only \
        a pick on one of them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Plot templates or not.
    :type debug: int
    :param debug: Level of debugging output, higher=more
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
        Whether to save the resulting party at every data step or not.
        Useful for long-running processes.

    :returns: List of :class:`obspy.core.stream.Stream` Templates
    :rtype: list

    .. note::
        *Method specific arguments:*

        - `from_client` requires:
            :param str client_id: string passable by obspy to generate Client
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
        - `from_meta_file` requires:
            :param str meta_file: Path to obspy-readable event file.
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
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> quakeml = 'eqcorrscan/tests/test_data/20130901T041115.xml'
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
    >>> sac_files = glob.glob('eqcorrscan/tests/test_data/SAC/2014p611252/*')
    >>> templates = template_gen(
    ...    method='from_sac', sac_files=sac_files, lowcut=2.0, highcut=10.0,
    ...    samp_rate=25.0, filt_order=4, length=2.0, swin='all', prepick=0.1,
    ...    all_horiz=True)
    >>> print(templates[0][0].stats.sampling_rate)
    25.0
    >>> print(len(templates[0]))
    15
    """
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
            client = FDSNClient(kwargs.get('client_id', None))
            available_stations = []
        else:
            client = SeisHubClient(kwargs.get('url', None), timeout=10)
            available_stations = client.waveform.get_station_ids()
    elif method == 'from_meta_file':
        if isinstance(kwargs.get('meta_file'), Catalog):
            catalog = kwargs.get('meta_file')
        else:
            catalog = read_events(kwargs.get('meta_file'))
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
        catalog = Catalog([sactoevent(st, debug=debug)])
        sub_catalogs = [catalog]

    temp_list = []
    process_lengths = []

    if "P_all" in swin or "S_all" in swin or all_horiz:
        all_channels = True
    else:
        all_channels = False
    for sub_catalog in sub_catalogs:
        if method in ['from_seishub', 'from_client']:
            debug_print("Downloading data", 1, debug)
            st = _download_from_client(
                client=client, client_type=client_map[method],
                catalog=sub_catalog, data_pad=data_pad,
                process_len=process_len, available_stations=available_stations,
                all_channels=all_channels, debug=debug)
        debug_print('Pre-processing data', 0, debug)
        st.merge()
        if process:
            data_len = max([len(tr.data) / tr.stats.sampling_rate
                            for tr in st])
            if 80000 < data_len < 90000:
                daylong = True
            else:
                daylong = False
            if daylong:
                starttime = min([tr.stats.starttime for tr in st])
                # Cope with the common starttime less than 1s before the
                #  start of day.
                if (starttime + 10).date > starttime.date:
                    starttime = (starttime + 10).date
                else:
                    starttime = starttime.date
                st = pre_processing.dayproc(
                    st=st, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, samp_rate=samp_rate, debug=debug,
                    parallel=parallel, starttime=UTCDateTime(starttime),
                    num_cores=num_cores)
            else:
                st = pre_processing.shortproc(
                    st=st, lowcut=lowcut, highcut=highcut,
                    filt_order=filt_order, parallel=parallel,
                    samp_rate=samp_rate, debug=debug, num_cores=num_cores)
        data_start = min([tr.stats.starttime for tr in st])
        data_end = max([tr.stats.endtime for tr in st])

        for event in sub_catalog:
            stations, channels, st_stachans = ([], [], [])
            if len(event.picks) == 0:
                debug_print('No picks for event {0}'.format(event.resource_id),
                            2, debug)
                continue
            use_event = True
            # Check that the event is within the data
            for pick in event.picks:
                if not data_start < pick.time < data_end:
                    debug_print("Pick outside of data span:\nPick time %s\n"
                                "Start time %s\nEnd time: %s" %
                                (str(pick.time), str(data_start),
                                 str(data_end)), 0, debug)
                    use_event = False
            if not use_event:
                debug_print('Event is not within data time-span', 2, debug)
                continue
            # Read in pick info
            debug_print("I have found the following picks", 0, debug)
            for pick in event.picks:
                if not pick.waveform_id:
                    debug_print(
                        'Pick not associated with waveforms, will not use:'
                        ' {0}'.format(pick), 1, debug)
                    continue
                debug_print(pick, 0, debug)
                stations.append(pick.waveform_id.station_code)
                channels.append(pick.waveform_id.channel_code)
            # Check to see if all picks have a corresponding waveform
            for tr in st:
                st_stachans.append('.'.join([tr.stats.station,
                                             tr.stats.channel]))
            # Cut and extract the templates
            template = _template_gen(
                event.picks, st, length, swin, prepick=prepick, plot=plot,
                debug=debug, all_horiz=all_horiz, delayed=delayed,
                min_snr=min_snr)
            process_lengths.append(len(st[0].data) / samp_rate)
            temp_list.append(template)
        if save_progress:
            if not os.path.isdir("eqcorrscan_temporary_templates"):
                os.makedirs("eqcorrscan_temporary_templates")
            for template in temp_list:
                template.write(
                    "eqcorrscan_temporary_templates{0}{1}.ms".format(
                        os.path.sep, template[0].stats.starttime),
                    format="MSEED")
        del st
    if return_event:
        return temp_list, catalog, process_lengths
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
    :param filt_order: If pre_processed=False then this is required, filter \
        order, defaults to False

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
            filt_order=filt_order, samp_rate=samp_rate, debug=0)
    # Loop through the stack and trim!
    out = Stream()
    for tr in new_template:
        # Find the matching delay
        delay = [d[2] for d in delays if d[0] == tr.stats.station and
                 d[1] == tr.stats.channel[-1]]
        if Z_include and len(delay) == 0:
            delay = [d[2] for d in delays if d[0] == tr.stats.station]
        if len(delay) == 0:
            debug_print("No matching template channel found for stack channel"
                        " {0}.{1}".format(tr.stats.station, tr.stats.channel),
                        2, 3)
            new_template.remove(tr)
        else:
            for d in delay:
                out += tr.copy().trim(
                    starttime=tr.stats.starttime + d + pre_pad - pre_pick,
                    endtime=tr.stats.starttime + d + pre_pad + length -
                    pre_pick)
    return out


def _download_from_client(client, client_type, catalog, data_pad, process_len,
                          available_stations=[], all_channels=False, debug=0):
    """
    Internal function to handle downloading from either seishub or fdsn client
    """
    st = Stream()
    catalog = Catalog(sorted(catalog, key=lambda e: e.origins[0].time))
    all_waveform_info = []
    for event in catalog:
        for pick in event.picks:
            if not pick.waveform_id:
                debug_print(
                    "Pick not associated with waveforms, will not use:"
                    " {0}".format(pick), 1, debug)
                continue
            if all_channels:
                channel_code = pick.waveform_id.channel_code[0:2] + "?"
            else:
                channel_code = pick.waveform_id.channel_code
            all_waveform_info.append((
                pick.waveform_id.network_code, pick.waveform_id.station_code,
                channel_code, pick.waveform_id.location_code))
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
            debug_print("Station not found in SeisHub DB", 0, debug)
            dropped_pick_stations += 1
            continue
        debug_print('Downloading for \n\tstart-time: %s\n\tend-time: %s' %
                    (str(starttime), str(endtime)), 0, debug)
        debug_print('.'.join([net, sta, loc, chan]), 0, debug)
        try:
            st += client.get_waveforms(net, sta, loc, chan,
                                       starttime, endtime)
        except Exception:
            debug_print('Found no data for this station', 2, debug)
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
            debug_print(
                "Data for {0}.{1} is {2} hours long, which is less than 80 "
                "percent of the desired length, will not pad".format(
                    tr.stats.station, tr.stats.channel,
                    (tr.stats.endtime - tr.stats.starttime) / 3600), 4, debug)
        elif not pre_processing._check_daylong(tr):
            print("Data are mostly zeros, removing trace: {0}".format(tr.id))
        else:
            final_channels.append(tr)
    st.traces = final_channels
    return st


def _template_gen(picks, st, length, swin='all', prepick=0.05,
                  all_horiz=False, delayed=True, plot=False, min_snr=None,
                  debug=0):
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
    :param prepick: Length in seconds to extract before the pick time \
            default is 0.05 seconds
    :type all_horiz: bool
    :param all_horiz: To use both horizontal channels even if there is only \
        a pick on one of them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: To plot the template or not, default is True
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.
    :type debug: int
    :param debug: Debug output level from 0-5.

    :returns: Newly cut template.
    :rtype: :class:`obspy.core.stream.Stream`

    .. note:: By convention templates are generated with P-phases on the \
        vertical channel and S-phases on the horizontal channels, normal \
        seismograph naming conventions are assumed, where Z denotes vertical \
        and N, E, R, T, 1 and 2 denote horizontal channels, either oriented \
        or not.  To this end we will **only** use Z channels if they have a \
        P-pick, and will use one or other horizontal channels **only** if \
        there is an S-pick on it.

    .. note:: swin argument: Setting to `P` will return only data for channels
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

    .. warning:: If there is no phase_hint included in picks, and swin=all, \
        all channels with picks will be used.
    """
    from eqcorrscan.utils.debug_log import debug_print
    from eqcorrscan.utils.plotting import pretty_template_plot as tplot
    from eqcorrscan.core.bright_lights import _rms
    picks_copy = copy.deepcopy(picks)  # Work on a copy of the picks and leave
    # the users picks intact.
    if not isinstance(swin, list):
        swin = [swin]
    for _swin in swin:
        assert _swin in ['P', 'all', 'S', 'P_all', 'S_all']
    for pick in picks_copy:
        if not pick.waveform_id:
            debug_print(
                "Pick not associated with waveform, will not use it: "
                "{0}".format(pick), 1, debug)
            picks_copy.remove(pick)
            continue
        if not pick.waveform_id.station_code or not \
                pick.waveform_id.channel_code:
            debug_print(
                "Pick not associated with a channel, will not use it:"
                " {0}".format(pick), 1, debug)
            picks_copy.remove(pick)
            continue
    for tr in st:
        # Check that the data can be represented by float16, and check they
        # are not all zeros
        if np.all(tr.data.astype(np.float16) == 0):
            debug_print("Trace is all zeros at float16 level, either gain or "
                        "check. Not using in template: {0}".format(tr), 4,
                        debug)
            st.remove(tr)
    # Get the earliest pick-time and use that if we are not using delayed.
    picks_copy.sort(key=lambda p: p.time)
    first_pick = picks_copy[0]
    if plot:
        stplot = st.slice(first_pick.time - 20,
                          first_pick.time + length + 90).copy()
    # Work out starttimes
    starttimes = []
    for _swin in swin:
        for tr in st:
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
        debug_print("Working on channel %s.%s" %
                    (_starttime['station'], _starttime['channel']),
                    debug_level=0, print_level=debug)
        tr = st.select(
            station=_starttime['station'], channel=_starttime['channel'])[0]
        debug_print("Found Trace %s" % tr.__str__(), debug_level=0,
                    print_level=debug)
        noise_amp = _rms(tr.data)
        used_tr = False
        for pick in _starttime['picks']:
            if not pick.phase_hint:
                warnings.warn(
                    "Pick for {0}.{1} has no phase hint given, you should not "
                    "use this template for cross-correlation"
                    " re-picking!".format(
                        pick.waveform_id.station_code,
                        pick.waveform_id.channel_code))
            starttime = pick.time - prepick
            debug_print(
                "Cutting " + tr.stats.station + '.' + tr.stats.channel, 0,
                debug)
            tr_cut = tr.slice(
                starttime=starttime, endtime=starttime + length,
                nearest_sample=False).copy()
            if len(tr_cut.data) == 0:
                debug_print(
                    "No data provided for {0}.{1} starting at {2}".format(
                        tr.stats.station, tr.stats.channel, starttime), 3,
                    debug)
                continue
            # Ensure that the template is the correct length
            if len(tr_cut.data) == (tr_cut.stats.sampling_rate *
                                    length) + 1:
                tr_cut.data = tr_cut.data[0:-1]
            debug_print(
                'Cut starttime = %s\nCut endtime %s' %
                (str(tr_cut.stats.starttime), str(tr_cut.stats.endtime)), 0,
                debug)
            if min_snr is not None and \
               max(tr_cut.data) / noise_amp < min_snr:
                debug_print(
                    "Signal-to-noise ratio below threshold for {0}.{1}".format(
                        tr_cut.stats.station, tr_cut.stats.channel), 3, debug)
                continue
            st1 += tr_cut
            used_tr = True
        if not used_tr:
            debug_print('No pick for ' + tr.stats.station + '.' +
                        tr.stats.channel, 0, debug)
    if plot:
        tplot(st1, background=stplot, picks=picks_copy,
              title='Template for ' + str(st1[0].stats.starttime))
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


# TODO: Remove these depreciated functions
def multi_template_gen(catalog, st, length, swin='all', prepick=0.05,
                       all_horiz=False, delayed=True, plot=False, debug=0,
                       return_event=False, min_snr=None):
    """
    Generate multiple templates from one stream of data.

    Thin wrapper around _template_gen to generate multiple templates from
    one stream of continuous data.  Takes processed (filtered and resampled)
    seismic data!

    :type catalog: obspy.core.event.Catalog
    :param catalog: Events to extract templates for
    :type st: obspy.core.stream.Stream
    :param st:
        Processed stream to extract from, e.g. filtered and re-sampled to what
        you want using pre_processing.dayproc.
    :type length: float
    :param length: Length of template in seconds
    :type swin: string
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type prepick: float
    :param prepick:
        Length in seconds to extract before the pick time default is
        0.05 seconds.
    :type all_horiz: bool
    :param all_horiz:
        To use both horizontal channels even if there is only a pick on one of
        them.  Defaults to False.
    :type delayed: bool
    :param delayed:
        If True, each channel will begin relative to it's own pick-time, if set
         to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: To plot the template or not, default is True
    :type debug: int
    :param debug: Debug output level from 0-5.
    :type return_event: bool
    :param return_event: Whether to return the event and process length or not.
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.

    :returns: List of :class:`obspy.core.stream.Stream` templates.
    :rtype: list

    .. warning::
        Data must be processed before using this function - highcut, lowcut and
        filt_order are only used to generate the meta-data for the templates.

    .. note:: By convention templates are generated with P-phases on the \
        vertical channel and S-phases on the horizontal channels, normal \
        seismograph naming conventions are assumed, where Z denotes vertical \
        and N, E, R, T, 1 and 2 denote horizontal channels, either oriented \
        or not.  To this end we will **only** use Z channels if they have a \
        P-pick, and will use one or other horizontal channels **only** if \
        there is an S-pick on it.

    .. warning:: If there is no phase_hint included in picks, and swin=all, \
        all channels with picks will be used.
    """
    EQcorrscanDeprecationWarning(
        "Function is depreciated and will be removed soon. Use "
        "template_gen.template_gen instead.")
    temp_list = template_gen(
        method="from_meta_file", process=False, meta_file=catalog, st=st,
        lowcut=None, highcut=None, samp_rate=st[0].stats.sampling_rate,
        filt_order=None, length=length, prepick=prepick,
        swin=swin, all_horiz=all_horiz, delayed=delayed, plot=plot,
        debug=debug, return_event=return_event, min_snr=min_snr,
        parallel=False)
    return temp_list


def from_client(catalog, client_id, lowcut, highcut, samp_rate, filt_order,
                length, prepick, swin, process_len=86400, data_pad=90,
                all_horiz=False, delayed=True, plot=False, debug=0,
                return_event=False, min_snr=None):
    """
    Generate multiplexed template from FDSN client.

    Function to generate templates from an FDSN client. Must be given \
    an obspy.Catalog class and the client_id as input. The function returns \
    a list of obspy.Stream classes containing steams for each desired \
    template.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Catalog class containing desired template events
    :type client_id: str
    :param client_id: Name of the client, either url, or Obspy \
        mappable (see the :mod:`obspy.clients.fdsn` documentation).
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will not apply a lowcut.
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will not apply a highcut.
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz.
    :type filt_order: int
    :param filt_order: Filter level (number of corners).
    :type length: float
    :param length: Extract length in seconds.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type process_len: int
    :param process_len: Length of data in seconds to download and process.
    :param data_pad: Length of data (in seconds) required before and after \
        any event for processing, use to reduce edge-effects of filtering on \
        the templates.
    :type data_pad: int
    :type all_horiz: bool
    :param all_horiz: To use both horizontal channels even if there is only \
        a pick on one of them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Plot templates or not.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type return_event: bool
    :param return_event: Whether to return the event and process length or not.
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.

    :returns: List of :class:`obspy.core.stream.Stream` Templates
    :rtype: list

    .. warning::
        This function is depreciated and will be removed in a forthcoming
        release. Please use `template_gen` instead.

    .. note::
        process_len should be set to the same length as used when computing
        detections using match_filter.match_filter, e.g. if you read
        in day-long data for match_filter, process_len should be 86400.

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from eqcorrscan.core.template_gen import from_client
    >>> client = Client('NCEDC')
    >>> catalog = client.get_events(eventid='72572665', includearrivals=True)
    >>> # We are only taking two picks for this example to speed up the
    >>> # example, note that you don't have to!
    >>> catalog[0].picks = catalog[0].picks[0:2]
    >>> templates = from_client(catalog=catalog, client_id='NCEDC',
    ...                         lowcut=2.0, highcut=9.0, samp_rate=20.0,
    ...                         filt_order=4, length=3.0, prepick=0.15,
    ...                         swin='all', process_len=300,
    ...                         all_horiz=True)
    >>> templates[0].plot(equal_scale=False, size=(800,600)) # doctest: +SKIP

    .. figure:: ../../plots/template_gen.from_client.png
    """
    EQcorrscanDeprecationWarning(
        "Function is depreciated and will be removed soon. Use "
        "template_gen.template_gen instead.")
    temp_list = template_gen(
        method="from_client", catalog=catalog, client_id=client_id,
        lowcut=lowcut, highcut=highcut, samp_rate=samp_rate,
        filt_order=filt_order, length=length, prepick=prepick,
        swin=swin, process_len=process_len, data_pad=data_pad,
        all_horiz=all_horiz, delayed=delayed, plot=plot, debug=debug,
        return_event=return_event, min_snr=min_snr)
    return temp_list


def from_seishub(catalog, url, lowcut, highcut, samp_rate, filt_order,
                 length, prepick, swin, process_len=86400, data_pad=90,
                 all_horiz=False, delayed=True, plot=False, debug=0,
                 return_event=False, min_snr=None):
    """
    Generate multiplexed template from SeisHub database.

    Function to generate templates from a SeisHub database. Must be given
    an obspy.Catalog class and the SeisHub url as input. The function returns
    a list of obspy.Stream classes containting steams for each desired
    template.

    :type catalog: :class:`obspy.core.event.Catalog`
    :param catalog: Catalog class containing desired template events
    :type url: str
    :param url: url of SeisHub database instance
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will not apply a lowcut.
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will not apply a highcut.
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz.
    :type filt_order: int
    :param filt_order: Filter level (number of corners).
    :type length: float
    :param length: Extract length in seconds.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type process_len: int
    :param process_len: Length of data in seconds to download and process.
    :param data_pad:
        Length of data (in seconds) required before and after any event for
        processing, use to reduce edge-effects of filtering on the templates.
    :type data_pad: int
    :type all_horiz: bool
    :param all_horiz:
        To use both horizontal channels even if there is only a pick on one
        of them.  Defaults to False.
    :type delayed: bool
    :param delayed:
        If True, each channel will begin relative to it's own pick-time, if
        set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Plot templates or not.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type return_event: bool
    :param return_event: Whether to return the event and process length or not.
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.

    :returns: List of :class:`obspy.core.stream.Stream` of newly cut templates
    :rtype: list

    .. note::
        process_len should be set to the same length as used when computing
        detections using match_filter.match_filter, e.g. if you read
        in day-long data fro match_filter, process_len should be 86400.

    .. warning::
        Not tested in continuous integration (due to lack of seishub client),
        let us know of any failures.
    """
    EQcorrscanDeprecationWarning(
        "Function is depreciated and will be removed soon. Use "
        "template_gen.template_gen instead.")
    temp_list = template_gen(
        method="from_seishub", catalog=catalog, url=url,
        lowcut=lowcut, highcut=highcut, samp_rate=samp_rate,
        filt_order=filt_order, length=length, prepick=prepick,
        swin=swin, process_len=process_len, data_pad=data_pad,
        all_horiz=all_horiz, delayed=delayed, plot=plot, debug=debug,
        return_event=return_event, min_snr=min_snr)
    return temp_list


def from_sac(sac_files, lowcut, highcut, samp_rate, filt_order, length, swin,
             prepick, all_horiz=False, delayed=True, plot=False, debug=0,
             return_event=False, min_snr=None):
    """
    Generate a multiplexed template from a list of SAC files.

    Function to read picks and waveforms from SAC data, and generate a \
    template from these. Usually sac_files is a list of all single-channel \
    SAC files for a given event, a single, multi-channel template will be \
    created from these traces.

    **All files listed in sac_files should be associated with a single event.**

    :type sac_files: list
    :param sac_files: osbpy.core.stream.Stream of sac waveforms, or
        list of paths to sac waveforms.
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will not apply a lowcut.
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will not apply a highcut.
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz.
    :type filt_order: int
    :param filt_order: Filter level.
    :type length: float
    :param length: Extract length in seconds.
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type prepick: float
    :param prepick: Length to extract prior to the pick in seconds.
    :type all_horiz: bool
    :param all_horiz: To use both horizontal channels even if there is only \
        a pick on one of them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Turns template plotting on or off.
    :type debug: int
    :param debug: Debug level, higher number=more output.
    :type return_event: bool
    :param return_event: Whether to return the event and process length or not.
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.

    :returns: Newly cut template.
    :rtype: :class:`obspy.core.stream.Stream`

    .. note:: This functionality is not supported for obspy versions below \
        1.0.0 as references times are not read in by SACIO, which are needed \
        for defining pick times.

    .. rubric:: Example

    >>> from eqcorrscan.core.template_gen import from_sac
    >>> import glob
    >>> # Get all the SAC-files associated with one event.
    >>> sac_files = glob.glob('eqcorrscan/tests/test_data/SAC/2014p611252/*')
    >>> templates = from_sac(sac_files=sac_files, lowcut=2.0, highcut=10.0,
    ...                      samp_rate=25.0, filt_order=4, length=2.0,
    ...                      swin='all', prepick=0.1, all_horiz=True)
    >>> print(templates[0][0].stats.sampling_rate)
    25.0
    >>> print(len(templates[0]))
    15
    """
    EQcorrscanDeprecationWarning(
        "Function is depreciated and will be removed soon. Use "
        "template_gen.template_gen instead.")
    temp_list = template_gen(
        method="from_sac", sac_files=sac_files,
        lowcut=lowcut, highcut=highcut, samp_rate=samp_rate,
        filt_order=filt_order, length=length, prepick=prepick,
        swin=swin, all_horiz=all_horiz, delayed=delayed, plot=plot,
        debug=debug, return_event=return_event, min_snr=min_snr,
        parallel=False)
    return temp_list


def from_meta_file(meta_file, st, lowcut, highcut, samp_rate, filt_order,
                   length, prepick, swin, all_horiz=False, delayed=True,
                   plot=False, parallel=True, debug=0, return_event=False,
                   min_snr=None):
    """
    Generate a multiplexed template from a local file.

    Function to generate a template from a local observation file
    and an obspy.Stream object.

    :type meta_file: str
    :param meta_file: File containing pick information, can contain \
        multiple events.  File must be formatted in a way readable by \
        :func:`obspy.core.event.read_events`.
    :type st: obspy.core.stream.Stream
    :param st: Stream containing waveform data for template (hopefully). \
        Note that this should be the same length of stream as you will use \
        for the continuous detection, e.g. if you detect in day-long files, \
        give this a day-long file!
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will not apply a lowcut.
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will not apply a highcut.
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz.
    :type filt_order: int
    :param filt_order: Filter level (number of corners).
    :type length: float
    :param length: Extract length in seconds.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin:
        P, S, P_all, S_all or all, defaults to all: see note in
        :func:`eqcorrscan.core.template_gen.template_gen`
    :type all_horiz: bool
    :param all_horiz: To use both horizontal channels even if there is only \
        a pick on one of them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Display template plots or not
    :type parallel: bool
    :param parallel: Whether to process data in parallel or not.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type return_event: bool
    :param return_event: Whether to return the event and process length or not.
    :type min_snr: float
    :param min_snr:
        Minimum signal-to-noise ratio for a channel to be included in the
        template, where signal-to-noise ratio is calculated as the ratio of
        the maximum amplitude in the template window to the rms amplitude in
        the whole window given.

    :returns: List of :class:`obspy.core.stream.Stream` newly cut templates
    :rtype: list

    .. Note::
        All picks must be associated with a station and channel, this is
        not the case for NonLinLoc HYP files, will not use any picks that
        do not have this association.

    .. warning:: We suggest giving this function a full day of data, to \
        ensure templates are generated with **exactly** the same processing \
        as the continuous data.  Not doing this will result in slightly \
        reduced cross-correlation values.

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.core.template_gen import from_meta_file
    >>> st = read('eqcorrscan/tests/test_data/WAV/TEST_/' +
    ...           '2013-09-01-0410-35.DFDPC_024_00')
    >>> quakeml = 'eqcorrscan/tests/test_data/20130901T041115.xml'
    >>> templates = from_meta_file(meta_file=quakeml, st=st, lowcut=2.0,
    ...                            highcut=9.0, samp_rate=20.0, filt_order=3,
    ...                            length=2, prepick=0.1, swin='S',
    ...                            all_horiz=True)
    >>> print(len(templates[0]))
    10
    """
    EQcorrscanDeprecationWarning(
        "Function is depreciated and will be removed soon. Use "
        "template_gen.template_gen instead.")
    temp_list = template_gen(
        method="from_meta_file", meta_file=meta_file, st=st,
        lowcut=lowcut, highcut=highcut, samp_rate=samp_rate,
        filt_order=filt_order, length=length, prepick=prepick,
        swin=swin, all_horiz=all_horiz, delayed=delayed, plot=plot,
        debug=debug, return_event=return_event, min_snr=min_snr,
        parallel=parallel)
    return temp_list


if __name__ == "__main__":
    import doctest
    doctest.testmod()
