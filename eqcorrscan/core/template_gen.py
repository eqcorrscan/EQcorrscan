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

import warnings
import numpy as np
import copy
import os
import glob

from obspy import Stream, read, Trace, UTCDateTime, read_events

from eqcorrscan.utils.debug_log import debug_print
from eqcorrscan.utils.sac_util import sactoevent
from eqcorrscan.utils import pre_processing, sfile_util


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
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
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
    >>> template = from_sac(sac_files=sac_files, lowcut=2.0, highcut=10.0,
    ...                     samp_rate=25.0, filt_order=4, length=2.0,
    ...                     swin='all', prepick=0.1, all_horiz=True)
    >>> print(template[0].stats.sampling_rate)
    25.0
    >>> print(len(template))
    15
    """
    # Check whether sac_files is a stream or a list
    if isinstance(sac_files, list):
        if isinstance(sac_files[0], Stream) or isinstance(sac_files[0], Trace):
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
    event = sactoevent(st, debug=debug)
    # Process the data
    st.merge(fill_value='interpolate')
    st = pre_processing.shortproc(
        st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
        samp_rate=samp_rate, debug=debug)
    template = template_gen(
        picks=event.picks, st=st, length=length, swin=swin, prepick=prepick,
        plot=plot, debug=debug, all_horiz=all_horiz, delayed=delayed,
        min_snr=min_snr)
    if return_event:
        return template, event, len(st[0].data) / samp_rate
    return template


def from_sfile(sfile, lowcut, highcut, samp_rate, filt_order, length, swin,
               prepick, all_horiz=False, delayed=True, plot=False, debug=0,
               return_event=False, min_snr=None):
    """
    Generate multiplexed template from a Nordic (Seisan) s-file.

    Function to read in picks from sfile then generate the template from \
    the picks within this and the wavefile found in the pick file.

    :type sfile: str
    :param sfile: sfilename must be the \
        path to a seisan nordic type s-file containing waveform and pick \
        information.
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
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
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

    .. warning:: This will use whatever data is pointed to in the s-file, if \
        this is not the coninuous data, we recommend using other functions. \
        Differences in processing between short files and day-long files \
        (inherent to resampling) will produce lower cross-correlations.

    .. rubric:: Example

    >>> from eqcorrscan.core.template_gen import from_sfile
    >>> import os
    >>> sfile = os.path.join('eqcorrscan', 'tests', 'test_data',
    ...                      'REA', 'TEST_', '01-0411-15L.S201309')
    >>> template = from_sfile(sfile=sfile, lowcut=5.0, highcut=15.0,
    ...                       samp_rate=50.0, filt_order=4, swin='P',
    ...                       prepick=0.2, length=6, all_horiz=True)
    >>> print(len(template))
    15
    >>> print(template[0].stats.sampling_rate)
    50.0
    >>> template.plot(equal_scale=False, size=(800,600)) # doctest: +SKIP

    .. plot::

        from eqcorrscan.core.template_gen import from_sfile
        import os
        sfile = os.path.realpath('../../..') + os.sep +\
            os.path.join('tests', 'test_data', 'REA',
                         'TEST_', '01-0411-15L.S201309')
        template = from_sfile(sfile=sfile, lowcut=5.0, highcut=15.0,
                              samp_rate=50.0, filt_order=4, swin='P',
                              prepick=0.2, length=6)
        template.plot(equal_scale=False, size=(800, 600))
    """
    wavefiles = sfile_util.readwavename(sfile)
    pathparts = sfile.split(os.sep)[0:-1]
    new_path_parts = []
    for part in pathparts:
        if part == 'REA':
            part = 'WAV'
        new_path_parts.append(part)
    main_wav_parts = []
    for part in new_path_parts:
        main_wav_parts.append(part)
        if part == 'WAV':
            break
    if len(main_wav_parts[0]) == 2 and main_wav_parts[0][-1] == ':':
        # Replace
        main_wav_parts[1] = main_wav_parts[0] + os.sep + main_wav_parts[1]
        new_path_parts[1] = new_path_parts[0] + os.sep + new_path_parts[1]
        main_wav_parts.remove(main_wav_parts[0])
        new_path_parts.remove(new_path_parts[0])
    mainwav = os.path.join(*main_wav_parts) + os.path.sep
    wavpath = os.path.join(*new_path_parts) + os.path.sep
    # In case of absolute paths (not handled with .split() --> .join())
    if sfile[0] == os.path.sep:
        wavpath = os.path.sep + wavpath
        mainwav = os.path.sep + mainwav
    # Read in waveform file
    st = Stream()
    for wavefile in wavefiles:
        debug_print(''.join(["I am going to read waveform data from: ",
                             wavpath, wavefile]), 0, debug)
        if os.path.isfile(wavpath + wavefile):
            # Read from a given directory
            st += read(wavpath + wavefile)
        elif os.path.isfile(wavefile):
            # Read from local directory
            st += read(wavefile)
        else:
            # Try to read from main waveform directory
            st += read(mainwav + wavefile)
    for tr in st:
        if tr.stats.sampling_rate < samp_rate:
            msg = ' '.join(['Will not upsample: Trace:', tr.stats.station,
                            'sampling rate:', str(tr.stats.sampling_rate)])
            raise TemplateGenError(msg)
    # Read in pick info
    event = sfile_util.readpicks(sfile)
    # Read the list of Picks for this event
    picks = event.picks
    if debug > 0:
        print("I have found the following picks")
        for pick in picks:
            print(pick)
    # Process waveform data
    st.merge(fill_value='interpolate')
    st = pre_processing.shortproc(
        st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
        samp_rate=samp_rate, debug=debug, seisan_chan_names=True)
    template = template_gen(
        picks=picks, st=st, length=length, swin=swin, prepick=prepick,
        all_horiz=all_horiz, plot=plot, debug=debug, delayed=delayed,
        min_snr=min_snr)
    if return_event:
        return template, event, len(st[0].data) / samp_rate
    return template


def from_contbase(sfile, contbase_list, lowcut, highcut, samp_rate, filt_order,
                  length, prepick, swin, all_horiz=False, delayed=True,
                  plot=False, debug=0, return_event=False, min_snr=None):
    """
    Generate multiplexed template from a Nordic file using continuous data.

    Function to read in picks from s-file then generate the template from \
    the picks within this and the wavefiles from the continuous database of \
    day-long files.  Included is a section to sanity check that the files are \
    daylong and that they start at the start of the day.  You should ensure \
    this is the case otherwise this may alter your data if your data are \
    daylong but the headers are incorrectly set.

    :type sfile: str
    :param sfile: sfilename must be the path to a seisan nordic type s-file \
            containing waveform and pick information, all other arguments can \
            be numbers save for swin which must be either P, S or all \
            (case-sensitive).
    :type contbase_list: list
    :param contbase_list: List of tuples of the form \
        ('path', 'type', 'network').  Where path is the path to the \
        continuous database, type is the directory structure, which can be \
        either Yyyyy/Rjjj.01, which is the standard IRIS Year, julian day \
        structure, or, yyyymmdd which is a single directory for every day.
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
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type all_horiz: bool
    :param all_horiz: To use both horizontal channels even if there is only \
        a pick on one of them.  Defaults to False.
    :type delayed: bool
    :param delayed: If True, each channel will begin relative to it's own \
        pick-time, if set to False, each channel will begin at the same time.
    :type plot: bool
    :param plot: Turns template plotting on or off.
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

    :returns: Newly cut template.
    :rtype: :class:`obspy.core.stream.Stream`
    """
    # Read in pick info
    event = sfile_util.readpicks(sfile)
    day = event.origins[0].time
    picks = event.picks
    pick_chans = []
    used_picks = []
    wavefiles = []
    for pick in picks:
        if not pick.waveform_id:
            print('Pick not associated with waveforms, will not use.')
            print(pick)
            continue
        station = pick.waveform_id.station_code
        channel = pick.waveform_id.channel_code
        phase = pick.phase_hint
        if station + channel not in pick_chans and phase in ['P', 'S']:
            pick_chans.append(station + channel)
            used_picks.append(pick)
            for contbase in contbase_list:
                if contbase[1] == 'yyyy/mm/dd':
                    daydir = os.path.join(str(day.year),
                                          str(day.month).zfill(2),
                                          str(day.day).zfill(2))
                elif contbase[1] == 'Yyyyy/Rjjj.01':
                    daydir = os.path.join('Y' + str(day.year),
                                          'R' + str(day.julday).zfill(3) +
                                          '.01')
                elif contbase[1] == 'yyyymmdd':
                    daydir = day.datetime.strftime('%Y%m%d')
                wavefiles += glob.glob(os.path.join(contbase[0], daydir,
                                                    '*' + station +
                                                    '.*'))
                wavefiles += glob.glob(os.path.join(contbase[0], daydir,
                                                    station + '.*'))
    picks = used_picks
    wavefiles = sorted(list(set(wavefiles)))

    # Read in waveform file
    st = Stream()
    for wavefile in wavefiles:
        st += read(wavefile)
    # Process waveform data
    st.merge(fill_value='interpolate')
    st = pre_processing.dayproc(
        st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
        samp_rate=samp_rate, starttime=day, debug=debug)
    # Cut and extract the templates
    template = template_gen(
        picks, st, length, swin, prepick=prepick, all_horiz=all_horiz,
        plot=plot, debug=debug, delayed=delayed, min_snr=min_snr)
    if return_event:
        return template, event, len(st[0].data) / samp_rate
    return template


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
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
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
    15
    """
    stations = []
    channels = []
    st_stachans = []
    templates = []
    process_lengths = []
    # Process waveform data
    st.merge(fill_value='interpolate')
    # Work out if the data are daylong or not...
    data_len = max([len(tr.data) / tr.stats.sampling_rate for tr in st])
    if 80000 < data_len < 90000:
        daylong = True
    else:
        daylong = False
    if daylong:
        starttime = min([tr.stats.starttime for tr in st])
        # Cope with the common starttime less than 1s before the start of day.
        if (starttime + 10).date > starttime.date:
            starttime = (starttime + 10).date
        else:
            starttime = starttime.date
        st = pre_processing.dayproc(
            st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
            samp_rate=samp_rate, debug=debug, parallel=parallel,
            starttime=UTCDateTime(starttime))
    else:
        st = pre_processing.shortproc(
            st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
            parallel=parallel, samp_rate=samp_rate, debug=debug)
    data_start = min([tr.stats.starttime for tr in st])
    data_end = max([tr.stats.endtime for tr in st])
    # Read QuakeML file into Catalog class
    catalog = read_events(meta_file)
    for event in catalog:
        if len(event.picks) == 0:
            warnings.warn('No picks for event %s' % event.resource_id)
            continue
        use_event = True
        # Check that the event is within the data
        for pick in event.picks:
            if not data_start < pick.time < data_end:
                debug_print("Pick outside of data span:\nPick time %s\n"
                            "Start time %s\nEnd time: %s" %
                            (str(pick.time), str(data_start), str(data_end)),
                            0, debug)
                use_event = False
        if not use_event:
            warnings.warn('Event is not within data time-span')
            continue
        # Read in pick info
        debug_print("I have found the following picks", 0, debug)
        for pick in event.picks:
            if not pick.waveform_id:
                print('Pick not associated with waveforms, will not use.')
                print(pick)
                continue
            debug_print(pick, 0, debug)
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code)
        # Check to see if all picks have a corresponding waveform
        for tr in st:
            st_stachans.append('.'.join([tr.stats.station, tr.stats.channel]))
        for i in range(len(stations)):
            if not '.'.join([stations[i], channels[i]]) in st_stachans:
                warnings.warn('No data provided for ' + stations[i] + '.' +
                              channels[i])
        st1 = st.copy()
        # Cut and extract the templates
        template = template_gen(
            event.picks, st1, length, swin, prepick=prepick, plot=plot,
            debug=debug, all_horiz=all_horiz, delayed=delayed, min_snr=min_snr)
        templates.append(template)
        process_lengths.append(len(st1[0].data) / samp_rate)
    if return_event:
        return templates, catalog, process_lengths
    return templates


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
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
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
    from obspy.clients.seishub import Client
    client = Client(url, timeout=10)
    temp_list = []
    process_lengths = []
    sub_catalogs = _group_events(
        catalog=catalog, process_len=process_len, data_pad=data_pad)
    for sub_catalog in sub_catalogs:
        # Figure out which picks we have
        all_waveform_info = []
        for event in sub_catalog:
            for pick in event.picks:
                if not pick.waveform_id:
                    print('Pick not associated with waveforms, will not use.')
                    print(pick)
                    continue
                all_waveform_info.append(pick.waveform_id)
        _all_waveform_info = []
        for w in all_waveform_info:
            _all_waveform_info.append((w.network_code,
                                       w.station_code,
                                       w.channel_code,
                                       w.location_code))
        all_waveform_info = list(set(_all_waveform_info))
        del _all_waveform_info
        all_waveform_info.sort()
        print("Fetching the following traces from SeisHub")
        for waveform_info in all_waveform_info:
            net = waveform_info[0]
            sta = waveform_info[1]
            chan = waveform_info[2]
            loc = waveform_info[3]
            if not loc:
                loc = ''
            starttime = UTCDateTime(sub_catalog[0].origins[0].time -
                                    data_pad)
            endtime = starttime + process_len
            if not endtime > sub_catalog[-1].origins[0].time + data_pad:
                raise IOError('Events do not fit in processing window')
            debug_print('start-time: %s\nend-time: %s\npick-time: %s' %
                        (str(starttime), str(endtime), str(pick.time)),
                        0, debug)
            debug_print('.'.join([net, sta, loc, chan]), 0, debug)
            if sta in client.waveform.get_station_ids(network=net):
                if 'st' not in locals():
                    st = client.waveform.get_waveform(net, sta, loc, chan,
                                                      starttime, endtime)
                else:
                    st += client.waveform.get_waveform(net, sta, loc, chan,
                                                       starttime, endtime)
            else:
                print('Station not found in SeisHub DB')
        if len(st) == 0:
            raise IOError('No waveforms found')
        if debug > 0:
            st.plot()
        print('Pre-processing data for event: %s' % event.resource_id)
        st.merge(fill_value='interpolate')
        # clients download chunks, we need to assert that the data are
        # the desired length
        for tr in st:
            tr.trim(starttime, endtime)
            print(len(tr))
        st1 = pre_processing.shortproc(
            st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
            samp_rate=samp_rate, debug=debug, parallel=True)
        for event in sub_catalog:
            template = template_gen(
                picks=event.picks, st=st1, length=length, swin=swin,
                prepick=prepick, all_horiz=all_horiz, plot=plot, debug=debug,
                delayed=delayed, min_snr=min_snr)
            process_lengths.append(len(st1[0].data) / samp_rate)
            temp_list.append(template)
            del st, st1
    if return_event:
        return temp_list, catalog, process_lengths
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
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
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
    Pre-processing data
    >>> templates[0].plot(equal_scale=False, size=(800,600)) # doctest: +SKIP

    .. figure:: ../../plots/template_gen.from_client.png
    """
    from obspy.clients.fdsn import Client
    from obspy.clients.fdsn.header import FDSNException

    client = Client(client_id)
    temp_list = []
    process_lengths = []
    # Group catalog into days and only download the data once per day
    sub_catalogs = _group_events(
        catalog=catalog, process_len=process_len, data_pad=data_pad)
    for sub_catalog in sub_catalogs:
        st = Stream()
        all_waveform_info = []
        for event in sub_catalog:
            for pick in event.picks:
                if not pick.waveform_id:
                    print('Pick not associated with waveforms, will not use.')
                    print(pick)
                    continue
                all_waveform_info.append(pick.waveform_id)
        all_waveform_info = list(set([(w.network_code, w.station_code,
                                      w.channel_code, w.location_code)
                                     for w in all_waveform_info]))
        all_waveform_info.sort()
        dropped_pick_stations = 0
        for waveform_info in all_waveform_info:
            net = waveform_info[0]
            sta = waveform_info[1]
            chan = waveform_info[2]
            loc = waveform_info[3]
            starttime = UTCDateTime(sub_catalog[0].origins[0].time -
                                    data_pad)
            endtime = starttime + process_len
            # Check that endtime is after the last event
            if not endtime > sub_catalog[-1].origins[0].time + data_pad:
                raise TemplateGenError('Events do not fit in '
                                       'processing window')
            debug_print('start-time: %s\nend-time: %s\npick-time: %s\n'
                        'pick phase: %s' %
                        (str(starttime), str(endtime), str(pick.time),
                         pick.phase_hint), 0, debug)
            debug_print('.'.join([net, sta, loc, chan]), 0, debug)
            try:
                st += client.get_waveforms(net, sta, loc, chan,
                                           starttime, endtime)
            except FDSNException:
                warnings.warn('Found no data for this station')
                dropped_pick_stations += 1
        if debug > 0:
            st.plot()
        if not st and dropped_pick_stations == len(event.picks):
            raise FDSNException('No data available, is the server down?')
        print('Pre-processing data')
        st.merge(fill_value='interpolate')
        # clients download chunks, we need to assert that the data are
        # the desired length
        for tr in st:
            tr.trim(starttime, endtime)
            if len(tr.data) == (process_len * tr.stats.sampling_rate) + 1:
                tr.data = tr.data[1:len(tr.data)]
        st1 = pre_processing.shortproc(
            st=st, lowcut=lowcut, highcut=highcut, filt_order=filt_order,
            samp_rate=samp_rate, debug=debug, parallel=True)
        if debug > 0:
            st1.plot()
        for event in sub_catalog:
            template = template_gen(
                picks=event.picks, st=st1, length=length, swin=swin,
                prepick=prepick, plot=plot, debug=debug, all_horiz=all_horiz,
                delayed=delayed, min_snr=min_snr)
            process_lengths.append(len(st1[0].data) / samp_rate)
            temp_list.append(template)
        del st, st1
    if return_event:
        return temp_list, catalog, process_lengths
    return temp_list


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
    :param swin: P, S or all, defaults to all
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
    templates = []
    process_lengths = []
    working_catalog = catalog.copy()
    # copy this here so we don't remove picks from the real catalog
    stachans = [(tr.stats.station, tr.stats.channel) for tr in st]
    for event in working_catalog:
        picks = event.picks
        for pick in picks:
            if not pick.waveform_id:
                print('Pick not associated with waveforms, will not use.')
                print(pick)
                picks.remove(pick)
                continue
            if st[0].stats.starttime < pick.time < st[0].stats.endtime:
                pick_stachan = (pick.waveform_id.station_code,
                                pick.waveform_id.channel_code)
                if pick_stachan in stachans:
                    continue
                else:
                    # Only keep a pick if there as data for it
                    picks.remove(pick)
            else:
                picks.remove(pick)
        if len(picks) > 0:
            st_clip = st.copy()
            template = template_gen(picks=picks, st=st_clip, length=length,
                                    swin=swin, prepick=prepick, plot=plot,
                                    debug=debug, all_horiz=all_horiz,
                                    delayed=delayed, min_snr=min_snr)
            process_lengths.append(st[0].stats.endtime - st[0].stats.starttime)
            templates.append(template)
    if return_event:
        return templates, catalog, process_lengths
    return templates


def template_gen(picks, st, length, swin='all', prepick=0.05,
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
    :param swin: P, S or all, defaults to all
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

    .. warning:: If there is no phase_hint included in picks, and swin=all, \
        all channels with picks will be used.
    """
    from eqcorrscan.utils.debug_log import debug_print
    from eqcorrscan.utils.plotting import pretty_template_plot as tplot
    from eqcorrscan.core.bright_lights import _rms
    stations = []
    channels = []
    st_stachans = []
    picks_copy = copy.deepcopy(picks)  # Work on a copy of the picks and leave
    # the users picks intact.
    if swin not in ['P', 'all', 'S']:
        raise IOError('Phase type is not in [all, P, S]')
    for pick in picks_copy:
        if not pick.waveform_id:
            print('Pick not associated with waveform, will not use it.')
            print(pick)
            picks_copy.remove(pick)
            continue
        if not pick.waveform_id.station_code:
            print('Pick not associated with a station, will not use it.')
            print(pick)
            picks_copy.remove(pick)
            continue
        if not pick.waveform_id.channel_code:
            print('Pick not associated with a station, will not use it.')
            print(pick)
            picks_copy.remove(pick)
            continue
        # Check to see that we are only taking the appropriate picks
        if swin == 'all':
            # Annoying comparability with seisan two channel codes
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code[0] +
                            pick.waveform_id.channel_code[-1])
        elif swin == 'P' and 'P' in pick.phase_hint.upper():
            # Use the 'in' statement to cope with phase names like 'PN' etc.
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code[0] +
                            pick.waveform_id.channel_code[-1])
        elif swin == 'S' and 'S' in pick.phase_hint.upper():
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code[0] +
                            pick.waveform_id.channel_code[-1])
    for tr in st:
        # Check that the data can be represented by float16, and check they
        # are not all zeros
        if np.all(tr.data.astype(np.float16) == 0):
            warnings.warn('Trace is all zeros at float16 level,'
                          'either gain or check. Not using in template.')
            print(tr)
            st.remove(tr)
        else:
            st_stachans.append('.'.join([tr.stats.station, tr.stats.channel]))
    for i, station in enumerate(stations):
        if '.'.join([station, channels[i]]) not in st_stachans and debug > 0:
            warnings.warn('No data provided for ' + station + '.' +
                          channels[i])
    if plot:
        stplot = st.copy()
    # Get the earliest pick-time and use that if we are not using delayed.
    event_start_time = min([pick.time for pick in picks_copy])
    event_start_time -= prepick
    # Cut the data
    st1 = Stream()
    for tr in st:
        noise_amp = _rms(tr.data)
        used_tr = False
        for pick in picks_copy:
            starttime = None
            if swin == 'all':
                if not pick.phase_hint:
                    msg = 'Pick for ' + pick.waveform_id.station_code + '.' +\
                        pick.waveform_id.channel_code + ' has no phase ' +\
                        'hint given, you should not use this template for ' +\
                        'cross-correlation re-picking!'
                    warnings.warn(msg)
                    if pick.waveform_id.station_code == tr.stats.station and \
                            pick.waveform_id.channel_code == \
                            tr.stats.channel:
                        starttime = pick.time - prepick
                else:
                    if pick.waveform_id.station_code == tr.stats.station and \
                            pick.waveform_id.channel_code ==\
                            tr.stats.channel:
                        starttime = pick.time - prepick
                    # Cope with taking all the horizontals for S-picks.
                    elif all_horiz and pick.waveform_id.station_code ==\
                            tr.stats.station:
                        if tr.stats.channel[-1] not in ['Z', 'U']\
                                and pick.phase_hint == 'S':
                            starttime = pick.time - prepick
            else:
                if pick.waveform_id.station_code == tr.stats.station and\
                        swin in pick.phase_hint.upper():
                    starttime = pick.time - prepick
            if starttime is not None:
                debug_print("Cutting " + tr.stats.station + '.' +
                            tr.stats.channel, 0, debug)
                if not delayed:
                    starttime = event_start_time
                tr_cut = tr.copy().trim(starttime=starttime,
                                        endtime=starttime + length,
                                        nearest_sample=False)
                if len(tr_cut.data) == 0:
                    print('No data provided for %s.%s starting at %s' %
                          (tr.stats.station, tr.stats.channel, str(starttime)))
                    continue
                # Ensure that the template is the correct length
                if len(tr_cut.data) == (tr_cut.stats.sampling_rate *
                                        length) + 1:
                    tr_cut.data = tr_cut.data[0:-1]
                debug_print('Cut starttime = %s\nCut endtime %s' %
                            (str(tr_cut.stats.starttime),
                             str(tr_cut.stats.endtime)), 0, debug)
                if min_snr is not None and \
                   max(tr_cut.data) / noise_amp < min_snr:
                    print('Signal-to-noise ratio below threshold for %s.%s' %
                          (tr_cut.stats.station, tr_cut.stats.channel))
                    continue
                st1 += tr_cut
                used_tr = True
        if not used_tr:
            debug_print('No pick for ' + tr.stats.station + '.' +
                        tr.stats.channel, 0, debug)
    if plot:
        background = stplot.trim(
            st1.sort(['starttime'])[0].stats.starttime - 10,
            st1.sort(['starttime'])[-1].stats.endtime + 10)
        tplot(st1, background=background, picks=picks_copy,
              title='Template for ' + str(st1[0].stats.starttime))
        del stplot
    del st
    return st1


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
        new_template = pre_processing.shortproc(st=new_template,
                                                lowcut=lowcut,
                                                highcut=highcut,
                                                filt_order=filt_order,
                                                samp_rate=samp_rate,
                                                debug=0)
    # Loop through the stack and trim!
    out = Stream()
    for tr in new_template:
        # Find the matching delay
        delay = [d[2] for d in delays if d[0] == tr.stats.station and
                 d[1] == tr.stats.channel[-1]]
        if Z_include and len(delay) == 0:
            delay = [d[2] for d in delays if d[0] == tr.stats.station]
        if len(delay) == 0:
            msg = ' '.join(['No matching template channel found for stack',
                            'channel', tr.stats.station, tr.stats.channel])
            warnings.warn(msg)
            new_template.remove(tr)
        else:
            for d in delay:
                out += tr.copy().trim(
                    starttime=tr.stats.starttime + d + pre_pad - pre_pick,
                    endtime=tr.stats.starttime + d + pre_pad + length -
                    pre_pick)
    return out


def _group_events(catalog, process_len, data_pad):
    """
    Internal function to group events into sub-catalogs based on process_len.

    :param catalog: Catalog to groups into sub-catalogs
    :type catalog: obspy.core.event.Catalog
    :param process_len: Length in seconds that data will be processed in
    :type process_len: int
    :param data_pad: Length of data (in seconds) required before and after \
        any event for processing, use to reduce edge-effects of filtering on \
        the templates.
    :type data_pad: int

    :return: List of catalogs
    :rtype: list
    """
    from obspy.core.event import Catalog
    # case for catalog only containing one event
    if len(catalog) == 1:
        return [catalog]
    sub_catalogs = []
    # Sort catalog by date
    cat_list = [(event, event.origins[0].time) for event in catalog]
    cat_list.sort(key=lambda tup: tup[1])
    catalog = Catalog([tup[0] for tup in cat_list])
    sub_catalog = Catalog([catalog[0]])
    for event in catalog[1:]:
        if (event.origins[0].time + data_pad) - \
                (sub_catalog[0].origins[0].time - data_pad) < process_len:
            sub_catalog.append(event)
        else:
            sub_catalogs.append(sub_catalog)
            sub_catalog = Catalog([event])
    sub_catalogs.append(sub_catalog)
    return sub_catalogs

if __name__ == "__main__":
    import doctest
    doctest.testmod()
