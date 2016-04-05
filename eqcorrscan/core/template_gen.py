#!/usr/bin/python
r"""Functions to generate template waveforms and information to go with them \
for the application of cross-correlation of seismic data for the detection of \
repeating events.

.. note:: All of these functions work for a single template, however all of \
    them call _template_gen, which takes care of pick association and \
    cutting.  If you have many templates in one day of data it would be \
    simple to write a wrapper that cuts multiple templates from one day \
    of processed data rather than re-processing the same day of data \
    for each template.

.. note:: All functions use obspy filters, which are implemented such that \
    if both highcut and lowcut are set a bandpass filter will be used, \
    but of highcut is not set (None) then a highpass filter will be used and \
    if only the highcut is set then a lowpass filter will be used.

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


def from_sac(sac_files, lowcut, highcut, samp_rate, filt_order, length, swin,
             prepick=0.05, debug=0, plot=False):
    """Function to read picks and waveforms from SAC data, and generate a \
    template from these. Usually sac_files is a list of all single-channel \
    SAC files for a given event, a single, multi-channel template will be \
    created from these traces.

    :type sac_files: list or stream
    :param sac_files: List or stream of sac waveforms, or list of paths to \
        sac waveforms.
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will look in template \
            defaults file
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will look in template \
            defaults file
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz, if set to None will look in \
            template defaults file
    :type filt_order: int
    :param filt_order: Filter level, if set to None will look in \
            template defaults file
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type length: float
    :param length: Extract length in seconds, if None will look in template \
            defaults file.
    :type prepick: float
    :param prepick: Length to extract prior to the pick in seconds.
    :type debug: int
    :param debug: Debug level, higher number=more output.
    :type plot: bool
    :param plot: Turns template plotting on or off.

    :returns: obspy.Stream Newly cut template

    .. note:: This functionality is not supported for obspy versions below \
        1.0.0 as references times are not read in by SACIO, which are needed \
        for defining pick times.
    """
    from obspy import read, Stream
    from eqcorrscan.utils.sac_util import sactoevent
    from eqcorrscan.utils import pre_processing
    # Check whether sac_files is a stream or a list
    if isinstance(sac_files, list):
        if isinstance(sac_files[0], str) or isinstance(sac_files[0], unicode):
            sac_files = [read(sac_file)[0] for sac_file in sac_files]
        if isinstance(sac_files[0], Stream):
            # This is a list of streams...
            st = sac_files[0]
            for sac_file in sac_files[1:]:
                st += sac_file
        st = Stream(sac_files)
    elif isinstance(sac_files, Stream):
        st = sac_files
    # Make an event object...
    event = sactoevent(st, debug=debug)
    # Process the data
    st.merge(fill_value='interpolate')
    st = pre_processing.shortproc(st, lowcut, highcut, filt_order,
                                  samp_rate, debug)
    template = _template_gen(picks=event.picks, st=st, length=length,
                             swin=swin, prepick=prepick, plot=plot,
                             debug=debug)
    return template


def from_sfile(sfile, lowcut, highcut, samp_rate, filt_order, length, swin,
               prepick=0.05, debug=0, plot=False):
    r"""Function to read in picks from sfile then generate the template from \
    the picks within this and the wavefile found in the pick file.

    :type sfile: string
    :param sfile: sfilename must be the \
        path to a seisan nordic type s-file containing waveform and pick \
        information.
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will look in template \
            defaults file
    :type highcut: float
    :param highcut: High cut (Hz), if set to None will look in template \
            defaults file
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz, if set to None will look in \
            template defaults file
    :type filt_order: int
    :param filt_order: Filter level, if set to None will look in \
            template defaults file
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type length: float
    :param length: Extract length in seconds, if None will look in template \
            defaults file.
    :type prepick: float
    :param prepick: Length to extract prior to the pick in seconds.
    :type debug: int
    :param debug: Debug level, higher number=more output.
    :type plot: bool
    :param plot: Turns template plotting on or off.

    :returns: obspy.Stream Newly cut template

    .. warning:: This will use whatever data is pointed to in the s-file, if \
        this is not the coninuous data, we recommend using other functions. \
        Differences in processing between short files and day-long files \
        (inherent to resampling) will produce lower cross-correlations.
    """
    # Perform some checks first
    import os
    if not os.path.isfile(sfile):
        raise IOError('sfile does not exist')

    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils import sfile_util
    from obspy import read as obsread
    # Read in the header of the sfile
    wavefiles = sfile_util.readwavename(sfile)
    pathparts = sfile.split('/')[0:-1]
    new_path_parts = []
    for part in pathparts:
        if part == 'REA':
            part = 'WAV'
        new_path_parts.append(part)
    # * argument to allow .join() to accept a list
    wavpath = os.path.join(*new_path_parts) + '/'
    # In case of absolute paths (not handled with .split() --> .join())
    if sfile[0] == '/':
        wavpath = '/' + wavpath
    # Read in waveform file
    for wavefile in wavefiles:
        print(''.join(["I am going to read waveform data from: ", wavpath,
                       wavefile]))
        if 'st' not in locals():
            st = obsread(wavpath + wavefile)
        else:
            st += obsread(wavpath + wavefile)
    for tr in st:
        if tr.stats.sampling_rate < samp_rate:
            print('Sampling rate of data is lower than sampling rate asked ' +
                  'for')
            print('Not good practice for correlations: I will not do this')
            raise ValueError("Trace: " + tr.stats.station +
                             " sampling rate: " + str(tr.stats.sampling_rate))
    # Read in pick info
    catalog = sfile_util.readpicks(sfile)
    # Read the list of Picks for this event
    picks = catalog[0].picks
    print("I have found the following picks")
    for pick in picks:
        print(' '.join([pick.waveform_id.station_code,
                        pick.waveform_id.channel_code, pick.phase_hint,
                        str(pick.time)]))

    # Process waveform data
    st.merge(fill_value='interpolate')
    st = pre_processing.shortproc(st, lowcut, highcut, filt_order,
                                  samp_rate, debug)
    st1 = _template_gen(picks=picks, st=st, length=length, swin=swin,
                        prepick=prepick, plot=plot, debug=debug)
    return st1


def from_contbase(sfile, contbase_list, lowcut, highcut, samp_rate, filt_order,
                  length, prepick, swin, debug=0, plot=False):
    r"""Function to read in picks from sfile then generate the template from \
    the picks within this and the wavefiles from the continous database of \
    day-long files.  Included is a section to sanity check that the files are \
    daylong and that they start at the start of the day.  You should ensure \
    this is the case otherwise this may alter your data if your data are \
    daylong but the headers are incorrectly set.

    :type sfile: string
    :param sfile: sfilename must be the path to a seisan nordic type s-file \
            containing waveform and pick information, all other arguments can \
            be numbers save for swin which must be either P, S or all \
            (case-sensitive).
    :type contbase_list: List of tuple of string
    :param contbase_list: List of tuples of the form \
        ['path', 'type', 'network'].  Where path is the path to the \
        continuous database, type is the directory structure, which can be \
        either Yyyyy/Rjjj.01, which is the standard IRIS Year, julian day \
        structure, or, yyyymmdd which is a single directory for every day.
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will look in template \
            defaults file
    :type highcut: float
    :param lowcut: High cut (Hz), if set to None will look in template \
            defaults file
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz, if set to None will look in \
            template defaults file
    :type filt_order: int
    :param filt_order: Filter level, if set to None will look in \
            template defaults file
    :type length: float
    :param length: Extract length in seconds, if None will look in template \
            defaults file.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type plot: bool
    :param plot: Turns template plotting on or off.

    :returns: obspy.Stream Newly cut template
    """
    # Perform some checks first
    import os
    if not os.path.isfile(sfile):
        raise IOError('sfile does not exist')

    # import some things
    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils import sfile_util
    import glob
    from obspy import read as obsread

    # Read in the header of the sfile
    event = sfile_util.readheader(sfile)
    day = event.origins[0].time

    # Read in pick info
    catalog = sfile_util.readpicks(sfile)
    picks = catalog[0].picks
    print("I have found the following picks")
    pick_chans = []
    used_picks = []
    for pick in picks:
        station = pick.waveform_id.station_code
        channel = pick.waveform_id.channel_code
        phase = pick.phase_hint
        pcktime = pick.time
        if station + channel not in pick_chans and phase in ['P', 'S']:
            pick_chans.append(station + channel)
            used_picks.append(pick)
            print(pick)
            # #########Left off here
            for contbase in contbase_list:
                if contbase[1] == 'yyyy/mm/dd':
                    daydir = os.path.join([str(day.year),
                                           str(day.month).zfill(2),
                                           str(day.day).zfill(2)])
                elif contbase[1] == 'Yyyyy/Rjjj.01':
                    daydir = os.path.join(['Y' + str(day.year),
                                           'R' + str(day.julday).zfill(3) +
                                           '.01'])
                elif contbase[1] == 'yyyymmdd':
                    daydir = day.datetime.strftime('%Y%m%d')
                if 'wavefiles' not in locals():
                    wavefiles = (glob.glob(os.path.join([contbase[0], daydir,
                                                         '*' + station +
                                                         '.*'])))
                else:
                    wavefiles += glob.glob(os.path.join([contbase[0], daydir,
                                                         '*' + station +
                                                         '.*']))
        elif phase in ['P', 'S']:
            print(' '.join(['Duplicate pick', station, channel,
                            phase, str(pcktime)]))
        elif phase == 'IAML':
            print(' '.join(['Amplitude pick', station, channel,
                            phase, str(pcktime)]))
    picks = used_picks
    wavefiles = list(set(wavefiles))

    # Read in waveform file
    wavefiles.sort()
    for wavefile in wavefiles:
        print("I am going to read waveform data from: " + wavefile)
        if 'st' not in locals():
            st = obsread(wavefile)
        else:
            st += obsread(wavefile)
    # Process waveform data
    st.merge(fill_value='interpolate')
    for tr in st:
        tr = pre_processing.dayproc(tr, lowcut, highcut, filt_order,
                                    samp_rate, debug, day)
    # Cut and extract the templates
    st1 = _template_gen(picks, st, length, swin, prepick=prepick, plot=plot,
                        debug=debug)
    return st1


def from_quakeml(quakeml, st, lowcut, highcut, samp_rate, filt_order,
                 length, prepick, swin, debug=0, plot=False):
    r"""Function to generate a template from a local quakeml file \
    and an obspy.Stream object.

    :type quakeml: string
    :param quakeml: QuakeML file containing pick information, can contain \
        multiple events.
    :type st: class: obspy.Stream
    :param st: Stream containing waveform data for template (hopefully). \
        Note that this should be the same length of stream as you will use \
        for the continuous detection, e.g. if you detect in day-long files, \
        give this a day-long file!
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will look in template \
            defaults file
    :type highcut: float
    :param lowcut: High cut (Hz), if set to None will look in template \
            defaults file
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz, if set to None will look in \
            template defaults file
    :type filt_order: int
    :param filt_order: Filter level, if set to None will look in \
            template defaults file
    :type length: float
    :param length: Extract length in seconds, if None will look in template \
            defaults file.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type plot: bool
    :param plot: Display template plots or not

    :returns: list of obspy.Stream Newly cut templates

    .. warning:: We suggest giving this function a full day of data, to \
        ensure templates are generated with **exactly** the same processing \
        as the continuous data.  Not doing this will result in slightly \
        reduced cross-correlation values.
    """
    # Perform some checks first
    import os
    import warnings
    if not os.path.isfile(quakeml):
        raise IOError('QuakeML file does not exist')
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy import read_events
    else:
        from obspy import readEvents as read_events
    from obspy import UTCDateTime
    from eqcorrscan.utils import pre_processing
    stations = []
    channels = []
    st_stachans = []
    # Process waveform data
    st.merge(fill_value='interpolate')
    for tr in st:
        tr = pre_processing.dayproc(tr, lowcut, highcut, filt_order,
                                    samp_rate, debug=debug,
                                    starttime=UTCDateTime(tr.stats.
                                                          starttime.date))
    # Read QuakeML file into Catalog class
    catalog = read_events(quakeml)
    templates = []
    for event in catalog:
        # Read in pick info
        print("I have found the following picks")
        for pick in event.picks:
            print(' '.join([pick.waveform_id.station_code,
                            pick.waveform_id.channel_code,
                            pick.phase_hint, str(pick.time)]))
            stations.append(pick.waveform_id.station_code)
            channels.append(pick.waveform_id.channel_code)
        # Check to see if all picks have a corresponding waveform
        for tr in st:
            st_stachans.append('.'.join([tr.stats.station, tr.stats.channel]))
        for i in xrange(len(stations)):
            if not '.'.join([stations[i], channels[i]]) in st_stachans:
                warnings.warn('No data provided for ' + stations[i] + '.' +
                              channels[i])
        st1 = st.copy()
        # Cut and extract the templates
        template = _template_gen(event.picks, st1, length, swin,
                                 prepick=prepick, plot=plot, debug=debug)
        templates.append(template)
    return templates


def from_seishub(catalog, url, lowcut, highcut, samp_rate, filt_order,
                 length, prepick, swin, debug=0, plot=False):
    r"""Function to generate templates from a SeisHub database.Must be given \
    an obspy.Catalog class and the SeisHub url as input. The function returns \
    a list of obspy.Stream classes containting steams for each desired \
    template.

    :type catalog: obspy.Catalog
    :param catalog: Catalog class containing desired template events
    :type url: string
    :param url: url of SeisHub database instance
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will look in template \
            defaults file
    :type highcut: float
    :param lowcut: High cut (Hz), if set to None will look in template \
            defaults file
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz, if set to None will look in \
            template defaults file
    :type filt_order: int
    :param filt_order: Filter level, if set to None will look in \
            template defaults file
    :type length: float
    :param length: Extract length in seconds, if None will look in template \
            defaults file.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type plot: bool
    :param plot: Plot templates or not.

    :returns: obspy.Stream Newly cut template
    """
    # This import section copes with namespace changes between obspy versions
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.seishub import Client
    else:
        from obspy.seishub import Client
    from eqcorrscan.utils import pre_processing
    from obspy import UTCDateTime
    client = Client(url)
    temp_list = []
    for event in catalog:
        # Figure out which picks we have
        day = event.origins[0].time
        picks = event.picks
        print("Fetching the following traces from SeisHub")
        for pick in picks:
            net = pick.waveform_id.network_code
            sta = pick.waveform_id.station_code
            chan = pick.waveform_id.channel_code
            loc = pick.waveform_id.location_code
            starttime = UTCDateTime(pick.time.date)
            endtime = starttime + 86400
            # Here we download a full day of data.  We do this so that minor
            # differences in processing during processing due to the effect
            # of resampling do not impinge on our cross-correaltions.

            if debug > 0:
                print('start-time: ' + str(starttime))
                print('end-time: ' + str(endtime))
                print('pick-time: ' + str(pick.time))
            print('.'.join([net, sta, loc, chan]))
            if sta in client.waveform.getStationIds(network=net):
                if 'st' not in locals():
                    st = client.waveform.getWaveform(net, sta, loc, chan,
                                                     starttime, endtime)
                else:
                    st += client.waveform.getWaveform(net, sta, loc, chan,
                                                      starttime, endtime)
            else:
                print('Station not found in SeisHub DB')
        if debug > 0:
            st.plot()
        print('Preprocessing data for event: '+str(event.resource_id))
        st.merge(fill_value='interpolate')
        st1 = pre_processing.dayproc(st, lowcut, highcut, filt_order,
                                     samp_rate, starttime=starttime,
                                     debug=debug)
        template = _template_gen(event.picks, st1, length, swin, prepick,
                                 plot=plot, debug=debug)
        del st, st1
        temp_list.append(template)
    return temp_list


def from_client(catalog, client_id, lowcut, highcut, samp_rate, filt_order,
                length, prepick, swin, debug=0, plot=False):
    r"""Function to generate templates from a SeisHub database.Must be given \
    an obspy.Catalog class and the SeisHub url as input. The function returns \
    a list of obspy.Stream classes containting steams for each desired \
    template.

    :type catalog: obspy.Catalog
    :param catalog: Catalog class containing desired template events
    :type url: string
    :param url: url of SeisHub database instance
    :type lowcut: float
    :param lowcut: Low cut (Hz), if set to None will look in template\
            defaults file
    :type highcut: float
    :param lowcut: High cut (Hz), if set to None will look in template\
            defaults file
    :type samp_rate: float
    :param samp_rate: New sampling rate in Hz, if set to None will look in\
            template defaults file
    :type filt_order: int
    :param filt_order: Filter level, if set to None will look in\
            template defaults file
    :type length: float
    :param length: Extract length in seconds, if None will look in template\
            defaults file.
    :type prepick: float
    :param prepick: Pre-pick time in seconds
    :type swin: str
    :param swin: Either 'all', 'P' or 'S', to select which phases to output.
    :type debug: int
    :param debug: Level of debugging output, higher=more
    :type plot: bool
    :param plot: Plot templates or not.

    :returns: obspy.Stream Newly cut template
    """
    # This import section copes with namespace changes between obspy versions
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
        from obspy.clients.fdsn.header import FDSNException
    else:
        from obspy.fdsn import Client
        from obspy.fdsn.header import FDSNException
    from eqcorrscan.utils import pre_processing
    from obspy import UTCDateTime
    import warnings

    client = Client(client_id)
    temp_list = []
    for event in catalog:
        # Figure out which picks we have
        day = event.origins[0].time
        print("Fetching the following traces from " + client_id)
        for pick in event.picks:
            net = pick.waveform_id.network_code
            sta = pick.waveform_id.station_code
            chan = pick.waveform_id.channel_code
            loc = pick.waveform_id.location_code
            starttime = UTCDateTime(pick.time.date)
            endtime = starttime + 86400
            # Here we download a full day of data.  We do this so that minor
            # differences in processing during processing due to the effect
            # of resampling do not impinge on our cross-correaltions.
            if debug > 0:
                print('start-time: ' + str(starttime))
                print('end-time: ' + str(endtime))
                print('pick-time: ' + str(pick.time))
                print('pick phase: ' + pick.phase_hint)
            print('.'.join([net, sta, loc, chan]))
            if 'st' not in locals():
                try:
                    st = client.get_waveforms(net, sta, loc, chan,
                                              starttime, endtime)
                except FDSNException:
                    warnings.warn('Found no data for this station')
            else:
                try:
                    st += client.get_waveforms(net, sta, loc, chan,
                                               starttime, endtime)
                except FDSNException:
                    warnings.warn('Found no data for this station')
        if debug > 0:
            st.plot()
        print('Pre-processing data for event: '+str(event.resource_id))
        st.merge(fill_value='interpolate')
        st1 = pre_processing.dayproc(st, lowcut, highcut, filt_order,
                                     samp_rate, starttime=starttime,
                                     debug=debug, parallel=True)
        if debug > 0:
            st1.plot()
        template = _template_gen(event.picks, st1, length, swin, prepick,
                                 plot=plot, debug=debug)
        del st, st1
        temp_list.append(template)
    return temp_list


def _template_gen(picks, st, length, swin='all', prepick=0.05, plot=False,
                  debug=0):
    r"""Function to generate a cut template in the obspy \
    Stream class from a given set of picks and data, also in an obspy stream \
    class.  Should be given pre-processed data (downsampled and filtered).

    :type picks: List of obspy.core.event.Pick
    :param picks: Picks to extract data around
    :type st: :class: 'obspy.Stream'
    :param st: Stream to etract templates from
    :type length: float
    :param length: Length of template in seconds
    :type swin: string
    :param swin: P, S or all, defaults to all
    :type prepick: float
    :param prepick: Length in seconds to extract before the pick time \
            default is 0.05 seconds
    :type plot: bool
    :param plot: To plot the template or not, default is True
    :type debug: int
    :param debug: Debug output level from 0-5.

    :returns: obspy.Stream Newly cut template.

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
    import copy
    from eqcorrscan.utils.plotting import pretty_template_plot as\
        tplot
    from obspy import Stream
    import warnings
    stations = []
    channels = []
    st_stachans = []
    if not swin in ['P', 'all', 'S']:
        raise IOError('Phase type is not in [all, P, S]')
    for pick in picks:
        # Check to see that we are only taking the appropriate picks
        if swin == 'all':
            # Annoying compatability with seisan two channel codes
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
        st_stachans.append('.'.join([tr.stats.station, tr.stats.channel]))
    for i, station in enumerate(stations):
        if not '.'.join([station, channels[i]]) in st_stachans:
            warnings.warn('No data provided for ' + station + '.' +
                          channels[i])
    # Select which channels we actually have picks for
    for tr in st:
        if tr.stats.station in stations:
            # This is used to cope with seisan handling channel codes as
            # two charectar codes, internally we will do the same.
            if len(tr.stats.channel) == 3:
                temp_channel = tr.stats.channel[0] + tr.stats.channel[2]
            elif len(tr.stats.channel) == 2:
                temp_channel = tr.stats.channel
            # Take all channels
            tr.stats.channel = temp_channel
            if 'st1' not in locals():
                st1 = Stream(tr)
            else:
                st1 += tr
    if 'st1' not in locals():
        msg = ('No data available for these picks or no picks match ' +
               'these data!  Will not error, but you should check yo self')
        warnings.warn(msg)
        return
    st = copy.deepcopy(st1)
    del st1
    if plot:
        stplot = st.copy()
    # Cut the data
    for tr in st:
        if 'starttime' in locals():
            del starttime
        if swin == 'all':
            for pick in picks:
                if not pick.phase_hint:
                    msg = 'Pick for ' + pick.waveform_id.station_code + '.' +\
                        pick.waveform_id.channel_code + ' has no phase ' +\
                        'hint given, you should not use this template for ' +\
                        'cross-correlation re-picking!'
                    warnings.warn(msg)
                    if pick.waveform_id.station_code == tr.stats.station and \
                            pick.waveform_id.channel_code[0] + \
                            pick.waveform_id.channel_code[-1] == \
                            tr.stats.channel:
                        starttime = pick.time - prepick
                else:
                    # If there is phase information then we should use our
                    # convention.
                    if pick.waveform_id.station_code == tr.stats.station and \
                            pick.waveform_id.channel_code[0] + \
                            pick.waveform_id.channel_code[-1] ==\
                            tr.stats.channel \
                            and 'P' in pick.phase_hint.upper():
                        starttime = pick.time - prepick
                    elif pick.waveform_id.station_code == tr.stats.station and\
                            tr.stats.channel[-1] in ['1', '2', 'N',
                                                     'E', 'R', 'T'] and\
                            'S' in pick.phase_hint.upper():
                        starttime = pick.time - prepick
        else:
            for pick in picks:
                if pick.waveform_id.station_code == tr.stats.station and\
                        swin in pick.phase_hint.upper():
                    starttime = pick.time - prepick
        if 'starttime' in locals():
            if debug > 0:
                print("Cutting " + tr.stats.station + '.' + tr.stats.channel)
            tr.trim(starttime=starttime, endtime=starttime + length,
                    nearest_sample=False)
            if debug > 0:
                print('Cut starttime = ' + str(tr.stats.starttime))
                print('Cut endtime = ' + str(tr.stats.endtime))
            if 'st1' not in locals():
                st1 = Stream(tr)
            else:
                st1 += tr
        elif debug > 0:
            print('No pick for ' + tr.stats.station + '.' + tr.stats.channel)
        # Ensure that the template is the correct length
        if len(tr.data) == (tr.stats.sampling_rate * length) + 1:
            tr.data = tr.data[0:-1]
    if plot:
        background = stplot.trim(st1.sort(['starttime'])[0].stats.starttime -
                                 10,
                                 st1.sort(['starttime'])[-1].stats.endtime +
                                 10)
        tplot(st1, background=background,
              title='Template for '+str(st1[0].stats.starttime),
              picks=picks)
        del stplot
    del st
    # st1.plot(size=(800,600))
    return st1


def extract_from_stack(stack, template, length, pre_pick, pre_pad,
                       Z_include=False, pre_processed=True, samp_rate=None,
                       lowcut=None, highcut=None, filt_order=3):
    r"""Function to extract a new template from a stack of previous detections.
    Requires the stack, the template used to make the detections for the \
    stack, and we need to know if the stack has been pre-processed.

    :type stack: :class:obspy.Stream
    :param stack: Waveform stack from detections.  Can be of any length and \
        can have delays already included, or not.
    :type template: :class:obspy.Stream
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

    :returns: obspy.Stream Newly cut template
    """
    from eqcorrscan.utils import pre_processing
    import warnings
    new_template = stack.copy()
    # Copy the data before we trim it to keep the stack safe
    # Get the earliest time in the template as this is when the detection is
    # taken.
    mintime = min([tr.stats.starttime for tr in template])
    # Generate a list of tuples of (station, channel, delay) with delay in
    # seconds
    delays = [(tr.stats.station, tr.stats.channel[-1],
               tr.stats.starttime - mintime) for tr in template]
    # Loop through the stack and trim!
    for tr in new_template:
        # Process the data if necessary
        if not pre_processed:
            new_template = pre_processing.shortproc(new_template, lowcut,
                                                    highcut, filt_order,
                                                    samp_rate, 0)
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
        elif len(delay) > 1:
            msg = ' '.join(['Multiple delays found for stack channel',
                            tr.stats.station, tr.stats.channel])
            warnings.warn(msg)
        else:
            tr.trim(starttime=tr.stats.starttime + delay[0] + pre_pad -
                    pre_pick,
                    endtime=tr.stats.starttime + delay[0] + pre_pad + length -
                    pre_pick)
    return new_template


if __name__ == "__main__":
    import doctest
    doctest.testmod()
