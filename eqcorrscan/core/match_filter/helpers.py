"""
Helper functions for network matched-filter detection of seismic data.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import contextlib
import os
import shutil
import tarfile
import tempfile
import logging

from typing import List, Set

import numpy as np
from obspy import Stream
from obspy.core.event import Event

from eqcorrscan.utils.correlate import get_array_xcorr
from eqcorrscan.core.match_filter import Template


Logger = logging.getLogger(__name__)


def _group_templates(
        templates: List[Template],
        st_seed_ids: Set[str],
        group_size: int,
        fill_groups: bool = True
) -> List[List[Template]]:
    """
    Group templates to reduce dissimilar traces

    :param templates:
        Templates to group together
    :param st_seed_ids:
        Seed ids in the stream to be matched with
    :param group_size:
        Maximum group size - will not exceed this size
    :param fill_groups:
        Whether to fill groups up to group size, or just use the
        most similar groupings below group_size.

    :return:
        List of lists of templates grouped.
    """
    # Get overlapping seed ids so that we only group based on the
    # channels we have in the data.
    template_seed_ids = {
        template.name:
            {tr.id for tr in template.st}.intersection(st_seed_ids)
        for template in templates}
    # Get initial groups with matched traces - sort by length
    sorted_templates = sorted(
        template_seed_ids, key=lambda key: len(template_seed_ids[key]),
        reverse=True)
    groups, group = [], [sorted_templates[0]]
    for i in range(1, len(sorted_templates)):
        # Check that we don't exceed the maximum group size first
        if len(group) >= group_size or template_seed_ids[sorted_templates[i]] != template_seed_ids[group[0]]:
            groups.append(group)
            group = [sorted_templates[i]]
        else:
            group.append(sorted_templates[i])
    # Get the final group
    groups.append(group)
    if len(groups) > 1:
        # Re-group to make the largest groups possible with minimal nan-traces
        out_groups, out_group = [], groups[0]
        for i in range(1, len(groups)):
            if len(out_group) + len(groups[i]) <= group_size:
                out_group.extend(groups[i])
            elif fill_groups:
                n_included = group_size - len(out_group)
                out_group.extend(groups[i][0:n_included])
                # This should now be a full group
                out_groups.append(out_group)
                out_group = groups[i][n_included:]
            else:
                out_groups.append(out_group)
                out_group = groups[i]
        # Get final group
        out_groups.append(groups[i])
    else:
        out_groups = groups
    # Convert from groups if template names to groups of templates
    template_dict = {t.name: t for t in templates}
    out_groups = [[template_dict[t] for t in out_group] for out_group in out_groups]
    return out_groups


@contextlib.contextmanager
def temporary_directory():
    """ make a temporary directory, yeild its name, cleanup on exit """
    dir_name = tempfile.mkdtemp()
    yield dir_name
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)


def get_waveform_client(waveform_client):
    """
    Bind a `get_waveforms_bulk` method to client if it doesn't already have one

    :param waveform_client: Obspy client with a `get_waveforms` method

    :returns: waveform_client with `get_waveforms_bulk`.
    """
    def _get_waveforms_bulk_naive(self, bulk_arg):
        """ naive implementation of get_waveforms_bulk that uses iteration. """
        st = Stream()
        for arg in bulk_arg:
            st += self.get_waveforms(*arg)
        return st

    # add waveform_bulk method dynamically if it doesn't exist already
    if not hasattr(waveform_client, "get_waveforms_bulk"):
        bound_method = _get_waveforms_bulk_naive.__get__(waveform_client)
        setattr(waveform_client, "get_waveforms_bulk", bound_method)

    return waveform_client


def _spike_test(stream, percent=0.99, multiplier=1e7):
    """
    Check for very large spikes in data and raise an error if found.

    :param stream: Stream to look for spikes in.
    :type stream: :class:`obspy.core.stream.Stream`
    :param percent: Percentage as a decimal to calculate range for.
    :type percent: float
    :param multiplier: Multiplier of range to define a spike.
    :type multiplier: float
    """
    from eqcorrscan.core.match_filter.matched_filter import MatchFilterError

    list_ids = []
    for tr in stream:
        if (tr.data > 2 * np.max(np.sort(
                np.abs(tr.data))[0:int(percent * len(tr.data))]
                                 ) * multiplier).sum() > 0:
            list_ids.append(tr.id)
    if list_ids != []:
        ids = ', '.join(list_ids)
        msg = ('Spikes above ' + str(multiplier) +
               ' of the range of ' + str(percent) +
               ' of the data present, check:\n' + ids + '.\n'
               'This would otherwise likely result in an issue during ' +
               'FFT prior to cross-correlation.\n' +
               'If you think this spike is real please report ' +
               'this as a bug.')
        raise MatchFilterError(msg)


def _total_microsec(t1, t2):
    """
    Calculate difference between two datetime stamps in microseconds.

    :type t1: :class: `datetime.datetime`
    :type t2: :class: `datetime.datetime`
    :return: int

    .. rubric:: Example

    >>> from obspy import UTCDateTime
    >>> print(_total_microsec(UTCDateTime(2013, 1, 1).datetime,
    ...                       UTCDateTime(2014, 1, 1).datetime))
    -31536000000000
    """
    td = t1 - t2
    return (td.seconds + td.days * 24 * 3600) * 10 ** 6 + td.microseconds


def _templates_match(t, family_file):
    """
    Return True if a tribe matches a family file path.

    :type t: Tribe
    :type family_file: str
    :return: bool
    """
    return t.name == family_file.split(os.sep)[-1].split('_detections.csv')[0]


def _test_event_similarity(event_1, event_2, verbose=False, shallow=False):
    """
    Check the similarity of the components of obspy events, discounting
    resource IDs, which are not maintained in nordic files.

    :type event_1: obspy.core.event.Event
    :param event_1: First event
    :type event_2: obspy.core.event.Event
    :param event_2: Comparison event
    :type verbose: bool
    :param verbose: If true and fails will output why it fails.

    :return: bool
    """
    if not isinstance(event_1, Event) or not isinstance(event_2, Event):
        raise NotImplementedError('Cannot compare things that are not Events')
    # Check origins
    if len(event_1.origins) != len(event_2.origins):
        return False
    for ori_1, ori_2 in zip(event_1.origins, event_2.origins):
        for key in ori_1.keys():
            if shallow and not hasattr(ori_2, key):
                continue
            if key not in ["resource_id", "comments", "arrivals",
                           "method_id", "origin_uncertainty", "depth_type",
                           "quality", "creation_info", "evaluation_mode",
                           "depth_errors", "time_errors"]:
                if ori_1[key] != ori_2[key]:
                    if verbose:
                        print('%s is not the same as %s for key %s' %
                              (ori_1[key], ori_2[key], key))
                    return False
            elif key == "arrivals" and not shallow:
                if len(ori_1[key]) != len(ori_2[key]):
                    if verbose:
                        print('Arrival: %i is not the same as %i for key %s' %
                              (len(ori_1[key]), len(ori_2[key]), key))
                    return False
                for arr_1, arr_2 in zip(ori_1[key], ori_2[key]):
                    for arr_key in arr_1.keys():
                        if arr_key not in ["resource_id", "pick_id",
                                           "distance"]:
                            if arr_1[arr_key] != arr_2[arr_key]:
                                if verbose:
                                    print('%s does not match %s for key %s' %
                                          (arr_1[arr_key], arr_2[arr_key],
                                           arr_key))
                                return False
                    if arr_1["distance"] and round(
                            arr_1["distance"]) != round(arr_2["distance"]):
                        if verbose:
                            print('%s does not match %s for key %s' %
                                  (arr_1[arr_key], arr_2[arr_key],
                                   arr_key))
                        return False
    # Check picks
    if len(event_1.picks) != len(event_2.picks):
        if verbose:
            print('Number of picks is not equal')
        return False
    event_1.picks.sort(key=lambda p: p.time)
    event_2.picks.sort(key=lambda p: p.time)
    for pick_1, pick_2 in zip(event_1.picks, event_2.picks):
        for key in pick_1.keys():
            if shallow and not hasattr(pick_2, key):
                continue
            if key not in ["resource_id", "waveform_id"]:
                if pick_1[key] != pick_2[key]:
                    # Cope with a None set not being equal to a set, but still
                    #  None quantity
                    if pick_1[key] is None and pick_2[key] is not None:
                        try:
                            if not any(pick_2[key].__dict__.values()):
                                continue
                        except AttributeError:
                            if verbose:
                                print('Pick: %s is not the same as %s for '
                                      'key %s' % (pick_1[key], pick_2[key],
                                                  key))
                            return False
                    if pick_2[key] is None and pick_1[key] is not None:
                        try:
                            if not any(pick_1[key].__dict__.values()):
                                continue
                        except AttributeError:
                            if verbose:
                                print('Pick: %s is not the same as %s for '
                                      'key %s' % (pick_1[key], pick_2[key],
                                                  key))
                            return False
                    if verbose:
                        print('Pick: %s is not the same as %s for key %s' %
                              (pick_1[key], pick_2[key], key))
                    return False
            elif key == "waveform_id":
                if pick_1[key].station_code != pick_2[key].station_code:
                    if verbose:
                        print('Station codes do not match')
                    return False
                if pick_1[key].channel_code[0] != pick_2[key].channel_code[0]:
                    if verbose:
                        print('Channel codes do not match')
                    return False
                if pick_1[key].channel_code[-1] != \
                        pick_2[key].channel_code[-1]:
                    if verbose:
                        print('Channel codes do not match')
                    return False
    # Check amplitudes
    if not len(event_1.amplitudes) == len(event_2.amplitudes):
        if verbose:
            print('Not the same number of amplitudes')
        return False
    event_1.amplitudes.sort(key=lambda a: a.generic_amplitude)
    event_2.amplitudes.sort(key=lambda a: a.generic_amplitude)
    for amp_1, amp_2 in zip(event_1.amplitudes, event_2.amplitudes):
        if shallow:
            continue
        # Assuming same ordering of amplitudes
        for key in amp_1.keys():
            if key not in ["resource_id", "pick_id", "waveform_id", "snr"]:
                if not amp_1[key] == amp_2[key] and amp_2[key] is not None:
                    if verbose:
                        print('Amplitude: %s is not the same as %s for key'
                              ' %s' % (amp_1[key], amp_2[key], key))
                    return False
            elif key == "waveform_id":
                if pick_1[key].station_code != pick_2[key].station_code:
                    if verbose:
                        print('Station codes do not match')
                    return False
                if pick_1[key].channel_code[0] != pick_2[key].channel_code[0]:
                    if verbose:
                        print('Channel codes do not match')
                    return False
                if pick_1[key].channel_code[-1] != \
                        pick_2[key].channel_code[-1]:
                    if verbose:
                        print('Channel codes do not match')
                    return False
    return True


def _par_read(dirname, compressed=True):
    """
    Internal write function to read a formatted parameter file.

    :type dirname: str
    :param dirname: Directory to read the parameter file from.
    :type compressed: bool
    :param compressed: Whether the directory is compressed or not.
    """
    from eqcorrscan.core.match_filter.matched_filter import MatchFilterError
    from eqcorrscan.core.match_filter.template import Template

    templates = []
    if compressed:
        arc = tarfile.open(dirname, "r:*")
        members = arc.getmembers()
        _parfile = [member for member in members
                    if member.name.split(os.sep)[-1] ==
                    'template_parameters.csv']
        if len(_parfile) == 0:
            arc.close()
            raise MatchFilterError(
                'No template parameter file in archive')
        parfile = arc.extractfile(_parfile[0])
    else:
        parfile = open(dirname + '/' + 'template_parameters.csv', 'r')
    for line in parfile:
        t_in = Template()
        for key_pair in line.rstrip().split(','):
            if key_pair.split(':')[0].strip() == 'name':
                t_in.__dict__[key_pair.split(':')[0].strip()] = \
                    key_pair.split(':')[-1].strip()
            elif key_pair.split(':')[0].strip() == 'filt_order':
                try:
                    t_in.__dict__[key_pair.split(':')[0].strip()] = \
                        int(key_pair.split(':')[-1])
                except ValueError:
                    pass
            else:
                try:
                    t_in.__dict__[key_pair.split(':')[0].strip()] = \
                        float(key_pair.split(':')[-1])
                except ValueError:
                    pass
        templates.append(t_in)
    parfile.close()
    if compressed:
        arc.close()
    return templates


def _resolved(x):
    return os.path.realpath(os.path.abspath(x))


def _badpath(path, base):
    """
    joinpath will ignore base if path is absolute.
    """
    return not _resolved(os.path.join(base, path)).startswith(base)


def _badlink(info, base):
    """
    Links are interpreted relative to the directory containing the link
    """
    tip = _resolved(os.path.join(base, os.path.dirname(info.name)))
    return _badpath(info.linkname, base=tip)


def _safemembers(members):
    """Check members of a tar archive for safety.
    Ensure that they do not contain paths or links outside of where we
    need them - this would only happen if the archive wasn't made by
    eqcorrscan.

    :type members: :class:`tarfile.TarFile`
    :param members: an open tarfile.
    """
    base = _resolved(".")

    for finfo in members:
        if _badpath(finfo.name, base):
            print(finfo.name, "is blocked (illegal path)")
        elif finfo.issym() and _badlink(finfo, base):
            print(finfo.name, "is blocked: Hard link to", finfo.linkname)
        elif finfo.islnk() and _badlink(finfo, base):
            print(finfo.name, "is blocked: Symlink to", finfo.linkname)
        else:
            yield finfo


def extract_from_stream(stream, detections, pad=5.0, length=30.0):
    """
    Extract waveforms for a list of detections from a stream.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream containing the detections.
    :type detections: list
    :param detections: list of eqcorrscan.core.match_filter.detection
    :type pad: float
    :param pad: Pre-detection extract time in seconds.
    :type length: float
    :param length: Total extracted length in seconds.

    :returns:
        list of :class:`obspy.core.stream.Stream`, one for each detection.
    :type: list
    """
    streams = []
    for detection in detections:
        cut_stream = Stream()
        for pick in detection.event.picks:
            tr = stream.select(station=pick.waveform_id.station_code,
                               channel=pick.waveform_id.channel_code)
            if len(tr) == 0:
                Logger.error('No data in stream for pick: {0}'.format(pick))
                continue
            cut_stream += tr.slice(
                starttime=pick.time - pad,
                endtime=pick.time - pad + length).copy()
        streams.append(cut_stream)
    return streams


def normxcorr2(template, image):
    """
    Thin wrapper to eqcorrscan.utils.correlate functions.

    :type template: numpy.ndarray
    :param template: Template array
    :type image: numpy.ndarray
    :param image:
        Image to scan the template through.  The order of these
        matters, if you put the template after the image you will get a
        reversed correlation matrix

    :return:
        New :class:`numpy.ndarray` of the correlation values for the
        correlation of the image with the template.
    :rtype: numpy.ndarray

    .. note::
        If your data contain gaps these must be padded with zeros before
        using this function. The `eqcorrscan.utils.pre_processing` functions
        will provide gap-filled data in the appropriate format.  Note that if
        you pad your data with zeros before filtering or resampling the gaps
        will not be all zeros after filtering. This will result in the
        calculation of spurious correlations in the gaps.
    """
    array_xcorr = get_array_xcorr()
    # Check that we have been passed numpy arrays
    if type(template) != np.ndarray or type(image) != np.ndarray:
        Logger.error(
            'You have not provided numpy arrays, I will not convert them')
        return 'NaN'
    if len(template) > len(image):
        ccc = array_xcorr(
            templates=np.array([image]).astype(np.float32),
            stream=template.astype(np.float32), pads=[0],
            threaded=False)[0][0]
    else:
        ccc = array_xcorr(
            templates=np.array([template]).astype(np.float32),
            stream=image.astype(np.float32), pads=[0], threaded=False)[0][0]
    ccc = ccc.reshape((1, len(ccc)))
    return ccc


if __name__ == "__main__":
    import doctest

    doctest.testmod()
