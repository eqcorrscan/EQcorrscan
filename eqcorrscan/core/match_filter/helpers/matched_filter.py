"""
Helper functions for network matched-filter detection of seismic data.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import os
import logging
import pickle

import numpy as np

from concurrent.futures import ThreadPoolExecutor
from obspy import Stream
from obspy.core.event import Event

from eqcorrscan.utils.correlate import get_array_xcorr


Logger = logging.getLogger(__name__)


def _tr_spike_test(data, percent, multiplier):
    data_len = data.shape[0]
    thresh = 2 * np.max(np.sort(
        np.abs(data))[0:np.int64(percent * data_len)]) * multiplier
    if (data > thresh).sum() > 0:
        return True
    return False


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

    to_check = ((tr.data, percent, multiplier) for tr in stream)
    list_ids = []
    with ThreadPoolExecutor() as executor:
        for tr, spiked in zip(stream, executor.map(
                lambda args: _tr_spike_test(*args), to_check)):
            if spiked:
                list_ids.append(tr.id)

    if len(list_ids):
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


def _remove_duplicates(party):
    for family in party:
        if family is not None:
            # Slow uniq:
            # family.detections = family._uniq().detections
            # Very quick uniq:
            det_tuples = [
                (det.id, str(det.detect_time), det.detect_val)
                for det in family]
            # Retrieve the indices for the first occurrence of each
            # detection in the family (so only unique detections will
            # remain).
            uniq_det_tuples, uniq_det_indices = np.unique(
                det_tuples, return_index=True, axis=0)
            uniq_detections = []
            for uniq_det_index in uniq_det_indices:
                uniq_detections.append(family[uniq_det_index])
            family.detections = uniq_detections
    return party


def _moveout(st: Stream) -> float:
    """ Maximum moveout across template in seconds. """
    return max(tr.stats.starttime for tr in st) - min(
        tr.stats.starttime for tr in st)


def _mad(cccsum):
    """
    Internal helper to compute MAD-thresholds in parallel.
    """
    return np.median(np.abs(cccsum))


def _pickle_stream(stream: Stream, filename: str):
    Logger.info(f"Pickling stream of {len(stream)} traces to {filename}")
    with open(filename, "wb") as f:
        pickle.dump(stream, f)
    Logger.info(f"Pickled to {filename}")
    return


def _unpickle_stream(filename: str) -> Stream:
    Logger.info(f"Unpickling from {filename}")
    with open(filename, "rb") as f:
        stream = pickle.load(f)
    assert isinstance(stream, Stream)
    Logger.info(f"Unpickled stream of {len(stream)} traces from {filename}")
    return stream


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
