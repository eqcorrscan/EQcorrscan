"""
Functions to enable simple energy-base triggering in a network setting where \
different stations have different noise parameters.

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

import getpass
import ast
import warnings
import numpy as np
from pprint import pprint

from multiprocessing import Pool, cpu_count
from obspy.core.util import AttribDict
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset, plot_trigger, recursive_sta_lta

import eqcorrscan
from eqcorrscan.utils.despike import median_filter


class TriggerParameters(AttribDict):
    """
    Base class for trigger parameters.

    >>> from eqcorrscan.utils.trigger import TriggerParameters
    >>> defaults = TriggerParameters()
    >>> defaults.station = 'TEST'
    >>> # Or use dictionaries
    >>> defaults['station'] = 'ALF'
    >>> defaults = TriggerParameters({'station': 'TEST',
    ...                               'channel': 'SHZ',
    ...                               'sta_len': 0.3,
    ...                               'lta_len': 10.0,
    ...                               'thr_on': 10,
    ...                               'thr_off': 3,
    ...                               'lowcut': 2,
    ...                               'highcut': 20})
    >>> print(defaults.station)
    TEST
    """
    defaults = {'station': '', 'channel': '', 'sta_len': 0.0, 'lta_len': 0.0,
                'thr_on': 0.0, 'thr_off': 0.0, 'lowcut': 0.0, 'highcut': 0.0}

    def __init__(self, header={}):
        """
        """
        super(TriggerParameters, self).__init__(header)

    def __setitem__(self, key, value):
        """
        """
        if key in ['station', 'channel']:
            # Ensure these are string type
            value = str(value)
        elif key in ['lowcut', 'highcut']:
            if value:
                value = float(value)
        else:
            value = float(value)
        if isinstance(value, dict):
            print('Setting dict')
            super(TriggerParameters, self).__setitem__(key, AttribDict(value))
        else:
            super(TriggerParameters, self).__setitem__(key, value)

    __setattr__ = __setitem__

    def __repr__(self):
        return "TriggerParameters()"

    def __str__(self):
        print_str = ' '.join(["TriggerParameters:",
                              "\n   station:", str(self.station),
                              "\n   channel:", str(self.channel),
                              "\n   sta_len:", str(self.sta_len),
                              "\n   lta_len:", str(self.lta_len),
                              "\n   thr_on:", str(self.thr_on),
                              "\n   thr_off:", str(self.thr_off),
                              "\n   lowcut:", str(self.lowcut),
                              "\n   highcut:", str(self.highcut),
                              ])
        return print_str

    def write(self, filename, append=True):
        """Write the parameters to a file as a human-readable series of dicts.

        :type filename: str
        :param filename: File to write to
        :type append: bool
        :param append: Append to already existing file or over-write.
        """
        header = ' '.join(['# User:', getpass.getuser(),
                           '\n# Creation date:', str(UTCDateTime()),
                           '\n# EQcorrscan version:',
                           str(eqcorrscan.__version__),
                           '\n\n\n'])
        if append:
            f = open(filename, 'a')
        else:
            f = open(filename, 'w')
            f.write(header)
        parameters = self.__dict__
        f.write(str(parameters))
        f.write('\n')
        f.close()
        return


def read_trigger_parameters(filename):
    """Read the trigger parameters into trigger_parameter classes.

    :type filename: str
    :param filename: Parameter file

    :returns: List of :class:`eqcorrscan.utils.trigger.TriggerParameters`
    :rtype: list

    .. rubric:: Example

    >>> from eqcorrscan.utils.trigger import read_trigger_parameters
    >>> parameters = read_trigger_parameters('parameters') # doctest: +SKIP
    """
    parameters = []
    f = open(filename, 'r')
    print('Reading parameters with the following header:')
    for line in f:
        if line[0] == '#':
            print(line.rstrip('\n').lstrip('\n'))
        else:
            parameter_dict = ast.literal_eval(line)
            # convert the dictionary to the class
            trig_par = TriggerParameters(parameter_dict)
            parameters.append(trig_par)
    f.close()
    return parameters


def _channel_loop(tr, parameters, max_trigger_length=60,
                  despike=False, debug=0):
    """
    Internal loop for parellel processing.

    :type tr: obspy.core.trace
    :param tr: Trace to look for triggers in.
    :type parameters: list
    :param parameters: List of TriggerParameter class for trace.
    :type max_trigger_length: float
    :type despike: bool
    :type debug: int

    :return: trigger
    :rtype: list
    """
    for par in parameters:
        if par['station'] == tr.stats.station and \
           par['channel'] == tr.stats.channel:
            parameter = par
            break
    else:
        msg = 'No parameters set for station ' + str(tr.stats.station)
        warnings.warn(msg)
        return []

    triggers = []
    if debug > 0:
        print(tr)
    tr.detrend('simple')
    if despike:
        median_filter(tr)
    if parameter['lowcut'] and parameter['highcut']:
        tr.filter('bandpass', freqmin=parameter['lowcut'],
                  freqmax=parameter['highcut'])
    elif parameter['lowcut']:
        tr.filter('highpass', freq=parameter['lowcut'])
    elif parameter['highcut']:
        tr.filter('lowpass', freq=parameter['highcut'])
    # find triggers for each channel using recursive_sta_lta
    df = tr.stats.sampling_rate
    cft = recursive_sta_lta(tr.data, int(parameter['sta_len'] * df),
                            int(parameter['lta_len'] * df))
    if max_trigger_length:
        trig_args = {'max_len_delete': True}
        trig_args['max_len'] = int(max_trigger_length *
                                   df + 0.5)
    if debug > 3:
        plot_trigger(tr, cft, parameter['thr_on'], parameter['thr_off'])
    tmp_trigs = trigger_onset(cft, float(parameter['thr_on']),
                              float(parameter['thr_off']),
                              **trig_args)
    for on, off in tmp_trigs:
        cft_peak = tr.data[on:off].max()
        cft_std = tr.data[on:off].std()
        on = tr.stats.starttime + \
            float(on) / tr.stats.sampling_rate
        off = tr.stats.starttime + \
            float(off) / tr.stats.sampling_rate
        triggers.append((on.timestamp, off.timestamp,
                         tr.id, cft_peak,
                         cft_std))
    return triggers


def network_trigger(st, parameters, thr_coincidence_sum, moveout,
                    max_trigger_length=60, despike=True, debug=0):
    """
    Main function to compute triggers for a network of stations.
    Computes single-channel characteristic functions using given parameters,
    then combines these to find network triggers on a number of stations
    within a set moveout window.

    :type st: obspy.core.stream.Stream
    :param st: Stream to compute detections within
    :type parameters: list
    :param parameters: List of parameter class
    :type thr_coincidence_sum: int
    :param thr_coincidence_sum:
        Minimum number of stations required to raise a network trigger.
    :type moveout: float
    :param moveout: Window to find triggers within in the network detection \
        stage.
    :type max_trigger_length: float
    :param max_trigger_length: Maximum trigger length in seconds, used to \
        remove long triggers - can set to False to not use.
    :type despike: bool
    :param despike: Whether to apply simple despiking routine or not
    :type debug: int
    :param debug: Debug output level, higher is more output.

    :returns: List of triggers
    :rtype: list

    .. rubric:: Example

    >>> from obspy import read
    >>> from eqcorrscan.utils.trigger import TriggerParameters, network_trigger
    >>> st = read("https://examples.obspy.org/" +
    ...           "example_local_earthquake_3stations.mseed")
    >>> parameters = []
    >>> for tr in st:
    ...    parameters.append(TriggerParameters({'station': tr.stats.station,
    ...                                          'channel': tr.stats.channel,
    ...                                          'sta_len': 0.5,
    ...                                          'lta_len': 10.0,
    ...                                          'thr_on': 10.0,
    ...                                          'thr_off': 3.0,
    ...                                          'lowcut': 2.0,
    ...                                          'highcut': 15.0}))
    >>> triggers = network_trigger(st=st, parameters=parameters,
    ...                            thr_coincidence_sum=5, moveout=30,
    ...                            max_trigger_length=60, despike=False)
    Looking for coincidence triggers ...
    Found 1 Coincidence triggers
    """
    triggers = []
    trace_ids = [tr.id for tr in st]
    trace_ids = dict.fromkeys(trace_ids, 1)
    if debug > 3:
        print('Not running in parallel')
        # Don't run in parallel
        for tr in st:
            triggers += _channel_loop(tr=tr, parameters=parameters,
                                      max_trigger_length=max_trigger_length,
                                      despike=despike, debug=debug)
    else:
        # Needs to be pickleable
        parameters = [par.__dict__ for par in parameters]
        pool = Pool(processes=cpu_count())
        results = [pool.apply_async(_channel_loop,
                                    args=(tr, parameters, max_trigger_length,
                                          despike, debug))
                   for tr in st]
        pool.close()
        triggers = [p.get() for p in results]
        pool.join()
        triggers = [item for sublist in triggers for item in sublist]
    triggers.sort()

    if debug > 0:
        details = True
    else:
        details = False

    print('Looking for coincidence triggers ...')
    # the coincidence triggering and coincidence sum computation
    coincidence_triggers = []
    last_off_time = 0.0
    while triggers:
        # remove first trigger from list and look for overlaps
        on, off, tr_id, cft_peak, cft_std = triggers.pop(0)
        event = {}
        event['time'] = UTCDateTime(on)
        event['stations'] = [tr_id.split(".")[1]]
        event['trace_ids'] = [tr_id]
        event['coincidence_sum'] = trace_ids[tr_id]
        if details:
            event['cft_peaks'] = [cft_peak]
            event['cft_stds'] = [cft_std]
        # compile the list of stations that overlap with the current trigger
        for trigger in triggers:
            tmp_on, tmp_off, tmp_tr_id, tmp_cft_peak, tmp_cft_std = trigger
            # skip retriggering of already present station in current
            # coincidence trigger
            if tmp_tr_id in event['trace_ids']:
                continue
            # check for overlapping trigger
            if tmp_on <= off + moveout:
                event['stations'].append(tmp_tr_id.split(".")[1])
                event['trace_ids'].append(tmp_tr_id)
                event['coincidence_sum'] += trace_ids[tmp_tr_id]
                if details:
                    event['cft_peaks'].append(tmp_cft_peak)
                    event['cft_stds'].append(tmp_cft_std)
                # allow sets of triggers that overlap only on subsets of all
                # stations (e.g. A overlaps with B and B overlaps w/ C => ABC)
                off = max(off, tmp_off)
            # break if there is a gap in between the two triggers
            else:
                break
        # skip if coincidence sum threshold is not met
        if event['coincidence_sum'] < thr_coincidence_sum:
            continue
        # skip coincidence trigger if it is just a subset of the previous
        # (determined by a shared off-time, this is a bit sloppy)
        if off == last_off_time:
            continue
        event['duration'] = off - on
        if details:
            weights = np.array([trace_ids[i] for i in event['trace_ids']])
            weighted_values = np.array(event['cft_peaks']) * weights
            event['cft_peak_wmean'] = weighted_values.sum() / weights.sum()
            weighted_values = np.array(event['cft_stds']) * weights
            event['cft_std_wmean'] = \
                (np.array(event['cft_stds']) * weights).sum() / weights.sum()
        coincidence_triggers.append(event)
        last_off_time = off

    if debug > 1:
        print('Coincidence triggers :')
        pprint(coincidence_triggers)

    print('Found %s Coincidence triggers' % len(coincidence_triggers))
    return coincidence_triggers


if __name__ == '__main__':
    import doctest
    doctest.testmod()
