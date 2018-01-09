"""
This script will construct a synthetic dataset of specified size and run the
various correlation functions available in EQcorrscan over the data for
parallelisms compatible with your machine.  In doing so timings and memory
usage will be computed to provide information to the user as to what the best
algorithm for their dataset and architecture is.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import time
import resource
import sys
import numpy as np
import matplotlib.pyplot as plt
import itertools

from memory_profiler import memory_usage
from obspy import Trace, Stream
from multiprocessing import cpu_count

from eqcorrscan.utils.correlate import (
    XCOR_FUNCS, XCORR_STREAM_METHODS, get_stream_xcorr)
from eqcorrscan.utils import correlate
from eqcorrscan.utils.parameters import EQcorrscanConfig, CorrelationDefaults

#TODO: Need to find a cross-platform way to limit memory
def memory_limit():
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024, hard))
    # resource is Unix specific


def get_memory():
    with open('/proc/meminfo', 'r') as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                free_memory += int(sline[1])
    return free_memory


def generate_dataset(n_templates, n_stations, n_channels, data_len,
                     template_len, sampling_rate):
    """

    :param n_templates:
    :param n_stations:
    :param n_channels:
    :param data_len:
    :param template_len:
    :param sampling_rate:
    :return:
    """
    print("=" * 80)
    print("Generating dataset".center(80))
    print("=" * 80)
    dataset = {'data': Stream(), 'templates': []}
    for i in range(n_stations * n_channels):
        dataset['data'] += Trace(np.random.randn(
            int(data_len * sampling_rate)))
    for t in range(n_templates):
        template = Stream()
        for i in range(n_stations * n_channels):
            template += Trace(np.random.randn(
                int(template_len * sampling_rate)))
        dataset['templates'].append(template)
    return dataset


def run_correlation(func, threads, dataset, loops=3):
    """
    Run the given correlation function and profile time and memory usage.

    :param func:
    :param threads:
    :param dataset:
    :param loops:
    :return:
    """
    max_mem = []
    timings = []
    for loop in range(loops):
        tic = time.time()
        max_mem = memory_usage((func, (dataset['templates'], dataset['data']),
                                {'cores': threads}), max_usage=True,
                               include_children=True, multiprocess=True)
        toc = time.time()
        timings.append(toc - tic)
    return np.mean(max_mem[0]), np.mean(timings)


def plot_profiles(times, memory_use):
    """

    :param times:
    :param memory_use:
    :return:
    """
    plt.xkcd()
    markers = itertools.cycle((',', '+', '.', 'o', '*'))
    for key in times.keys():
        plt.scatter(memory_use[key], times[key], s=10**2, marker=next(markers),
                    label=key)
    plt.xlabel('Memory use (MB)')
    plt.ylabel('Time (s)')
    plt.legend()
    plt.show()


def run_profiling(n_templates, n_stations, n_channels, data_len,
                  template_len, sampling_rate, loops=3):
    """

    :param n_templates:
    :param n_stations:
    :param n_channels:
    :param data_len:
    :param template_len:
    :param sampling_rate:
    :return:
    """
    print("Testing correlation functions in %s" % correlate.__file__)
    MAXTHREADS = cpu_count()
    dataset_size = {
        'n_templates': n_templates, 'n_stations': n_stations,
        'n_channels': n_channels, 'data_len': data_len,
        'template_len': template_len, 'sampling_rate': sampling_rate}
    # Check if this test has already been run
    config = EQcorrscanConfig()
    if config.get("correlation") is not None:
        for default in config.get("correlation"):
            if default.dataset_size == dataset_size:
                print("Already ran this dataset size, default set to %s" %
                      default.corr_func)
    dataset = generate_dataset(
        n_templates=n_templates, n_stations=n_stations, n_channels=n_channels,
        data_len=data_len, template_len=template_len,
        sampling_rate=sampling_rate)
    times = {}
    memory_use = {}
    best_time = {'None': np.inf}
    # Limit the memory
    memory_limit()
    for corr_func in XCOR_FUNCS.keys():
        if corr_func == 'default':
            continue
        print("=" * 80)
        print(("Running %s" % corr_func).center(80))
        print("=" * 80)
        for method in XCORR_STREAM_METHODS:
            try:
                print(("Testing %s method" % method).center(80))
                func = get_stream_xcorr(corr_func, method)
                mem, avtime = run_correlation(func, MAXTHREADS, dataset, loops)
                times.update({'.'.join([corr_func, method]): avtime})
                memory_use.update({'.'.join([corr_func, method]): mem})
                print(("Average time from %i loops: %f seconds" %
                       (loops, avtime)).center(80))
                print(("Average Max Memory: %f MB" % (mem)).center(80))
                print('-' * 80)
                if avtime < list(best_time.values())[0]:
                    best_time = {'.'.join([corr_func, method]): avtime}
            except MemoryError:
                print("Exceeded maximum memory allowed")
    plot_profiles(times, memory_use)
    # Write config to file
    best_correlation = CorrelationDefaults(corr_func=list(best_time.keys())[0])
    best_correlation.__dict__.update(dataset_size)
    config.defaults.append(best_correlation)
    config.write()


if __name__ == '__main__':
    #TODO: Use something nicer for getting args and giving feedback/help
    # Get arguments / print help
    help_msg = ("Profile correlation functions. Needs an ordered list of "
                "arguments as follows:\n"
                "Usage: python correlation_speeds.py n_templates n_stations "
                "n_channels data_len template_len sampling_rate loops")
    if len(sys.argv) == 1 or sys.argv[1] == "-h":
        print(help_msg)
        exit()
    try:
        run_profiling(
            n_templates=int(sys.argv[1]), n_stations=int(sys.argv[2]),
            n_channels=int(sys.argv[3]), data_len=float(sys.argv[4]),
            template_len=float(sys.argv[5]), sampling_rate=float(sys.argv[6]),
            loops=int(sys.argv[7]))
    except IndexError:
        print("Insufficient arguments.")
        print(help_msg)
