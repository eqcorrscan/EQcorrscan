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
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools

from obspy import Trace, Stream
from multiprocessing import cpu_count

from eqcorrscan.utils.correlate import (
    XCOR_FUNCS, XCORR_STREAM_METHODS, get_stream_xcorr)
from eqcorrscan.utils import correlate
from eqcorrscan.utils.parameters import EQcorrscanConfig, CorrelationDefaults
from eqcorrscan.helpers.memory_managerment import MemoryChecker


MIN_MEM = 102.4 ** 3  # Limit minimum memory available to 0.1 GB


# TODO: consider whether randn is good approx of data for timing purposes
def generate_dataset(n_templates, n_stations, n_channels, data_len,
                     template_len, sampling_rate):
    """
    Generate a synthetic dataset to test correlations within

    :param n_templates: Number of templates to create
    :param n_stations: Number of stations to create data for
    :param n_channels: Number of channels for each station to create data for
    :param data_len: Length of continuous data in seconds
    :param template_len: Length of templates in seconds
    :param sampling_rate: Sampling-rate for all data in Hz
    :return: Dictionary of dataset keyed by "templates" and "data"
    """
    print("=" * 80)
    print("Generating dataset".center(80))
    print("=" * 80)
    dataset = {'data': Stream(), 'templates': []}
    for i in range(n_stations * n_channels):
        _tr = np.random.randn(int(data_len * sampling_rate))
        _tr = _tr * _tr ** 10
        dataset['data'] += Trace(_tr)
    for t in range(n_templates):
        template = Stream()
        for i in range(n_stations * n_channels):
            _tr = np.random.randn(int(template_len * sampling_rate))
            _tr = _tr * _tr ** 4
            template += Trace(_tr)
        dataset['templates'].append(template)
    return dataset


def run_correlation(func, threads, dataset, loops=3, timeout=600):
    """
    Run the given correlation function and profile time and memory usage.

    :param func:
        Callable correlation function with the standard EQcorrscan signature.
    :param threads: Number of threads to run on for concurrent methods
    :param dataset: Dictionary of dataset including templates and data
    :param loops: Number of loops to average over.
    :param timeout: Timeout limit for correlations - defaults to 600

    :return: tuple of average memory and average time in GB and S respectively.
    """
    max_mem = []
    timings = []
    for loop in range(loops):
        overrun = False
        overtime = False
        tic = time.time()
        try:
            mem_checker = MemoryChecker(
                interval=0.05, min_mem=MIN_MEM, timeout=timeout,
                pid=os.getpid(), target=func,
                target_args=(dataset['templates'], dataset['data']),
                target_kwargs={'cores': threads}, verbose=False)
            mem_checker.stop()
        except MemoryError as e:
            overrun = True
            print("MemoryError: " + str(e))
        except RuntimeError as e:
            overtime = True
            print("RuntimeError: " + str(e))
        toc = time.time()
        max_mem.append(mem_checker.max_mem / (1024 ** 3))
        if not overrun and not overtime:
            timings.append(toc - tic)
        else:
            timings.append(np.nan)
    return np.mean(max_mem), np.mean(timings)


def plot_profiles(times, memory_use):
    """
    Plot profiling information.

    :param times: Dictionary of times in seconds, keyed by method
    :param memory_use: Dictionary of memory usage in GB keyed by method
    """
    plt.xkcd()
    colors = itertools.cycle(('k', 'r', 'b', 'g', 'c', 'y'))
    keys = sorted(list(times.keys()))
    methods = list(set([key.split('.')[0] for key in keys]))
    concurrencies = list(set([key.split('.')[1] for key in keys]))
    for method in methods:
        color = next(colors)
        markers = itertools.cycle((',', '+', '.', 'o', '*'))
        for concurrency in concurrencies:
            key = method + '.' + concurrency
            plt.scatter(memory_use[key], times[key], s=10 ** 2, c=color,
                        marker=next(markers), label=key)
    plt.xlabel('Memory use (GB)')
    plt.ylabel('Time (s)')
    plt.legend()
    plt.show()


def run_profiling(n_templates, n_stations, n_channels, data_len,
                  template_len, sampling_rate, loops=3, timeout=600):
    """
    Run profiling for available correlation functions and write config file.

    :param n_templates: Number of templates to create
    :param n_stations: Number of stations to create data for
    :param n_channels: Number of channels for each station to create data for
    :param data_len: Length of continuous data in seconds
    :param template_len: Length of templates in seconds
    :param sampling_rate: Sampling-rate for all data in Hz
    :param loops: Number of loops to average times and memory over.
    :param timeout: Timeout limit for correlations - defaults to 600
    """
    print("Found EQcorrscan correlation functions in %s" % correlate.__file__)
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
                mem, avtime = run_correlation(
                    func, MAXTHREADS, dataset, loops, timeout)
                times.update({'.'.join([corr_func, method]): avtime})
                memory_use.update({'.'.join([corr_func, method]): mem})
                print(("Average time from %i loops: %f seconds" %
                       (loops, avtime)).center(80))
                print(("Average Max Memory: %f GB" % (mem)).center(80))
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
    config.uniq().write(append=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Profile correlation functions to find out the fastest for"
                    " your system and dataset size.")
    parser.add_argument(
        '-t', '--n-templates', help="Number of templates for test dataset",
        required=True)
    parser.add_argument(
        '-s', '--n-stations', help="Number of stations for test dataset",
        required=True)
    parser.add_argument(
        '-c', '--n-channels',
        help="Number of channels for each station in test dataset",
        required=True)
    parser.add_argument(
        '-d', '--data-len', help="Length of continuous test data in seconds",
        required=True)
    parser.add_argument(
        '-l', '--template-len', help="Length of templates in seconds",
        required=True)
    parser.add_argument(
        '-r', '--sampling-rate', help="Sampling-rate in Hz for test dataset",
        required=True)
    parser.add_argument(
        '-p', '--loops', help="Number of loops to average over", required=True)
    parser.add_argument(
        '--timeout', help="Timeout limit for correlations - defaults to 600",
        required=False)
    args = vars(parser.parse_args())
    if args['timeout'] is not None:
        timeout = float(args['timeout'])
    else:
        timeout = 600
    run_profiling(
        n_templates=int(args['n_templates']),
        n_stations=int(args['n_stations']),
        n_channels=int(args['n_channels']), data_len=float(args['data_len']),
        template_len=float(args['template_len']),
        sampling_rate=float(args['sampling_rate']), loops=int(args['loops']),
        timeout=timeout)
