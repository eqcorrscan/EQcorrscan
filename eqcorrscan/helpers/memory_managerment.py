"""
Helpers to manage memory usage.

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

from threading import Thread
from multiprocessing import Process
import psutil
import os
import time


def _get_child_memory(process):
    """
    Returns a generator that yields memory for all child processes.

    Adapted from https://github.com/pythonprofilers/memory_profiler/blob/\
    master/memory_profiler.py
    """
    # Convert a pid to a process
    if isinstance(process, int):
        if process == -1:
            process = os.getpid()
        process = psutil.Process(process)

    # Select the psutil function get the children similar to how we selected
    # the memory_info attr (a change from excepting the AttributeError).
    children_attr = (
        'children' if hasattr(process, 'children') else 'get_children')

    # Loop over the child processes and yield their memory
    try:
        for child in getattr(process, children_attr)(recursive=True):
            yield child.memory_info()[0]
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        # https://github.com/fabianp/memory_profiler/issues/71
        yield 0.0


class MemoryChecker(object):
    """
    Process to run a memory check in the background.

    Inspired by: https://github.com/pythonprofilers/memory_profiler/blob/\
    master/memory_profiler.py
    """
    def __init__(self, interval=1, min_mem=1024, timeout=None, pid=os.getpid(),
                 target=None, target_args=(), target_kwargs={}, verbose=False):
        """
        Initializer

        :param interval: Interval to check memory available in seconds
        :param min_mem:
            Minimum memory available before raising an error in Bytes
        :param timeout: Timeout period in seconds
        :param target: function to run and check memory and timing
        :param target_args: args for function
        :param target_kwargs: kwargs for function
        :param verbose: Output memory and timing as it goes.
        """
        self.interval = interval
        self.min_mem = min_mem
        self._stop = False
        self.process = psutil.Process(pid)
        self.exception = None  # set to None when the threads starts
        self.starttime = time.time()
        self.runtime = 0
        if timeout is None:
            timeout = float("inf")
        self.timeout = timeout
        self.verbose = verbose

        self.max_mem = self.process.memory_info()[0]  # Gives RSS memory
        self.max_mem += sum(_get_child_memory(self.process))

        self.checker_thread = Thread(target=self.run)
        self.checker_thread.start()
        self.worker_process = Process(
            target=target, args=target_args, kwargs=target_kwargs)
        self.worker_process.start()

    def run(self):
        """
        Run the memory checker.
        """
        try:
            while not self._stop:
                self.runtime = time.time() - self.starttime
                available_memory = psutil.virtual_memory().available
                used_mem = self.process.memory_info()[0] + sum(
                    _get_child_memory(self.process))
                if used_mem > self.max_mem:
                    self.max_mem = used_mem

                if self.verbose:
                    print("Runtime: {0:4.2f}S\tMemory used: {1:4.2f}GB\t"
                          " Memory available {2:4.2f}GB".format(
                           self.runtime, self.max_mem / 1024 ** 3,
                           available_memory / 1024 ** 3),
                          end="\n")
                if available_memory <= self.min_mem:
                    self.exception = MemoryError(
                        "%f <= %f: Memory limit exceeded." %
                        (available_memory, self.min_mem))
                    if self.verbose:
                        print("Memory error - returning. Memory "
                              "available=%fB" % available_memory)
                    self.worker_process.terminate()
                if self.runtime >= self.timeout:
                    self.exception = RuntimeError(
                        "%f >= %f: Timeout limit exceeded." %
                        (self.runtime, self.timeout))
                    if self.verbose:
                        print("Timeout error - returning. Runtime=%fs"
                              % self.runtime)
                    self.worker_process.terminate()
                    return

                time.sleep(self.interval)
        except BaseException as e:
            self.exception = e
            self.worker_process.terminate()
            return

    def stop(self):
        """
        Stop the memory checker
        """
        self.worker_process.join()
        self._stop = True
        self.checker_thread.join()
        if self.exception:
            raise self.exception
        if self.verbose:
            print("\n")


if __name__ == '__main__':
    import doctest

    doctest.testmod()
