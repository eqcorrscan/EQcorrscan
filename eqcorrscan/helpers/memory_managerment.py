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


class BaseThread(Thread):
    """
    Taken from:
        https://gist.github.com/amirasaran/e91c7253c03518b8f7b7955df0e954bb
    """
    def __init__(self, callback=None, callback_args=None, *args, **kwargs):
        target = kwargs.pop('target')
        super(BaseThread, self).__init__(
            target=self.target_with_callback, *args, **kwargs)
        self.callback = callback
        self.method = target
        self.callback_args = callback_args

    def target_with_callback(self):
        self.method()
        if self.callback is not None:
            self.callback(*self.callback_args)


def callback_error(exception):
    if exception:
        raise exception


class MemoryChecker(object):
    """
    Process to run a memory check in the background.

    Inspired by: https://github.com/pythonprofilers/memory_profiler/blob/\
    master/memory_profiler.py
    """
    def __init__(self, interval=1, min_mem=1024, pid=os.getpid()):
        """
        Initializer

        :param interval: Interval to check memory available in seconds
        :param min_mem:
            Minimum memory available before raising an error in Bytes
        """
        self.interval = interval
        self.min_mem = min_mem
        self._stop = False
        self.process = psutil.Process(pid)
        self.exception = None

        self.max_mem = self.process.memory_info()[0]  # Gives RSS memory
        self.max_mem += sum(_get_child_memory(self.process))

        self.thread = BaseThread(target=self.run, callback=callback_error,
                                 callback_args=(self.exception, ), args=())
        self.thread.daemon = True
        self.thread.start()

    def run(self):
        """
        Run the memory checker.

        If min_mem remaining is reached then this will error.
        """
        while not self._stop:
            available_memory = psutil.virtual_memory().available
            used_mem = self.process.memory_info()[0] + sum(
                _get_child_memory(self.process))
            if used_mem > self.max_mem:
                self.max_mem = used_mem
            if available_memory <= self.min_mem:
                self.exception = MemoryError("Memory limit exceeded.")
                return

            time.sleep(self.interval)

    def stop(self):
        """
        Stop the memory checker
        """
        self._stop = True
        self.thread.join()
        if self.exception:
            raise self.exception


if __name__ == '__main__':
    import doctest

    doctest.testmod()
