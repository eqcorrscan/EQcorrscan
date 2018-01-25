"""
Functions to test the memory_management helper
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest

import os
import time
import psutil
import sys

from contextlib import contextmanager
from io import StringIO

from eqcorrscan.helpers.memory_managerment import MemoryChecker


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestMemoryChecker(unittest.TestCase):
    """Test all the options for the memory checker."""
    def test_timeout(self):
        """Check that timeout raises an error"""
        with self.assertRaises(RuntimeError):
            mem_checker = MemoryChecker(
                interval=0.05, min_mem=102.4 ** 3, timeout=0.1,
                pid=os.getpid(), target=time.sleep,
                target_args=(2.0,), target_kwargs={}, verbose=False)
            mem_checker.stop()

    def test_memory_limit(self):
        """Check that MemoryError is raised if not enough RAM is left."""
        available_memory = psutil.virtual_memory().available
        with self.assertRaises(MemoryError):
            mem_checker = MemoryChecker(
                interval=0.05, min_mem=available_memory - 100, timeout=20,
                pid=os.getpid(), target=time.sleep,
                target_args=(1.0,), target_kwargs={}, verbose=False)
            mem_checker.stop()

    def test_no_error(self):
        mem_checker = MemoryChecker(
                interval=0.01, min_mem=102.4 ** 3, timeout=1.0,
                pid=os.getpid(), target=time.sleep,
                target_args=(0.1,), target_kwargs={}, verbose=False)
        mem_checker.stop()
        self.assertAlmostEqual(mem_checker.runtime, 0.1, 1)

    def test_timeout_verbose(self):
        """Check that timeout raises an error"""
        with captured_output() as (out, err):
            try:
                mem_checker = MemoryChecker(
                    interval=0.05, min_mem=102.4 ** 3, timeout=0.1,
                    pid=os.getpid(), target=time.sleep,
                    target_args=(2.0,), target_kwargs={}, verbose=True)
                mem_checker.stop()
            except RuntimeError:
                pass
        self.assertTrue("Timeout error" in out.getvalue().strip())

    def test_memory_limit_verbose(self):
        """Check that MemoryError is raised if not enough RAM is left."""
        available_memory = psutil.virtual_memory().available
        with captured_output() as (out, err):
            try:
                mem_checker = MemoryChecker(
                    interval=0.05, min_mem=available_memory - 100, timeout=20,
                    pid=os.getpid(), target=time.sleep,
                    target_args=(1.0,), target_kwargs={}, verbose=True)
                mem_checker.stop()
            except MemoryError:
                pass
        self.assertTrue("Memory error" in out.getvalue().strip())

    def test_no_error_verbose(self):
        with captured_output() as (out, err):
            mem_checker = MemoryChecker(
                    interval=0.01, min_mem=102.4 ** 3, timeout=1.0,
                    pid=os.getpid(), target=time.sleep,
                    target_args=(0.1,), target_kwargs={}, verbose=True)
            mem_checker.stop()
        self.assertTrue("Runtime:" in out.getvalue().strip())


if __name__ == '__main__':
    unittest.main()
