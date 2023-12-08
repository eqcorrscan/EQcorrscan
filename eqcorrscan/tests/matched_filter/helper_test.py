"""
Test the helpers for EQcorrscan's matched-filter codes.
"""

import abc
import copy
import os.path
import time
import shutil
import tempfile
import logging
import unittest

from multiprocessing import Process, Queue

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.clients.earthworm import Client as EWClient

from eqcorrscan.core.match_filter.helpers import get_waveform_client
from eqcorrscan.core.match_filter.helpers.processes import (
    _get_detection_stream, Poison)


Logger = logging.getLogger(__name__)


class TestHelperFunctions(unittest.TestCase):
    def test_monkey_patching(self):
        """ Test that monkey patching a client works. """
        client = EWClient("pubavo1.wr.usgs.gov", 16022)
        self.assertFalse(hasattr(client, "get_waveforms_bulk"))
        client = get_waveform_client(client)
        self.assertTrue(hasattr(client, "get_waveforms_bulk"))


class ProcessTests(abc.ABC, unittest.TestCase):
    directories_to_nuke = []
    wait_time = 5
    kwargs = dict()
    process = None

    @classmethod
    def tearDownClass(cls):
        for _dir in cls.directories_to_nuke:
            try:
                shutil.rmtree(_dir)
            except FileNotFoundError:
                pass


class TestGetDetectionStreamProcess(ProcessTests):
    directories_to_nuke = []
    wait_time = 5

    @classmethod
    def setUpClass(cls):
        process_length = 360
        cls.global_kwargs = dict(
            template_channel_ids=[("NZ", "WVZ", "10", "HHZ")],
            client=Client("GEONET"),
            retries=3,
            min_gap=0.0,
            buff=3,
            full_stream_dir=None,
            pre_process=False,
            parallel_process=True,
            process_cores=1,
            daylong=False,
            overlap=3.0,
            ignore_length=False,
            ignore_bad_data=False,
            filt_order=4,
            highcut=20,
            lowcut=2,
            samp_rate=50,
            process_length=process_length,
        )

    def setUp(self):
        # Make a copy of the class-wide kwargs
        self.kwargs = copy.copy(self.global_kwargs)
        self.kwargs.update(
            dict(input_time_queue=Queue(),
                 poison_queue=Queue(),
                 output_filename_queue=Queue(),
                 temp_stream_dir=tempfile.mkdtemp()))
        Logger.info(self.kwargs)
         # Cleanup
        self.directories_to_nuke.append(self.kwargs['temp_stream_dir'])

    def test_poisoning(self):
        self.process = Process(
            target=_get_detection_stream, kwargs=self.kwargs,
            name="TestProcess")
        poisoning(obj=self)

    def test_poisoning_while_waiting_on_output(self):
        self.process = Process(
            target=_get_detection_stream, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_while_waiting_on_output(obj=self)

    def test_poisoning_from_input(self):
        self.process = Process(
            target=_get_detection_stream, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_from_input(obj=self)

    def test_normal_operation(self):
        self.process = Process(
            target=_get_detection_stream, kwargs=self.kwargs,
            name="TestProcess")
        self.process.start()

        # Populate time queue
        self.kwargs['input_time_queue'].put(
            (UTCDateTime(2021, 1, 1),
             UTCDateTime(2021, 1, 1, 0, 10)))
        self.kwargs['input_time_queue'].put(None)

        # Get the output
        filename = self.kwargs['output_filename_queue'].get()
        self.assertTrue(os.path.isfile(filename))
        self.assertEqual(self.kwargs['output_filename_queue'].get(), None)

        # Wait for the process to end
        time.sleep(self.wait_time)
        self.assertFalse(self.process.is_alive())

    def test_full_stream_operation(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            full_stream_dir=tempfile.mkdtemp(),
            pre_process=True,))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs,
            name="TestProcess")
        process.start()

        # Populate time queue
        kwargs['input_time_queue'].put(
            (UTCDateTime(2021, 1, 1),
             UTCDateTime(2021, 1, 1, 0, 10)))
        kwargs['input_time_queue'].put(None)

        # Get the output
        filename = kwargs['output_filename_queue'].get()
        self.assertTrue(os.path.isfile(filename))
        self.assertEqual(kwargs['output_filename_queue'].get(), None)

        # Wait for the process to end
        time.sleep(self.wait_time)
        self.assertFalse(process.is_alive())

        # Check for full stream
        full_st = read(f"{kwargs['full_stream_dir']}/*")
        st = read(filename)

        self.assertEqual(st[0].id, full_st[0].id)
        self.assertEqual(st[0].stats.sampling_rate, kwargs['samp_rate'])
        self.assertEqual(full_st[0].stats.sampling_rate, 100.0)

        # Cleanup
        self.directories_to_nuke.append(kwargs['full_stream_dir'])

    def test_exception_handling(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            overlap="bob",  # This need to be a float! Should raise exception
            pre_process=True,
        ))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        process.start()

        # Populate time queue
        kwargs['input_time_queue'].put(
            (UTCDateTime(2021, 1, 1),
             UTCDateTime(2021, 1, 1, 0, 10)))
        kwargs['input_time_queue'].put(None)

        time.sleep(self.wait_time)
        poison = kwargs['poison_queue'].get()
        self.assertIsInstance(poison, Poison)
        time.sleep(self.wait_time)
        self.assertFalse(process.is_alive())


################################################################################
#  STANDARD PROCESS DEATH TESTS
################################################################################


def poisoning(obj: ProcessTests):
    obj.process.start()
    # Test death
    obj.kwargs['poison_queue'].put(Exception("TestException"))
    time.sleep(obj.wait_time)
    obj.assertFalse(obj.process.is_alive())
    obj.directories_to_nuke.append(obj.kwargs['temp_stream_dir'])


def poisoning_while_waiting_on_output(obj: ProcessTests):
    obj.process.start()
    # Test death
    obj.kwargs['poison_queue'].put(Exception("TestException"))
    time.sleep(obj.wait_time)
    obj.assertFalse(obj.process.is_alive())
    obj.directories_to_nuke.append(obj.kwargs['temp_stream_dir'])


def poisoning_from_input(obj):
    # Test death
    obj.kwargs['input_time_queue'].put(Poison(Exception("TestException")))
    time.sleep(obj.wait_time)
    obj.assertFalse(obj.process.is_alive())
    # Cleanup
    shutil.rmtree(obj.kwargs['temp_stream_dir'])


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
