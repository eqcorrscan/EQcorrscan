"""
Test the helpers for EQcorrscan's matched-filter codes.
"""

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


class TestGetDetectionStreamProcess(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wait_time = 1
        process_length = 360
        cls.kwargs = dict(
            template_channel_ids=[("NZ", "WVZ", "10", "HHZ")],
            client=Client("GEONET"),
            # input_time_queue=input_queue,
            retries=3,
            min_gap=0.0,
            buff=3,
            # output_filename_queue=output_queue,
            # poison_queue=poison_queue,
            # temp_stream_dir=temp_stream_dir,
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
            process_length=process_length)
        # Queues need to be managed by individual tests
        # input_queue = Queue()
        # poison_queue = Queue()
        # output_queue = Queue()
        # temp_stream_dir = tempfile.mkdtemp()

    def test_poisoning(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            input_time_queue=Queue(),
            poison_queue=Queue(),
            output_filename_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp()))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")

        process.start()

        # Test death
        kwargs['poison_queue'].put(Exception("TestException"))
        time.sleep(self.wait_time)

        self.assertFalse(process.is_alive())

        # Cleanup
        shutil.rmtree(kwargs['temp_stream_dir'])

    def test_poisoning_while_waiting_on_output(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            input_time_queue=Queue(),
            poison_queue=Queue(),
            output_filename_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp()))
        kwargs['output_filename_queue'].put("this_is_a_test")
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        process.start()
        # Test death
        kwargs['poison_queue'].put(Exception("TestException"))
        time.sleep(self.wait_time)
        self.assertFalse(process.is_alive())
        # Cleanup
        shutil.rmtree(kwargs['temp_stream_dir'])

    def test_poisoning_from_input(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            input_time_queue=Queue(),
            poison_queue=Queue(),
            output_filename_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp()))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        process.start()
        # Test death
        kwargs['input_time_queue'].put(Poison(Exception("TestException")))
        time.sleep(self.wait_time)
        self.assertFalse(process.is_alive())
        # Cleanup
        shutil.rmtree(kwargs['temp_stream_dir'])

    def test_normal_operation(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            input_time_queue=Queue(),
            poison_queue=Queue(),
            output_filename_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp()))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        process.start()

        # Populate time queue
        kwargs['input_time_queue'].put((UTCDateTime(2021, 1, 1),
                                        UTCDateTime(2021, 1, 1, 0, 10)))
        kwargs['input_time_queue'].put(None)

        # Get the output
        filename = kwargs['output_filename_queue'].get()
        self.assertTrue(os.path.isfile(filename))
        self.assertEqual(kwargs['output_filename_queue'].get(), None)

        # Wait for the process to end
        time.sleep(self.wait_time)
        self.assertFalse(process.is_alive())
        # Cleanup
        shutil.rmtree(kwargs['temp_stream_dir'])

    def test_full_stream_operation(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            input_time_queue=Queue(),
            poison_queue=Queue(),
            output_filename_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp(),
            full_stream_dir=tempfile.mkdtemp(),
            pre_process=True,
        ))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        process.start()

        # Populate time queue
        kwargs['input_time_queue'].put((UTCDateTime(2021, 1, 1),
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
        shutil.rmtree(kwargs['temp_stream_dir'])
        shutil.rmtree(kwargs['full_stream_dir'])

    def test_exception_handling(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            input_time_queue=Queue(),
            poison_queue=Queue(),
            output_filename_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp(),
            overlap="bob",  # This need to be a float! Should raise exception
            pre_process=True,
        ))
        process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        process.start()

        # Populate time queue
        kwargs['input_time_queue'].put((UTCDateTime(2021, 1, 1),
                                        UTCDateTime(2021, 1, 1, 0, 10)))
        kwargs['input_time_queue'].put(None)

        time.sleep(self.wait_time)
        poison = kwargs['poison_queue'].get()
        self.assertIsInstance(poison, Poison)
        self.assertFalse(process.is_alive())
        # Cleanup
        shutil.rmtree(kwargs['temp_stream_dir'])


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
