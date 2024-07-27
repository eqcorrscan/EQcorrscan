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
import pytest

from multiprocessing import Process, Queue

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
from obspy.clients.earthworm import Client as EWClient

from eqcorrscan.core.match_filter import Party
from eqcorrscan.core.match_filter.helpers import (
    get_waveform_client, _test_event_similarity)
from eqcorrscan.core.match_filter.helpers.processes import (
    _get_detection_stream, Poison, _prepper, _pre_processor,
    _make_detections)

try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:
    pass
else:
    cleanup_on_sigterm()


Logger = logging.getLogger(__name__)
MAX_WAIT = 15  # Maximum wait time for process to close


def get_test_templates():
    party = Party().read(
        filename=os.path.join(
            os.path.abspath(os.path.dirname(os.path.dirname(__file__))),
            'test_data', 'test_party.tgz'))
    return [f.template for f in party]


class TestHelperFunctions(unittest.TestCase):
    def test_monkey_patching(self):
        """ Test that monkey patching a client works. """
        client = EWClient("pubavo1.wr.usgs.gov", 16022)
        self.assertFalse(hasattr(client, "get_waveforms_bulk"))
        client = get_waveform_client(client)
        self.assertTrue(hasattr(client, "get_waveforms_bulk"))

    @pytest.mark.network
    def test_event_similarity_quiet(self):
        self._event_similarity(verbose=False)

    @pytest.mark.network
    def test_event_similarity_loud(self):
        self._event_similarity(verbose=True)

    def _event_similarity(self, verbose: bool = False):
        client = Client("GEONET")
        event = client.get_events(eventid="2023p923930")[0]
        self.assertTrue(_test_event_similarity(
            event, event, verbose=verbose))
        event2 = event.copy()
        self.assertTrue(_test_event_similarity(
            event, event2, verbose=verbose))
        with self.assertRaises(NotImplementedError):
            _test_event_similarity(event, "bob")
        event2.origins = event2.origins[0:-2]
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.origins[-1].arrivals = event2.origins[-1].arrivals[0:-2]
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.origins[-1].arrivals[-1].time_residual = \
            event2.origins[-1].arrivals[-1].time_residual - .5
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.origins[-1].arrivals[-1].distance = \
            event2.origins[-1].arrivals[-1].distance - 1.5
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.origins[-1].time -= 60
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        # Picks
        event2 = event.copy()
        event2.picks = event2.picks[0:-2]
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.picks[0].time += 20
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.picks[0].waveform_id.station_code = "BOB"
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.picks[0].waveform_id.channel_code = "BOB"
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        # Amplitudes
        event2 = event.copy()
        event2.amplitudes = event2.amplitudes[0:-2]
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.amplitudes[0].generic_amplitude += 100
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.amplitudes[0].waveform_id.station_code = "BOB"
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))
        event2 = event.copy()
        event2.amplitudes[0].waveform_id.channel_code = "BOB"
        self.assertFalse(_test_event_similarity(
            event, event2, verbose=verbose))


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

    def wait_for_close(self):
        total_wait = 0
        while total_wait <= MAX_WAIT:
            if self.process.is_alive():
                time.sleep(self.wait_time)
            else:
                break
        self.assertFalse(self.process.is_alive())
        self.process.join()


class TestMakeDetections(ProcessTests):
    @classmethod
    def setUpClass(cls):
        templates = get_test_templates()
        cls.global_kwargs = dict(
            delta=0.01,
            templates=templates,
            threshold=8.0,
            threshold_type="MAD",
            save_progress=False,
        )

    def setUp(self):
        self.kwargs = copy.copy(self.global_kwargs)
        self.kwargs.update(dict(
            input_queue=Queue(),
            output_queue=Queue(),
            poison_queue=Queue()
        ))

    def test_poisoning(self):
        self.process = Process(
            target=_make_detections, kwargs=self.kwargs,
            name="TestProcess")
        poisoning(obj=self)

    def test_poisoning_from_input(self):
        self.process = Process(
            target=_make_detections, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_from_input(
            obj=self, input_queue=self.kwargs['input_queue'])


class TestPrepper(ProcessTests):
    @classmethod
    def setUpClass(cls):
        templates = get_test_templates()
        cls.global_kwargs = dict(
            templates=templates,
            group_size=5,
            groups=None,
            xcorr_func="fmf",
        )

    def setUp(self):
        self.kwargs = copy.copy(self.global_kwargs)
        self.kwargs.update(dict(
            input_stream_filename_queue=Queue(),
            output_queue=Queue(maxsize=1),
            poison_queue=Queue()
        ))

    def test_poisoning(self):
        self.process = Process(
            target=_prepper, kwargs=self.kwargs,
            name="TestProcess")
        poisoning(obj=self)

    def test_poisoning_while_waiting_on_output(self):
        self.process = Process(
            target=_prepper, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_while_waiting_on_output(
            obj=self, output_queue=self.kwargs["output_queue"])

    def test_poisoning_from_input(self):
        self.process = Process(
            target=_prepper, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_from_input(
            obj=self, input_queue=self.kwargs['input_stream_filename_queue'])


class TestPreProcessor(ProcessTests):
    @classmethod
    def setUpClass(cls):
        process_length = 360
        cls.global_kwargs = dict(
            template_ids={("NZ", "WVZ", "10", "HHZ")},
            pre_processed=False,
            ignore_length=False,
            ignore_bad_data=False,
            filt_order=4,
            highcut=20,
            lowcut=2,
            samp_rate=50,
            process_length=process_length,
            parallel=False,
            cores=1,
            overlap=3.0,
        )

    def setUp(self):
        self.kwargs = copy.copy(self.global_kwargs)
        self.kwargs.update(dict(
            input_stream_queue=Queue(),
            temp_stream_dir=tempfile.mkdtemp(),
            output_filename_queue=Queue(maxsize=1),
            poison_queue=Queue()
        ))
        Logger.info(self.kwargs)
        self.directories_to_nuke.append(self.kwargs['temp_stream_dir'])

    def test_poisoning(self):
        self.process = Process(
            target=_pre_processor, kwargs=self.kwargs,
            name="TestProcess")
        poisoning(obj=self)

    def test_poisoning_while_waiting_on_output(self):
        self.process = Process(
            target=_pre_processor, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_while_waiting_on_output(
            obj=self, output_queue=self.kwargs["output_filename_queue"])

    def test_poisoning_from_input(self):
        self.process = Process(
            target=_pre_processor, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_from_input(
            obj=self, input_queue=self.kwargs['input_stream_queue'])


class TestGetDetectionStreamProcess(ProcessTests):
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
                 output_filename_queue=Queue(maxsize=1),
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
        poisoning_while_waiting_on_output(
            obj=self, output_queue=self.kwargs['output_filename_queue'])

    def test_poisoning_from_input(self):
        self.process = Process(
            target=_get_detection_stream, kwargs=self.kwargs,
            name="TestProcess")
        poisoning_from_input(
            obj=self, input_queue=self.kwargs['input_time_queue'])

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
        self.wait_for_close()

    def test_full_stream_operation(self):
        kwargs = copy.copy(self.kwargs)
        kwargs.update(dict(
            full_stream_dir=tempfile.mkdtemp(),
            pre_process=True,))
        self.process = Process(
            target=_get_detection_stream, kwargs=kwargs,
            name="TestProcess")
        self.process.start()

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
        self.wait_for_close()

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
        self.process = Process(
            target=_get_detection_stream, kwargs=kwargs, name="ProcessProcess")
        self.process.start()

        # Populate time queue
        kwargs['input_time_queue'].put(
            (UTCDateTime(2021, 1, 1),
             UTCDateTime(2021, 1, 1, 0, 10)))
        kwargs['input_time_queue'].put(None)

        time.sleep(self.wait_time)
        poison = kwargs['poison_queue'].get()
        self.assertIsInstance(poison, Poison)
        self.wait_for_close()


###############################################################################
#  STANDARD PROCESS DEATH TESTS
###############################################################################


def poisoning(obj: ProcessTests):
    obj.process.start()
    # Test death
    obj.kwargs['poison_queue'].put(Exception("TestException"))
    obj.wait_for_close()


def poisoning_while_waiting_on_output(obj: ProcessTests, output_queue: Queue):
    # Fill output queue
    output_queue.put("This is dog")
    obj.process.start()
    # Test death
    obj.kwargs['poison_queue'].put(Exception("TestException"))
    obj.wait_for_close()


def poisoning_from_input(obj: ProcessTests, input_queue: Queue):
    obj.process.start()
    # Test death
    input_queue.put(Poison(Exception("TestException")))
    obj.wait_for_close()


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
