"""
Functions for testing the utils.pre_processing functions.
"""
import unittest
import os
import numpy as np
import glob

from random import shuffle
from copy import deepcopy
from obspy import read, Trace, UTCDateTime, Stream

from eqcorrscan.utils.pre_processing import (
    multi_process, _check_daylong, _prep_data_for_correlation,
    _multi_detrend, _multi_resample, _multi_filter,
)


class TestMultiThreadMethods(unittest.TestCase):
    """ Compare internal methods to obspy results. """
    @classmethod
    def setUpClass(cls):
        cls.real_st = read(os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data',
            'day_vols', 'Y2012', 'R086.01', '*'))
        random_st = Stream([Trace(np.random.randn(86401)) for _ in range(9)])
        for tr in random_st:
            tr.stats.sampling_rate = 100
        cls.random_st = random_st
        short_st = read(os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data',
            "WAV", "TEST_", "2013-09-11-2208-45.DFDPC_030_00"))
        cls.short_st = Stream(
            [tr for tr in short_st if tr.stats.sampling_rate == 200.0])
        # At the moment processing info is not retained in EQcorrscan
        cls.headers_to_compare = {
            "network", "station", "location", "channel", "starttime",
            "endtime", "sampling_rate", "delta", "npts", "calib"}

    def test_resample(self):
        samp_rates = [v / 10 for v in range(1, 10)]
        for st, name in zip([self.real_st, self.random_st, self.short_st],
                            ["real", "random", "short"]):
            print(f"Running for {name}")
            for samp_frac in samp_rates:
                print(f"Checking for {samp_frac} fractional sampling_rate")
                samp_rate = st[0].stats.sampling_rate * samp_frac
                acc_resample = _multi_resample(st.copy(), samp_rate)
                obspy_resample = st.copy().resample(samp_rate)
                # Order should not be changed so we can loop
                for acc_tr, obspy_tr in zip(acc_resample, obspy_resample):
                    for head in self.headers_to_compare:
                        assert acc_tr.stats[head] == obspy_tr.stats[head]
                    assert np.allclose(acc_tr.data, obspy_tr.data)

    def test_detrend(self):
        for st in [self.real_st, self.random_st, self.short_st]:
            acc_detrend = _multi_detrend(st.copy())
            obspy_detrend = st.copy().detrend()
            # Order should not be changed so we can loop
            for acc_tr, obspy_tr in zip(acc_detrend, obspy_detrend):
                for head in self.headers_to_compare:
                    assert acc_tr.stats[head] == obspy_tr.stats[head]
                assert np.allclose(acc_tr.data, obspy_tr.data)

    def test_bandpass(self):
        lows = [v / 20 for v in range(1, 9)]
        highs = [v / 20 for v in range(2, 10)]
        for st in [self.real_st, self.random_st, self.short_st]:
            for low, high in zip(lows, highs):
                lowcut = st[0].stats.sampling_rate * low
                highcut = st[0].stats.sampling_rate * high
                acc_filter = _multi_filter(
                    st.copy(), highcut=highcut, lowcut=lowcut, filt_order=4)
                obspy_filter = st.copy().filter(
                    "bandpass", freqmin=lowcut, freqmax=highcut, corners=4,
                    zerophase=True)
                # Order should not be changed so we can loop
                for acc_tr, obspy_tr in zip(acc_filter, obspy_filter):
                    for head in self.headers_to_compare:
                        assert acc_tr.stats[head] == obspy_tr.stats[head]
                    assert np.allclose(acc_tr.data, obspy_tr.data)

    def test_lowpass(self):
        highs = [v / 20 for v in range(2, 10)]
        for st in [self.real_st, self.random_st, self.short_st]:
            for high in highs:
                highcut = st[0].stats.sampling_rate * high
                acc_filter = _multi_filter(
                    st.copy(), highcut=highcut, lowcut=None, filt_order=4)
                obspy_filter = st.copy().filter(
                    "lowpass", freq=highcut, corners=4,
                    zerophase=True)
                # Order should not be changed so we can loop
                for acc_tr, obspy_tr in zip(acc_filter, obspy_filter):
                    for head in self.headers_to_compare:
                        assert acc_tr.stats[head] == obspy_tr.stats[head]
                    assert np.allclose(acc_tr.data, obspy_tr.data)

    def test_highpass(self):
        lows = [v / 20 for v in range(1, 9)]
        for st in [self.real_st, self.random_st, self.short_st]:
            for low in lows:
                lowcut = st[0].stats.sampling_rate * low
                acc_filter = _multi_filter(
                    st.copy(), lowcut=lowcut, highcut=None, filt_order=4)
                obspy_filter = st.copy().filter(
                    "highpass", freq=lowcut, corners=4,
                    zerophase=True)
                # Order should not be changed so we can loop
                for acc_tr, obspy_tr in zip(acc_filter, obspy_filter):
                    for head in self.headers_to_compare:
                        assert acc_tr.stats[head] == obspy_tr.stats[head]
                    assert np.allclose(acc_tr.data, obspy_tr.data)


class TestPreProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data',
            'day_vols', 'Y2012', 'R086.01', '*')
        cls.st = read(testing_path)
        cls.day_start = cls.st[0].stats.starttime.date
        cls.short_stream = cls.st.copy().trim(cls.st[0].stats.starttime,
                                              cls.st[0].stats.starttime + 3600)
        cls.instart = cls.short_stream[0].stats.starttime
        cls.inend = cls.short_stream[0].stats.endtime
        cls.nchans = len(cls.short_stream)
        cls.gap_starttime = cls.st[0].stats.starttime + 1800
        cls.gap_endtime = cls.st[0].stats.starttime + 1900
        tr1 = cls.st[0].copy().trim(
            cls.st[0].stats.starttime, cls.gap_starttime)
        tr2 = cls.st[0].copy().trim(
            cls.gap_endtime, cls.st[0].stats.starttime + 3600)
        cls.gappy_trace = tr1 + tr2

    def test_daylong_checks(self):
        """Test that the data are day-long."""
        self.assertTrue(_check_daylong(self.st[0].data))
        not_daylong = self.st[0].copy().trim(self.st[0].stats.starttime,
                                             self.st[0].stats.starttime + 3600)
        not_daylong.data = np.append(
            not_daylong.data, np.zeros(
                3602 * int(self.st[0].stats.sampling_rate)))
        self.assertFalse(_check_daylong(not_daylong.data))

    def test_shortproc(self):
        """Test the short-proc processing method."""
        processed = multi_process(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=False, num_cores=False,
            starttime=None, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_filter_error(self):
        """Check that we don't allow filtering above the nyquist."""
        with self.assertRaises(IOError):
            multi_process(
                self.short_stream.copy(), lowcut=0.1, highcut=0.6,
                filt_order=4, samp_rate=1, parallel=False,
                num_cores=False, starttime=None, endtime=None)

    def test_shortproc_set_start(self):
        """Check that shortproc trims properly."""
        processed = multi_process(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=False, num_cores=False,
            starttime=self.short_stream[0].stats.starttime + 2, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart + 2, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_shortproc_set_end(self):
        """Check that shortproc trims properly."""
        processed = multi_process(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=False, num_cores=False,
            starttime=None, endtime=self.short_stream[0].stats.endtime - 2)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend - 2, tr.stats.endtime)

    def test_shortproc_set_start_and_end(self):
        """Check that shortproc trims properly."""
        processed = multi_process(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=False, num_cores=False,
            starttime=self.short_stream[0].stats.starttime + 2,
            endtime=self.short_stream[0].stats.endtime - 2)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart + 2, tr.stats.starttime)
            self.assertEqual(self.inend - 3, tr.stats.endtime)

    def test_trace_as_argument(self):
        """
        Check that we can cope with a trace, and that a trace is returned.
        """
        processed = multi_process(
            self.short_stream.copy()[0], lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=False, num_cores=False,
            starttime=None, endtime=None)
        self.assertTrue(isinstance(processed, Trace))
        self.assertEqual(self.instart, processed.stats.starttime)
        self.assertEqual(self.inend, processed.stats.endtime)

    def test_parallel(self):
        """Test the parallel implementation."""
        processed = multi_process(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=True, num_cores=2,
            starttime=None, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_parallel_core_unset(self):
        """Test the parallel implementation without num_cores set."""
        processed = multi_process(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, parallel=True, num_cores=False,
            starttime=None, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_dayproc(self):
        """Test a straight-forward day processing implementation."""
        processed = multi_process(
            st=self.st.copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, parallel=True,
            num_cores=2, daylong=True)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(UTCDateTime(self.day_start), tr.stats.starttime)
            self.assertEqual(tr.stats.npts, 86400)

    def test_dayproc_trace(self):
        """
        Test a straight-forward day processing implementation with a Trace.
        """
        processed = multi_process(
            st=self.st[0].copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, parallel=True,
            num_cores=2, daylong=True)
        self.assertTrue(isinstance(processed, Trace))
        self.assertEqual(UTCDateTime(self.day_start),
                         processed.stats.starttime)
        self.assertEqual(processed.stats.npts, 86400)

    def test_dayproc_nyquist_error(self):
        """Test a failing day processing."""
        with self.assertRaises(IOError):
            multi_process(
                st=self.st.copy(), lowcut=0.1, highcut=0.6, filt_order=3,
                samp_rate=1, starttime=self.day_start,
                parallel=True, num_cores=2, daylong=True)

    def test_dayproc_serial(self):
        """Test the serial implementation of dayproc."""
        processed = multi_process(
            st=self.st.copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, parallel=False,
            num_cores=2, daylong=True)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(UTCDateTime(self.day_start), tr.stats.starttime)
            self.assertEqual(tr.stats.npts, 86400)

    def test_dayproc_parallel_cores_unset(self):
        """Test a straight-forward day processing implementation."""
        processed = multi_process(
            st=self.st.copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, parallel=True,
            num_cores=False, daylong=True)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(UTCDateTime(self.day_start), tr.stats.starttime)
            self.assertEqual(tr.stats.npts, 86400)

    def test_process(self):
        """Test a basic process implementation."""
        processed = multi_process(
            st=self.st[0].copy(), lowcut=0.1, highcut=0.4,
            filt_order=3, samp_rate=1, starttime=False,
            daylong=True, seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_process_datetime(self):
        """Test a basic process implementation."""
        processed = multi_process(
            st=self.st[0].copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, daylong=True,
            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_process_nyquist_fail(self):
        """Test a nyquist error is raised."""
        with self.assertRaises(IOError):
            multi_process(
                st=self.st[0].copy(), lowcut=0.1, highcut=0.6,
                filt_order=3, samp_rate=1, starttime=False, daylong=True,
                seisan_chan_names=True, ignore_length=False)

    def test_process_bad_data(self):
        """Check that we won't allow data that are mostly zeros."""
        not_daylong = self.st[0].copy().trim(self.st[0].stats.starttime,
                                             self.st[0].stats.starttime + 3600)
        not_daylong.data = np.append(
            not_daylong.data, np.zeros(
                3602 * int(self.st[0].stats.sampling_rate)))
        with self.assertRaises(ValueError):
            multi_process(
                st=not_daylong, lowcut=0.1, highcut=0.4,
                filt_order=3, samp_rate=1, starttime=False, daylong=True,
                seisan_chan_names=True, ignore_length=False)

    def test_short_data_fail(self):
        """Check that we don't allow too much missing data."""
        with self.assertRaises(NotImplementedError):
            multi_process(
                st=self.st[0].copy().trim(
                    endtime=self.st[0].stats.endtime - 18000),
                lowcut=0.1, highcut=0.4, filt_order=3, samp_rate=1,
                starttime=self.day_start, daylong=True,
                seisan_chan_names=True, ignore_length=False)

    def test_short_data_pass(self):
        """Check that we do allow missing data if ignore_length is True."""
        processed = multi_process(
            st=self.st[0].copy().trim(endtime=self.
                                      st[0].stats.endtime - 18000), lowcut=0.1,
            highcut=0.4, filt_order=3, samp_rate=1,
            starttime=self.day_start, daylong=True,
            seisan_chan_names=True, ignore_length=True)
        self.assertEqual(processed.stats.npts, 86400)

    def test_short_data_empty_return(self):
        """
        Check that we do not include data that is too short even if
        ignore_bad_data is True.
        """
        processed = multi_process(
            st=self.st[0].copy().trim(endtime=self.
                                      st[0].stats.endtime - 28000), lowcut=0.1,
            highcut=0.4, filt_order=3, samp_rate=1,
            starttime=self.day_start, daylong=True,
            seisan_chan_names=True, ignore_bad_data=True)
        self.assertEqual(processed.stats.npts, 0)

    def test_highcut_debug(self):
        """Test a basic process implementation with just a highcut"""
        processed = multi_process(
            st=self.st[0].copy(), lowcut=None, highcut=0.4,
            filt_order=3, samp_rate=1, starttime=False, daylong=True,
            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_lowcut_debug(self):
        """Test a basic process implementation with just a highcut"""
        processed = multi_process(
            st=self.st[0].copy(), lowcut=0.1, highcut=None,
            filt_order=3, samp_rate=1, starttime=False, daylong=True,
            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_masked_trace(self):
        """Test that processing a masked array works."""
        tr = self.gappy_trace.copy()
        processed = multi_process(
            st=tr, lowcut=0.1, highcut=0.4,
            filt_order=3, samp_rate=1, starttime=False,
            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 3601)
        self.assertFalse(isinstance(processed.data, np.ma.MaskedArray))
        self.assertTrue(np.all(
            processed.trim(self.gap_starttime, self.gap_endtime).data) == 0)

    def test_masked_trace_no_fill(self):
        """Test that processing a masked array without filling gaps works."""
        tr = self.gappy_trace.copy()
        processed = multi_process(
            st=tr, lowcut=0.1, highcut=0.4, filt_order=3, samp_rate=1,
            starttime=False, seisan_chan_names=True, ignore_length=False,
            fill_gaps=False)
        self.assertEqual(processed.stats.npts, 3601)
        self.assertTrue(isinstance(processed.data, np.ma.MaskedArray))

    def test_masked_array_resample(self):
        """Test that processing and resampling a masked array works."""
        tr = self.gappy_trace.copy()
        processed = multi_process(
            st=tr, lowcut=0.1, highcut=0.2, filt_order=3, samp_rate=0.5,
            starttime=False, seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 1800)
        self.assertTrue(np.all(
            processed.trim(self.gap_starttime, self.gap_endtime).data) == 0)

    def test_gap_overlength(self):
        """ Test that gappy data that are too long are processed correctly."""
        tr_after = self.gappy_trace.copy()
        tr_after.stats.starttime += 3600
        tr_before = self.gappy_trace.copy()
        tr_before.stats.starttime -= 3600
        tr = tr_before + self.gappy_trace + tr_after
        processed = multi_process(
            st=tr, lowcut=0.1, highcut=0.4, filt_order=3, samp_rate=1,
            starttime=self.gappy_trace.stats.starttime,
            endtime=self.gappy_trace.stats.starttime + 3600,
            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 3600)
        self.assertFalse(isinstance(processed.data, np.ma.MaskedArray))
        self.assertTrue(np.all(
            processed.trim(self.gap_starttime, self.gap_endtime).data) == 0)
        # Check that there is actually data there!
        self.assertEqual(np.isnan(processed.data).sum(), 0)


class TestDataPrep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data',
            'similar_events_processed')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file)
                       for stream_file in stream_files]
        cls.stream_list = stream_list

    def test_order_irrevelant(self):
        """ Early on the order of streams input was affecting output. """
        st1 = self.stream_list[0].copy()

        ordered_streams = sorted(
            self.stream_list, key=lambda s: s[0].stats.starttime)
        shuffled_streams = deepcopy(ordered_streams)
        shuffle(shuffled_streams)

        order_prepped_st1, order_prepped_streams = _prep_data_for_correlation(
            stream=st1, templates=ordered_streams, force_stream_epoch=False)
        shuff_prepped_st1, shuff_prepped_streams = _prep_data_for_correlation(
            stream=st1, templates=shuffled_streams, force_stream_epoch=False)

        self.assertEqual(order_prepped_st1, shuff_prepped_st1)
        shuff_prepped_streams = sorted(
            shuff_prepped_streams, key=lambda s: s[0].stats.starttime)

        for ordered, shuffled in zip(order_prepped_streams,
                                     shuff_prepped_streams):
            for tr_ordered, tr_shuffled in zip(ordered, shuffled):
                self.assertEqual(tr_ordered, tr_shuffled)

    def test_correct_padding(self):
        """check that nans are put where they are supposed to be """
        _, prepped_streams = _prep_data_for_correlation(
            stream=self.stream_list[0], templates=deepcopy(self.stream_list),
            force_stream_epoch=False)

        for original, prepped in zip(self.stream_list, prepped_streams):
            for tr in prepped:
                original_tr = original.select(id=tr.id)
                if len(original_tr) == 0:
                    self.assertTrue(np.all(np.isnan(tr.data)))
                else:
                    self.assertEqual(tr, original_tr[0])

    def test_duplicate_template_channels(self):
        """
        Check that duplicate template channels result in duplicate template
        channels but not stream channels.
        """
        templates = deepcopy(self.stream_list)
        templates.sort(key=lambda t: t[0].stats.starttime)
        stream = deepcopy(self.stream_list[0])
        additional_channel = templates[0][0]
        additional_channel.stats.starttime += 30
        templates[0] += additional_channel
        assert len(templates[0]) == 10
        continuous_data, templates = _prep_data_for_correlation(
            stream=stream, templates=templates, force_stream_epoch=True)
        for template in templates:
            assert len(template) == 10
        assert len(continuous_data) == 9

    def test_continuous_data_removal(self):
        """Check that data that should be removed are."""
        st = read()
        templates = []
        for i in range(20):
            template = Stream()
            for _tr in st[1:]:
                tr = _tr.copy()
                tr.stats.starttime += i * 100
                tr.data = np.random.randn(200)
                template += tr
            templates.append(template)
        continuous_data = st
        for tr in continuous_data:
            tr.data = np.random.randn(200 * 20)

        continuous_data, templates = _prep_data_for_correlation(
            stream=continuous_data, templates=templates,
            force_stream_epoch=True)
        self.assertEqual(len(continuous_data), 2)
        for template in templates:
            self.assertEqual(len(template), 2)

    def test_length_checking(self):
        """Check that lengths are fixed"""
        st = read()
        templates = []
        for i in range(20):
            template = Stream()
            for _tr in st:
                tr = _tr.copy()
                tr.stats.starttime += i * 100
                tr.data = np.random.randn(200)
                template += tr
            templates.append(template)
        continuous_data = st
        for i, tr in enumerate(continuous_data):
            tr.data = np.random.randn((200 * 20) + ((i + 1) * 5))
            tr.stats.starttime += (i + 1) * 5
        min_start = min([tr.stats.starttime for tr in continuous_data])
        max_end = max([tr.stats.endtime for tr in continuous_data])
        data_len = (
            (max_end - min_start) * continuous_data[0].stats.sampling_rate) + 1

        continuous_data_prepped, templates = _prep_data_for_correlation(
            stream=continuous_data.copy(), templates=templates,
            force_stream_epoch=True)

        for tr in continuous_data_prepped:
            self.assertEqual(tr.stats.npts, data_len)
            self.assertEqual(tr.stats.starttime, min_start)
            self.assertEqual(tr.stats.endtime, max_end)

        continuous_data_prepped, templates = _prep_data_for_correlation(
            stream=continuous_data, templates=templates,
            force_stream_epoch=False)

        for tr in continuous_data_prepped:
            self.assertEqual(
                tr.stats.npts, (200 * 20) + (len(continuous_data) * 5))


if __name__ == '__main__':
    unittest.main()
