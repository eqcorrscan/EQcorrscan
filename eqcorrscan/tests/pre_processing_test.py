"""
Functions for testing the utils.pre_processing functions.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import os
import numpy as np

from obspy import read, Trace, UTCDateTime

from eqcorrscan.utils.pre_processing import process, dayproc, shortproc
from eqcorrscan.utils.pre_processing import _check_daylong


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
        self.assertTrue(_check_daylong(self.st[0]))
        not_daylong = self.st[0].copy().trim(self.st[0].stats.starttime,
                                             self.st[0].stats.starttime + 3600)
        not_daylong.data = np.append(
            not_daylong.data, np.zeros(
                3602 * int(self.st[0].stats.sampling_rate)))
        self.assertFalse(_check_daylong(not_daylong))

    def test_shortproc(self):
        """Test the short-proc processing method."""
        processed = shortproc(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=False, num_cores=False,
            starttime=None, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_filter_error(self):
        """Check that we don't allow filtering above the nyquist."""
        with self.assertRaises(IOError):
            shortproc(self.short_stream.copy(), lowcut=0.1, highcut=0.6,
                      filt_order=4, samp_rate=1, debug=0, parallel=False,
                      num_cores=False, starttime=None, endtime=None)

    def test_shortproc_set_start(self):
        """Check that shortproc trims properly."""
        processed = shortproc(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=False, num_cores=False,
            starttime=self.short_stream[0].stats.starttime + 2, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart + 2, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_shortproc_set_end(self):
        """Check that shortproc trims properly."""
        processed = shortproc(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=False, num_cores=False,
            starttime=None, endtime=self.short_stream[0].stats.endtime - 2)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend - 2, tr.stats.endtime)

    def test_shortproc_set_start_and_end(self):
        """Check that shortproc trims properly."""
        processed = shortproc(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=False, num_cores=False,
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
        processed = shortproc(
            self.short_stream.copy()[0], lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=False, num_cores=False,
            starttime=None, endtime=None)
        self.assertTrue(isinstance(processed, Trace))
        self.assertEqual(self.instart, processed.stats.starttime)
        self.assertEqual(self.inend, processed.stats.endtime)

    def test_parallel(self):
        """Test the parallel implementation."""
        processed = shortproc(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=True, num_cores=2,
            starttime=None, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)

    def test_parallel_core_unset(self):
        """Test the parallel implementation without num_cores set."""
        processed = shortproc(
            self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
            samp_rate=1, debug=0, parallel=True, num_cores=False,
            starttime=None, endtime=None)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(self.instart, tr.stats.starttime)
            self.assertEqual(self.inend, tr.stats.endtime)
    # This tries to plot, can't run.
    # def test_high_debug(self):
    #     """Test the debug=5 level"""
    #     processed = shortproc(
    #         self.short_stream.copy(), lowcut=0.1, highcut=0.4, filt_order=4,
    #         samp_rate=1, debug=5, parallel=True, num_cores=False,
    #         starttime=None, endtime=None)
    #     self.assertEqual(len(processed), self.nchans)
    #     for tr in processed:
    #         self.assertEqual(self.instart, tr.stats.starttime)
    #         self.assertEqual(self.inend, tr.stats.endtime)

    def test_dayproc(self):
        """Test a straight-forward day processing implementation."""
        processed = dayproc(
            st=self.st.copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, debug=0, parallel=True,
            num_cores=2)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(UTCDateTime(self.day_start), tr.stats.starttime)
            self.assertEqual(tr.stats.npts, 86400)

    def test_dayproc_trace(self):
        """
        Test a straight-forward day processing implementation with a Trace.
        """
        processed = dayproc(
            st=self.st[0].copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, debug=0, parallel=True,
            num_cores=2)
        self.assertTrue(isinstance(processed, Trace))
        self.assertEqual(UTCDateTime(self.day_start),
                         processed.stats.starttime)
        self.assertEqual(processed.stats.npts, 86400)

    def test_dayproc_nyquist_error(self):
        """Test a failing day processing."""
        with self.assertRaises(IOError):
            dayproc(st=self.st.copy(), lowcut=0.1, highcut=0.6, filt_order=3,
                    samp_rate=1, starttime=self.day_start, debug=0,
                    parallel=True, num_cores=2)

    def test_dayproc_serial(self):
        """Test the serial implementation of dayproc."""
        processed = dayproc(st=self.st.copy(), lowcut=0.1, highcut=0.4,
                            filt_order=3, samp_rate=1,
                            starttime=self.day_start, debug=0, parallel=False,
                            num_cores=2)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(UTCDateTime(self.day_start), tr.stats.starttime)
            self.assertEqual(tr.stats.npts, 86400)

    def test_dayproc_parallel_cores_unset(self):
        """Test a straight-forward day processing implementation."""
        processed = dayproc(
            st=self.st.copy(), lowcut=0.1, highcut=0.4, filt_order=3,
            samp_rate=1, starttime=self.day_start, debug=0, parallel=True,
            num_cores=False)
        self.assertEqual(len(processed), self.nchans)
        for tr in processed:
            self.assertEqual(UTCDateTime(self.day_start), tr.stats.starttime)
            self.assertEqual(tr.stats.npts, 86400)

    def test_process(self):
        """Test a basic process implementation."""
        processed = process(tr=self.st[0].copy(), lowcut=0.1, highcut=0.4,
                            filt_order=3, samp_rate=1, debug=0,
                            starttime=False, clip=False, length=86400,
                            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_process_with_debug(self):
        """Test the debug option of 2 for process."""
        processed = process(tr=self.st[0].copy(), lowcut=0.1, highcut=0.4,
                            filt_order=3, samp_rate=1, debug=2,
                            starttime=False, clip=False, length=86400,
                            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_process_datetime(self):
        """Test a basic process implementation."""
        processed = process(tr=self.st[0].copy(), lowcut=0.1, highcut=0.4,
                            filt_order=3, samp_rate=1, debug=0,
                            starttime=self.day_start, clip=False, length=86400,
                            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_process_nyquist_fail(self):
        """Test a nyquist error is raised."""
        with self.assertRaises(IOError):
            process(tr=self.st[0].copy(), lowcut=0.1, highcut=0.6,
                    filt_order=3, samp_rate=1, debug=0,
                    starttime=False, clip=False, length=86400,
                    seisan_chan_names=True, ignore_length=False)

    def test_process_bad_data(self):
        """Check that we won't allow data that are mostly zeros."""
        not_daylong = self.st[0].copy().trim(self.st[0].stats.starttime,
                                             self.st[0].stats.starttime + 3600)
        not_daylong.data = np.append(
            not_daylong.data, np.zeros(
                3602 * int(self.st[0].stats.sampling_rate)))
        with self.assertRaises(ValueError):
            process(tr=not_daylong, lowcut=0.1, highcut=0.4,
                    filt_order=3, samp_rate=1, debug=0,
                    starttime=False, clip=False, length=86400,
                    seisan_chan_names=True, ignore_length=False)

    def test_short_data_fail(self):
        """Check that we don't allow too much missing data."""
        with self.assertRaises(NotImplementedError):
            process(tr=self.st[0].copy().
                    trim(endtime=self.st[0].stats.endtime - 18000), lowcut=0.1,
                    highcut=0.4, filt_order=3, samp_rate=1, debug=0,
                    starttime=self.day_start, clip=True, length=86400,
                    seisan_chan_names=True, ignore_length=False)

    def test_short_data_pass(self):
        """Check that we do allow missing data if ignore_length is True."""
        processed = process(
            tr=self.st[0].copy().trim(endtime=self.
                                      st[0].stats.endtime - 18000), lowcut=0.1,
            highcut=0.4, filt_order=3, samp_rate=1, debug=0,
            starttime=self.day_start, clip=True, length=86400,
            seisan_chan_names=True, ignore_length=True)
        self.assertEqual(processed.stats.npts, 86400)

    def test_highcut_debug(self):
        """Test a basic process implementation with just a highcut"""
        processed = process(tr=self.st[0].copy(), lowcut=None, highcut=0.4,
                            filt_order=3, samp_rate=1, debug=2,
                            starttime=False, clip=False, length=86400,
                            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_lowcut_debug(self):
        """Test a basic process implementation with just a highcut"""
        processed = process(tr=self.st[0].copy(), lowcut=0.1, highcut=None,
                            filt_order=3, samp_rate=1, debug=2,
                            starttime=False, clip=False, length=86400,
                            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 86400)

    def test_masked_trace(self):
        """Test that processing a masked array works."""
        tr = self.gappy_trace
        processed = process(tr=tr, lowcut=0.1, highcut=0.4,
                            filt_order=3, samp_rate=1, debug=0,
                            starttime=False, clip=False, length=3600,
                            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 3601)
        self.assertFalse(isinstance(processed.data, np.ma.MaskedArray))
        self.assertTrue(np.all(
            processed.trim(self.gap_starttime, self.gap_endtime).data) == 0)

    def test_masked_trace_no_fill(self):
        """Test that processing a masked array without filling gaps works."""
        tr = self.gappy_trace
        processed = process(tr=tr, lowcut=0.1, highcut=0.4,
                            filt_order=3, samp_rate=1, debug=0,
                            starttime=False, clip=False, length=3600,
                            seisan_chan_names=True, ignore_length=False,
                            fill_gaps=False)
        self.assertEqual(processed.stats.npts, 3601)
        self.assertTrue(isinstance(processed.data, np.ma.MaskedArray))

    def test_masked_array_resample(self):
        """Test that processing and resampling a masked array works."""
        tr = self.gappy_trace
        processed = process(tr=tr, lowcut=0.1, highcut=0.2,
                            filt_order=3, samp_rate=0.5, debug=0,
                            starttime=False, clip=False, length=3600,
                            seisan_chan_names=True, ignore_length=False)
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
        processed = process(
            tr=tr, lowcut=0.1, highcut=0.4, filt_order=3, samp_rate=1, debug=0,
            starttime=self.gappy_trace.stats.starttime, clip=True, length=3600,
            seisan_chan_names=True, ignore_length=False)
        self.assertEqual(processed.stats.npts, 3600)
        self.assertFalse(isinstance(processed.data, np.ma.MaskedArray))
        self.assertTrue(np.all(
            processed.trim(self.gap_starttime, self.gap_endtime).data) == 0)



if __name__ == '__main__':
    unittest.main()
