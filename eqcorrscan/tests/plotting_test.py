"""
Plotting tests for EQcorrscan.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import unittest
import pytest
import datetime as dt
import os
import glob

from obspy import read, read_events, UTCDateTime, Catalog
from obspy.io.nordic.core import readwavename

from eqcorrscan.utils.plotting import (
    _check_save_args, chunk_data, xcorr_plot, triple_plot, peaks_plot,
    cumulative_detections, threeD_gridplot, threeD_seismplot,
    multi_event_singlechan, multi_trace_plot, detection_multiplot,
    interev_mag, obspy_3d_plot, noise_plot, pretty_template_plot, plot_repicked,
    NR_plot, SVD_plot, plot_synth_real, freq_mag, spec_trace,
    subspace_detector_plot, subspace_fc_plot, _match_filter_plot,
    _plotting_decimation)
from eqcorrscan.utils.stacking import align_traces
from eqcorrscan.utils import findpeaks
from eqcorrscan.core.match_filter import normxcorr2


class MultiStreamMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sfiles = glob.glob(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'test_data/REA/TEST_/*.S??????'))
        cls.catalog = Catalog()
        cls.streams = []
        for sfile in sfiles:
            cls.catalog += read_events(sfile)
            wavefile = readwavename(sfile)[0]
            stream_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'test_data/WAV/TEST_', wavefile)
            stream = read(stream_path)
            for tr in stream:
                tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
            cls.streams.append(stream)

    @pytest.mark.mpl_image_compare
    def test_multi_event_singlechan(self):
        _, _, fig = multi_event_singlechan(
            streams=self.streams, catalog=self.catalog, station="GCSZ",
            channel="EZ", show=False)
        return fig

    @pytest.mark.mpl_image_compare
    def test_multi_event_singlechan_realign(self):
        _, _, fig = multi_event_singlechan(
            streams=self.streams, catalog=self.catalog, station="GCSZ",
            channel="EZ", show=False, realign=True, freqmin=2, freqmax=20,
            PWS=True)
        return fig


class StreamPlottingMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.st = read()

    def test_chunk_data(self):
        desired_rate = 10
        trout = chunk_data(self.st[0], samp_rate=desired_rate, state="Mean")
        assert trout.stats.sampling_rate == desired_rate
        assert np.allclose(trout.data.mean(), self.st[0].data.mean(),
                           atol=.00000000001)
        trout = chunk_data(self.st[0], samp_rate=desired_rate, state="Min")
        assert trout.data.mean() < self.st[0].data.mean()
        assert trout.stats.sampling_rate == desired_rate
        trout = chunk_data(self.st[0], samp_rate=desired_rate, state="Max")
        assert trout.stats.sampling_rate == desired_rate
        assert trout.data.mean() > self.st[0].data.mean()
        trout = chunk_data(self.st[0], samp_rate=desired_rate, state="Maxabs")
        assert trout.stats.sampling_rate == desired_rate

    @pytest.mark.mpl_image_compare
    def test_xcorr_plot(self):
        st = self.st.copy().detrend('simple').filter(
            'bandpass', freqmin=2, freqmax=15)
        shifts, ccs = align_traces([st[0], st[1]], 40)
        shift = shifts[1] * st[1].stats.sampling_rate
        cc = ccs[1]
        fig = xcorr_plot(template=st[1].data, image=st[0].data, shift=shift,
                         cc=cc, show=False)
        return fig

    @pytest.mark.mpl_image_compare
    def test_xcorr_plot_cc_vec(self):
        st = self.st.copy().detrend('simple').filter(
            'bandpass', freqmin=2, freqmax=15)
        cc_vec = normxcorr2(
            template=st[0].data.astype(np.float32)[40:-40],
            image=st[1].data.astype(np.float32))
        cc_vec = cc_vec[0]
        fig = xcorr_plot(
            template=st[1].data[40:-40], image=st[0].data, cc_vec=cc_vec,
            show=False)
        return fig

    @pytest.mark.mpl_image_compare
    def test_triple_plot(self):
        template = self.st[0].copy().trim(self.st[0].stats.starttime + 8,
                                          self.st[0].stats.starttime + 12)
        tr = self.st[0]
        ccc = normxcorr2(template=template.data, image=tr.data)
        tr.data = tr.data[0:len(ccc[0])]
        fig = triple_plot(cccsum=ccc[0], cccsum_hist=ccc[0], trace=tr,
                          threshold=0.8, show=False)
        return fig


class DataPlottingMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = np.random.randn(200)

    @pytest.mark.mpl_image_compare
    def test_peaks_plot(self):
        data = self.data.copy()
        data[30] = 100
        data[60] = 40
        threshold = 10
        peaks = findpeaks.find_peaks2_short(data, threshold, 3)
        fig = peaks_plot(data=data, starttime=UTCDateTime("2008001"),
                         samp_rate=10, peaks=peaks, show=False)
        return fig

    @pytest.mark.mpl_image_compare
    def test_cumulative_detections(self):
        dates = []
        for i in range(3):
            dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
                          for n in np.random.randn(100)])
        fig = cumulative_detections(dates, ['a', 'b', 'c'], show=False)
        return fig

    @pytest.mark.mpl_image_compare
    def test_cumulative_detections_rate(self):
        dates = []
        for i in range(3):
            dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
                          for n in np.random.randn(100)])
        fig = cumulative_detections(dates, ['a', 'b', 'c'], show=False,
                                    plot_grouped=True, rate=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_threeD_gridplot(self):
        nodes = [(-43.5, 170.4, 4), (-43.3, 170.8, 12), (-43.4, 170.3, 8)]
        fig = threeD_gridplot(nodes=nodes, show=False)
        return fig


if __name__ == "__main__":
    unittest.main()

