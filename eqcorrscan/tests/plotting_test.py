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
import os

from obspy import read, read_events

from eqcorrscan.utils.plotting import (
    _check_save_args, chunk_data, xcorr_plot, triple_plot, peaks_plot,
    cumulative_detections, threeD_gridplot, threeD_seismplot,
    multi_event_singlechan, multi_trace_plot, detection_multiplot,
    interev_mag, obspy_3d_plot, noise_plot, pretty_template_plot, plot_repicked,
    NR_plot, SVD_plot, plot_synth_real, freq_mag, spec_trace,
    subspace_detector_plot, subspace_fc_plot, _match_filter_plot,
    _plotting_decimation)
from eqcorrscan.utils.stacking import align_traces
from eqcorrscan.core.match_filter import normxcorr2


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
                         cc=cc)
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
            template=st[1].data[40:-40], image=st[0].data, cc_vec=cc_vec)
        return fig


if __name__ == "__main__":
    unittest.main()

