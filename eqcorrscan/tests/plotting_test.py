"""
Plotting tests for EQcorrscan.
"""
import numpy as np
import unittest
import pytest
import datetime as dt
import os
import glob

from obspy import read, read_events, UTCDateTime, Catalog, Stream, Trace
from obspy.io.nordic.core import readwavename

from eqcorrscan.utils.plotting import (
    chunk_data, xcorr_plot, triple_plot, peaks_plot,
    cumulative_detections, threeD_gridplot, multi_event_singlechan,
    detection_multiplot, interev_mag, obspy_3d_plot, noise_plot,
    pretty_template_plot, plot_repicked, svd_plot, plot_synth_real,
    freq_mag, spec_trace, subspace_detector_plot, subspace_fc_plot,
    twoD_seismplot)
from eqcorrscan.utils.stacking import align_traces
from eqcorrscan.utils import findpeaks
from eqcorrscan.core.match_filter import normxcorr2
from eqcorrscan.core import template_gen, subspace


class SeimicityPlottingMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from obspy.clients.fdsn import Client
        client = Client("IRIS")
        starttime = UTCDateTime("2000-01-01")
        endtime = UTCDateTime("2020-05-16")
        cls.catalog = client.get_events(
            starttime=starttime, endtime=endtime, latitude=32.5,
            longitude=47.5, maxradius=0.7)

    @pytest.mark.mpl_image_compare
    def test_twoD_seismplot_depth_catalog(self):
        fig = twoD_seismplot(
            catalog=self.catalog, method='depth',
            show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_twoD_seismplot_time_catalog(self):
        fig = twoD_seismplot(
            catalog=self.catalog, method='time',
            show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_twoD_seismplot_sequence_catalog(self):
        fig = twoD_seismplot(
            catalog=self.catalog, method='sequence',
            show=False, return_figure=True)
        return fig


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
            channel="EZ", show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_multi_event_singlechan_realign(self):
        _, _, fig = multi_event_singlechan(
            streams=self.streams, catalog=self.catalog, station="GCSZ",
            channel="EZ", show=False, realign=True, freqmin=2, freqmax=20,
            PWS=True, return_figure=True)
        return fig


class EventPlottingMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'test_data', 'REA',
            'TEST_', '01-0411-15L.S201309')
        test_wavefile = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'test_data', 'WAV',
            'TEST_', '2013-09-01-0410-35.DFDPC_024_00')
        cls.event = read_events(test_file)[0]
        cls.st = read(test_wavefile)
        cls.st.detrend().filter("bandpass", freqmin=2, freqmax=20, corners=4)
        for tr in cls.st:
            tr = tr.trim(tr.stats.starttime + 30, tr.stats.endtime - 30)
            # Hack around seisan 2-letter channel naming
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        cls.template = template_gen._template_gen(cls.event.picks, cls.st, 2)

    @pytest.mark.mpl_image_compare
    def test_detection_multiplot(self):
        times = [min([pk.time - 0.05 for pk in self.event.picks])]
        times.append(times[0] + 10)
        fig = detection_multiplot(
            stream=self.st, template=self.template, times=times,
            show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_pretty_template_plot(self):
        fig = pretty_template_plot(
            self.template, background=self.st, event=self.event,
            show=False, return_figure=True, title="test template")
        return fig

    @pytest.mark.mpl_image_compare
    def test_pretty_template_plot_basic(self):
        fig = pretty_template_plot(
            self.template, show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_pretty_template_plot_sort(self):
        fig = pretty_template_plot(
            self.template, background=self.st, event=self.event,
            show=False, return_figure=True, title="sorted test template")
        return fig

    @pytest.mark.mpl_image_compare
    def test_pretty_template_plot_sort_by_picktime(self):
        fig = pretty_template_plot(
            self.template, background=self.st, event=self.event,
            sort_by="pick_time", show=False, return_figure=True,
            title="sorted test template")
        return fig

    @pytest.mark.mpl_image_compare
    def test_plot_repicked(self):
        _picks = [pick.copy() for pick in self.event.picks]
        picked_channels = []
        picks = []
        for pick in _picks:
            if pick.waveform_id not in picked_channels:
                pick.time += 3
                picks.append(pick)
                picked_channels.append(pick.waveform_id)
        fig = plot_repicked(
            template=self.template, picks=picks, det_stream=self.st,
            show=False, return_figure=True, title="test_detection")
        return fig


class StreamPlottingMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.st = read()

    @pytest.mark.mpl_image_compare
    def test_xcorr_plot(self):
        st = self.st.copy().detrend('simple').filter(
            'bandpass', freqmin=2, freqmax=15)
        shifts, ccs = align_traces([st[0], st[1]], 40)
        shift = shifts[1] * st[1].stats.sampling_rate
        cc = ccs[1]
        fig = xcorr_plot(template=st[1].data, image=st[0].data, shift=shift,
                         cc=cc, show=False, return_figure=True)
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
            show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_noise_plot(self):
        signal = self.st.copy().trim(
            self.st[0].stats.starttime + 10,
            self.st[0].stats.starttime + 15)
        noise = self.st.copy().trim(self.st[0].stats.starttime + 20)
        fig = noise_plot(
            signal=signal, noise=noise, show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_noise_plot_normalized(self):
        signal = self.st.copy().trim(
            self.st[0].stats.starttime + 10,
            self.st[0].stats.starttime + 15)
        noise = self.st.copy().trim(self.st[0].stats.starttime + 20)
        fig = noise_plot(
            signal=signal, noise=noise, normalise=True, show=False,
            return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_triple_plot(self):
        template = self.st[0].copy().trim(self.st[0].stats.starttime + 8,
                                          self.st[0].stats.starttime + 12)
        tr = self.st[0].copy()
        ccc = normxcorr2(template=template.data, image=tr.data)
        tr.data = tr.data[0:len(ccc[0])]
        fig = triple_plot(cccsum=ccc[0], cccsum_hist=ccc[0], trace=tr,
                          threshold=0.8, show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_plot_synth_real(self):
        from eqcorrscan.utils.synth_seis import seis_sim

        synth = Stream(Trace(seis_sim(sp=100, flength=200)))
        synth[0].stats.station = 'RJOB'
        synth[0].stats.channel = 'EHZ'
        synth[0].stats.sampling_rate = 100
        synth = synth.filter('bandpass', freqmin=2, freqmax=8)
        real = self.st.select(
           station='RJOB', channel='EHZ').detrend('simple').filter(
              'bandpass', freqmin=2, freqmax=8)
        real.trim(starttime=real[0].stats.starttime + 4.9,
                  endtime=real[0].stats.starttime + 6.9)
        fig = plot_synth_real(real_template=real, synthetic=synth,
                              size=(7, 4), show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_spec_trace(self):
        fig = spec_trace(
            traces=self.st, trc="white", show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_spec_trace_list(self):
        traces = [tr for tr in self.st]
        fig = spec_trace(
            traces=traces, trc="white", show=False, return_figure=True)
        return fig

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
                         samp_rate=10, peaks=peaks, show=False,
                         return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_cumulative_detections(self):
        dates = []
        for i in range(3):
            dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
                          for n in np.random.randn(100)])
        fig = cumulative_detections(
            dates, ['a', 'b', 'c'], show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_cumulative_detections_rate(self):
        dates = []
        for i in range(3):
            dates.append([dt.datetime(2012, 3, 26) + dt.timedelta(n)
                          for n in np.random.randn(100)])
        fig = cumulative_detections(
            dates, ['a', 'b', 'c'], show=False, plot_grouped=True, rate=True,
            return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_threeD_gridplot(self):
        nodes = [(-43.5, 170.4, 4), (-43.3, 170.8, 12), (-43.4, 170.3, 8)]
        fig = threeD_gridplot(nodes=nodes, show=False, return_figure=True)
        return fig


class MultiStreamPlottingMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        wavefiles = glob.glob(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'test_data', 'WAV',
            'TEST_', '2013-*'))
        streams = [read(w) for w in wavefiles[1:10]]
        cls.stream_list = []
        for st in streams:
            tr = st.select(station='GCSZ', channel='EHZ')
            tr = tr.detrend('simple').resample(100).filter(
                'bandpass', freqmin=2, freqmax=8)
            cls.stream_list.append(tr)
        cls.detector = subspace.Detector().read(os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            '..', 'tests', 'test_data', 'subspace', 'stat_test_detector.h5'))

    @pytest.mark.mpl_image_compare
    def test_svd_plot(self):
        from eqcorrscan.utils.clustering import svd, svd_to_stream
        uvec, sval, svec, stachans = svd(stream_list=self.stream_list)
        svstreams = svd_to_stream(uvectors=uvec, stachans=stachans, k=3,
                                  sampling_rate=100)
        fig = svd_plot(svstreams=svstreams, svalues=sval, stachans=stachans,
                       show=False, return_figure=True)
        return fig[0]

    @pytest.mark.mpl_image_compare
    def test_subspace_detector_plot(self):
        fig = subspace_detector_plot(
            detector=self.detector, stachans='all', size=(10, 7), show=False,
            return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_subspace_fc_plot(self):
        fig = subspace_fc_plot(
            detector=self.detector, stachans='all', size=(10, 7), show=False,
            return_figure=True)
        return fig


class NetworkPlottingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from obspy.clients.fdsn import Client
        client = Client('IRIS')
        t1 = UTCDateTime('2012-03-26T00:00:00')
        t2 = t1 + (3 * 86400)
        cls.catalog = client.get_events(
            starttime=t1, endtime=t2, latitude=0,
            longitude=170, maxradius=40)
        cls.inventory = client.get_stations(
            starttime=t1, endtime=t2, latitude=0, longitude=170,
            maxradius=40)

    @pytest.mark.mpl_image_compare
    def test_interev_mag(self):
        magnitudes = [event.preferred_magnitude().mag
                      for event in self.catalog]
        times = [event.preferred_origin().time for event in self.catalog]
        fig = interev_mag(times, magnitudes, show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_obspy_3d_plot(self):
        fig = obspy_3d_plot(
            inventory=self.inventory, catalog=self.catalog, show=False,
            return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_freq_mag(self):
        magnitudes = [event.preferred_magnitude().mag
                      for event in self.catalog]
        fig = freq_mag(
            magnitudes, completeness=3.5, max_mag=7, show=False,
            return_figure=True)
        return fig


if __name__ == "__main__":
    unittest.main()
