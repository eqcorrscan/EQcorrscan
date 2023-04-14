"""
A series of test functions for the core functions in EQcorrscan.
"""
import copy
import os
import shutil
import glob
import unittest
import pytest
import logging

import numpy as np
from obspy import read, UTCDateTime, read_events, Catalog, Stream, Trace
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException
from obspy.clients.earthworm import Client as EWClient
from obspy.core.event import Pick, Event
from obspy.core.util.base import NamedTemporaryFile

from eqcorrscan.core.match_filter import (
    normxcorr2, Detection, read_detections, get_catalog,
    write_catalog, extract_from_stream, Tribe, Template, Party, Family,
    read_party, read_tribe, _spike_test)
from eqcorrscan.core.match_filter.matched_filter import (
    match_filter, MatchFilterError)
from eqcorrscan.core.match_filter.helpers import get_waveform_client
from eqcorrscan.core.match_filter.template import (
    quick_group_templates, group_templates_by_seedid)

from eqcorrscan.utils import pre_processing, catalog_utils
from eqcorrscan.utils.correlate import fftw_normxcorr, numpy_normxcorr
from eqcorrscan.utils.catalog_utils import filter_picks


Logger = logging.getLogger(__name__)


class TestHelpers(unittest.TestCase):
    def test_monkey_patching(self):
        """ Test that monkey patching a client works. """
        client = EWClient("pubavo1.wr.usgs.gov", 16022)
        self.assertFalse(hasattr(client, "get_waveforms_bulk"))
        client = get_waveform_client(client)
        self.assertTrue(hasattr(client, "get_waveforms_bulk"))
        # TODO: This should test that the method actually works as expected.


class TestCoreMethods(unittest.TestCase):
    """
    Tests for internal _template_loop and normxcorr2 functions.
    """
    # changed to overflowerror in 0.5.0
    def test_detection_assertion(self):
        with self.assertRaises(OverflowError):
            Detection(
                template_name='a', detect_time=UTCDateTime(), threshold=1.2,
                threshold_input=8.0, threshold_type="MAD", typeofdet="corr",
                no_chans=3, detect_val=20)

    def test_perfect_normxcorr2(self):
        """
        Simple test of normxcorr2 to ensure data are detected
        """
        template = np.random.randn(100).astype(np.float32)
        image = np.zeros(1000).astype(np.float32)
        image[200] = 1.0
        image = np.convolve(template, image)
        ccc = normxcorr2(template, image).astype(np.float16)
        self.assertEqual(ccc.max(), 1.0)

    def test_fail_normxcorr2(self):
        """
        Ensure if template is nan then return is nan
        """
        template = np.array([np.nan] * 100)
        image = np.zeros(1000)
        image[200] = 1.0
        image = np.convolve(template, image)
        ccc = normxcorr2(template, image)
        self.assertTrue(np.all(ccc == 0.0))

    def test_normal_normxcorr2(self):
        """
        Check that if match is not perfect correlation max isn't unity
        """
        template = np.random.randn(100) * 10.0
        image = np.zeros(1000)
        image[200] = 1.0
        image = np.convolve(template, image)
        image += np.random.randn(1099)  # Add random noise
        ccc = normxcorr2(template, image)
        self.assertNotEqual(ccc.max(), 1.0)

    def test_set_normxcorr2(self):
        """
        Check that correlations output are the same irrespective of version.
        """
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        template = read(os.path.join(testing_path, 'test_template.ms'))
        template = template[0].data.astype(np.float32)
        image = read(os.path.join(testing_path, 'test_image.ms'))
        image = image[0].data.astype(np.float32)
        ccc = normxcorr2(template, image)[0]
        expected_ccc = np.load(os.path.join(testing_path, 'test_ccc.npy'))
        # We know that conda installs give a slightly different result
        self.assertTrue(np.allclose(expected_ccc, ccc, atol=0.003))
        # Differences occur for low correlation values, peak should be the same
        self.assertTrue(expected_ccc.max(), ccc.max())
        self.assertTrue(expected_ccc.argmax(), ccc.argmax())

    def test_failed_normxcorr(self):
        """Send it the wrong type."""
        ccc = normxcorr2(template=[0, 1, 2, 3, 4], image='bob')
        self.assertEqual(ccc, 'NaN')

    def test_spike_test(self):
        """Check that an error is raised!"""
        stream = read()
        stream[0].data[100] = 1e20
        with self.assertRaises(MatchFilterError):
            _spike_test(stream)


@pytest.mark.serial
class TestSynthData(unittest.TestCase):
    def test_threshold_methods(self):
        # Test other threshold methods
        for threshold_type, threshold in [('absolute', 2),
                                          ('av_chan_corr', 0.5)]:
            kfalse, ktrue = test_match_filter(threshold_type=threshold_type,
                                              threshold=threshold)
            self.assertTrue(kfalse / ktrue < 0.25)

    def test_missing_data(self):
        # Test case where there are non-matching streams in the template
        for conc_proc in [True, False]:
            test_match_filter(stream_excess=True, conc_proc=conc_proc)

    def test_extra_templates(self):
        # Test case where there are non-matching streams in the data
        for conc_proc in [True, False]:
            test_match_filter(template_excess=True, conc_proc=conc_proc)

    def test_onesamp_diff(self):
        """Tests to check that traces in stream are set to same length."""
        stream = Stream(traces=[
            Trace(data=np.random.randn(100)),
            Trace(data=np.random.randn(101))])
        stream[0].stats.sampling_rate = 40
        stream[0].stats.station = 'A'
        stream[1].stats.sampling_rate = 40
        stream[1].stats.station = 'B'
        templates = [Stream(traces=[Trace(data=np.random.randn(20)),
                                    Trace(data=np.random.randn(20))])]
        templates[0][0].stats.sampling_rate = 40
        templates[0][0].stats.station = 'A'
        templates[0][1].stats.sampling_rate = 40
        templates[0][1].stats.station = 'B'
        for conc_proc in [True, False]:
            match_filter(template_names=['1'], template_list=templates,
                         st=stream, threshold=8, threshold_type='MAD',
                         trig_int=1, plot=False,
                         concurrent_processing=conc_proc)

    def test_half_samp_diff(self):
        """
        Check that traces with different start-times by less than a sample
        are handled as expected.
        """
        stream = Stream(traces=[
            Trace(data=np.random.randn(100)),
            Trace(data=np.random.randn(101))])
        stream[0].stats.sampling_rate = 40
        stream[0].stats.station = 'A'
        stream[1].stats.sampling_rate = 40
        stream[1].stats.station = 'B'
        # Add some fraction of a sample to the starttime
        stream[0].stats.starttime += 0.25 * stream[0].stats.delta
        templates = [Stream(traces=[Trace(data=np.random.randn(20)),
                                    Trace(data=np.random.randn(20))])]
        templates[0][0].stats.sampling_rate = 40
        templates[0][0].stats.station = 'A'
        templates[0][1].stats.sampling_rate = 40
        templates[0][1].stats.station = 'B'
        for conc_proc in [True, False]:
            match_filter(template_names=['1'], template_list=templates,
                         st=stream, threshold=8, threshold_type='MAD',
                         trig_int=1, plot=False,
                         concurrent_processing=conc_proc)


@pytest.mark.network
class TestGeoNetCase(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        client = Client('GEONET')
        cls.t1 = UTCDateTime(2016, 9, 4)
        cls.t2 = cls.t1 + 86400
        catalog = client.get_events(
            starttime=cls.t1, endtime=cls.t2, minmagnitude=4, minlatitude=-49,
            maxlatitude=-35, minlongitude=175.0, maxlongitude=-175.0)
        catalog = catalog_utils.filter_picks(
            catalog, channels=['EHZ'], top_n_picks=5)
        for event in catalog:
            extra_pick = Pick()
            extra_pick.phase_hint = 'S'
            extra_pick.time = event.picks[0].time + 10
            extra_pick.waveform_id = event.picks[0].waveform_id
            event.picks.append(extra_pick)
        cls.tribe = Tribe()
        cls.tribe.construct(
            method='from_client', catalog=catalog, client_id='GEONET',
            lowcut=2.0, highcut=9.0, samp_rate=50.0, filt_order=4,
            length=3.0, prepick=0.15, swin='all', process_len=3600)
        cls.templates = [t.st for t in cls.tribe.templates]
        # Download and process the day-long data
        bulk_info = [(tr.stats.network, tr.stats.station, '*',
                      tr.stats.channel, cls.t1 + (4 * 3600),
                      cls.t1 + (5 * 3600))
                     for tr in cls.templates[0]]
        # Just downloading an hour of data
        print('Downloading data')
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        st.trim(cls.t1 + (4 * 3600), cls.t1 + (5 * 3600)).sort()
        # This is slow?
        print('Processing continuous data')
        cls.st = pre_processing.multi_process(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=50.0,
            num_cores=1)
        cls.st.trim(cls.t1 + (4 * 3600), cls.t1 + (5 * 3600)).sort()
        cls.template_names = [str(template[0].stats.starttime)
                              for template in cls.templates]

    def test_duplicate_channels_in_template(self):
        """
        Test using a template with duplicate channels.
        """
        templates = copy.deepcopy(self.templates)
        # Do this to test an extra condition in match_filter
        templates[0].remove(templates[0].select(station='CNGZ')[0])
        for conc_proc in [False, True]:
            Logger.info(f"Running for conc_proc={conc_proc}")
            detections = match_filter(
                template_names=self.template_names, template_list=templates,
                st=self.st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, plot=False, plotdir='.', cores=1,
                concurrent_processing=conc_proc)
            self.assertEqual(len(detections), 1)
            self.assertEqual(detections[0].no_chans, 6)

    def test_duplicate_cont_data(self):
        """ Check that error is raised if duplicate channels are present in
        the continuous data."""
        tr = self.st[0].copy()
        tr.data = np.random.randn(100)
        st = self.st.copy() + tr
        for conc_proc in [False, True]:
            Logger.info(f"Running for conc_proc={conc_proc}")
            with self.assertRaises(NotImplementedError):
                match_filter(template_names=self.template_names,
                             template_list=self.templates, st=st, threshold=8.0,
                             threshold_type='MAD', trig_int=6.0, plot=False,
                             plotdir='.', cores=1,
                             concurrent_processing=conc_proc)

    def test_missing_cont_channel(self):
        """ Remove one channel from continuous data and check that everything
        still works. """
        st = self.st.copy()
        st.remove(st[-1])
        for conc_proc in [False, True]:
            Logger.info(f"Running for conc_proc={conc_proc}")
            detections, det_cat = match_filter(
                template_names=self.template_names, template_list=self.templates,
                st=st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
                plot=False, plotdir='.', cores=1, output_cat=True,
                concurrent_processing=conc_proc)
            self.assertEqual(len(detections), 1)
            self.assertEqual(detections[0].no_chans, 5)
            self.assertEqual(len(detections), len(det_cat))

    def test_no_matching_data(self):
        """ No matching data between continuous and templates."""
        st = self.st.copy()
        for tr, staname in zip(st, ['a', 'b', 'c', 'd', 'e']):
            tr.stats.station = staname
        for conc_proc in [False, True]:
            Logger.info(f"Running for conc_proc={conc_proc}")
            with self.assertRaises(IndexError):
                match_filter(
                    template_names=self.template_names,
                    template_list=self.templates, st=st, threshold=8.0,
                    threshold_type='MAD', trig_int=6.0, plot=False,
                    plotdir='.', cores=1, concurrent_processing=conc_proc)

    @pytest.mark.flaky(reruns=2)
    def test_geonet_tribe_detect(self):
        client = Client('GEONET')
        # Try to force issues with starting samples on wrong day for geonet
        # data
        # TODO: This does nothing
        # tribe = self.tribe.copy()
        # for template in tribe.templates:
        #     template.process_length = 86400
        #     template.st = Stream(template.st[0])
        #     # Only run one channel templates
        for conc_proc in [False, True]:
            Logger.info(f"Running for conc_proc={conc_proc}")
            party = self.tribe.copy().client_detect(
                client=client, starttime=self.t1, endtime=self.t2,
                threshold=8.0, threshold_type='MAD', trig_int=6.0,
                daylong=False, plot=False, concurrent_processing=conc_proc)
            self.assertEqual(len(party), 16)


class TestGappyData(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        cls.client = Client("GEONET")
        cls.starttime = UTCDateTime(2016, 11, 13, 23)
        cls.endtime = UTCDateTime(2016, 11, 14)
        catalog = cls.client.get_events(
            starttime=UTCDateTime(2016, 11, 13, 23, 30),
            endtime=UTCDateTime(2016, 11, 13, 23, 31))
        catalog = filter_picks(catalog=catalog, stations=['KHZ'])
        cls.tribe = Tribe().construct(
            method="from_client", lowcut=2, highcut=8, samp_rate=20,
            filt_order=4, length=10, prepick=0.5, catalog=catalog,
            client_id="GEONET", process_len=3600, swin="P")
        cls.st = cls.client.get_waveforms(
            station="KHZ", network="NZ", channel="HHZ", location="10",
            starttime=cls.starttime, endtime=cls.endtime)

    def test_gappy_data(self):
        gaps = self.st.get_gaps()
        self.assertEqual(len(gaps), 1)
        start_gap = gaps[0][4]
        end_gap = gaps[0][5]
        for conc_proc in [True, False]:
            party = self.tribe.client_detect(
                client=self.client, starttime=self.starttime,
                endtime=self.endtime, threshold=0.6,
                threshold_type="absolute", trig_int=2, plot=False,
                parallel_process=False, cores=1,
                concurrent_processing=conc_proc)
            for family in party:
                print(family)
                for detection in family:
                    self.assertFalse(
                        start_gap <= detection.detect_time <= end_gap)
            for family in party:
                self.assertTrue(len(family) in [5, 1])

    def test_gappy_data_removal(self):
        for conc_proc in [True, False]:
            party = self.tribe.client_detect(
                client=self.client, starttime=self.starttime,
                endtime=self.endtime, threshold=8,
                threshold_type="MAD", trig_int=2, plot=False,
                parallel_process=False, min_gap=1,
                concurrent_processing=conc_proc)
            self.assertEqual(len(party), 0)


@pytest.mark.network
class TestNCEDCCases(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        client = Client('NCEDC')
        t1 = UTCDateTime(2004, 9, 28, 17)
        t2 = t1 + 3600
        process_len = 3600
        # t1 = UTCDateTime(2004, 9, 28)
        # t2 = t1 + 80000
        # process_len = 80000
        catalog = client.get_events(
            starttime=t1, endtime=t2, minmagnitude=4, minlatitude=35.7,
            maxlatitude=36.1, minlongitude=-120.6, maxlongitude=-120.2,
            includearrivals=True)
        catalog = catalog_utils.filter_picks(
            catalog, channels=['EHZ'], top_n_picks=5)
        cls.tribe = Tribe()
        cls.tribe.construct(
            method='from_client', catalog=catalog, client_id='NCEDC',
            lowcut=2.0, highcut=9.0, samp_rate=50.0, filt_order=4,
            length=3.0, prepick=0.15, swin='all', process_len=process_len)
        cls.templates = [t.st.copy() for t in cls.tribe]
        for template in cls.templates:
            template.sort()
        # Download and process the day-long data
        template_stachans = []
        for template in cls.templates:
            for tr in template:
                template_stachans.append(
                    (tr.stats.network, tr.stats.station, tr.stats.channel))
        template_stachans = list(set(template_stachans))
        bulk_info = [(stachan[0], stachan[1], '*', stachan[2],
                      t1, t1 + process_len + 1)
                     for stachan in template_stachans]
        # Just downloading an hour of data
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        cls.unproc_st = st.copy()
        cls.st = pre_processing.multi_process(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=50.0,
            num_cores=1, starttime=st[0].stats.starttime,
            endtime=st[0].stats.starttime + process_len)
        cls.template_names = [str(template[0].stats.starttime)
                              for template in cls.templates]

    def test_detection_extraction(self):
        # Test outputting the streams works
        for conc_proc in [True, False]:
            detections, detection_streams = match_filter(
                template_names=self.template_names,
                template_list=self.templates, st=self.st,
                threshold=8.0, threshold_type='MAD',
                trig_int=6.0, plot=False, plotdir='.',
                cores=1, extract_detections=True,
                concurrent_processing=conc_proc)
            self.assertEqual(len(detections), 4)
            self.assertEqual(len(detection_streams), len(detections))

    def test_normxcorr(self):
        # Test a known issue with early normalisation methods
        template_array = np.array(
            [t.select(station='PHOB', channel='EHZ')[0].data.astype(np.float32)
             for t in self.templates])
        # template_array = np.array([template_array[0]])
        stream = self.st.select(
            station='PHOB', channel='EHZ')[0].data.astype(np.float32)
        pads = [0 for _ in range(len(template_array))]
        ccc_numpy, no_chans = numpy_normxcorr(template_array, stream, pads)
        ccc, no_chans = fftw_normxcorr(template_array, stream, pads)
        self.assertTrue(np.allclose(ccc, ccc_numpy, atol=0.04))

    def test_catalog_extraction(self):
        for conc_proc in [True, False]:
            detections, det_cat, detection_streams = match_filter(
                template_names=self.template_names,
                template_list=self.templates, st=self.st,
                threshold=8.0, threshold_type='MAD',
                trig_int=6.0, plot=False, plotdir='.',
                cores=1, extract_detections=True, output_cat=True,
                concurrent_processing=conc_proc)
            self.assertEqual(len(detections), 4)
            self.assertEqual(len(detection_streams), len(detections))
            self.assertEqual(len(detection_streams), len(det_cat))

    def test_same_detections_individual_and_parallel(self):
        """
        Check that the same detections are made regardless of whether templates
        are run together or separately.
        """
        individual_detections = []
        for template, template_name in zip(self.templates,
                                           self.template_names):
            individual_detections += match_filter(
                template_names=[template_name], template_list=[template],
                st=self.st.copy(), threshold=8.0, threshold_type='MAD',
                trig_int=6.0, plot=False, plotdir='.', cores=1)
        individual_dict = []
        for detection in individual_detections:
            individual_dict.append({'template_name': detection.template_name,
                                    'time': detection.detect_time,
                                    'cccsum': detection.detect_val.round(4)})
        detections = match_filter(template_names=self.template_names,
                                  template_list=self.templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plot=False, plotdir='.',
                                  cores=1)
        self.assertEqual(len(individual_detections), len(detections))
        for detection in detections:
            detection_dict = {'template_name': detection.template_name,
                              'time': detection.detect_time,
                              'cccsum': detection.detect_val.round(4)}
            if detection_dict not in individual_dict:
                print(f"Detection:\n{detection_dict}\nnot found in:"
                      f"\n{individual_dict}")
            self.assertTrue(detection_dict in individual_dict)

    def test_read_write_detections(self):
        """Check that we can read and write detections accurately."""
        if os.path.isfile("dets_out.txt"):
            os.remove("dets_out.txt")
        detections = match_filter(template_names=self.template_names,
                                  template_list=self.templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plot=False, plotdir='.',
                                  cores=1)
        detection_dict = []
        for detection in detections:
            detection_dict.append(
                {'template_name': detection.template_name,
                 'time': detection.detect_time,
                 'cccsum': detection.detect_val,
                 'no_chans': detection.no_chans,
                 'chans': detection.chans,
                 'threshold': detection.threshold,
                 'typeofdet': detection.typeofdet})
            detection.write(fname='dets_out.txt')
        detections_back = read_detections('dets_out.txt')
        self.assertEqual(len(detections), len(detections_back))
        for detection in detections_back:
            d = {'template_name': detection.template_name,
                 'time': detection.detect_time,
                 'cccsum': detection.detect_val,
                 'no_chans': detection.no_chans,
                 'chans': detection.chans,
                 'threshold': detection.threshold,
                 'typeofdet': detection.typeofdet}
            self.assertTrue(d in detection_dict)
        os.remove('dets_out.txt')

    def test_get_catalog(self):
        """Check that catalog objects are created properly."""
        detections = match_filter(
            template_names=self.template_names, template_list=self.templates,
            st=self.st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
            plot=False, plotdir='.', cores=1)
        cat = get_catalog(detections)
        self.assertEqual(len(cat), len(detections))
        for det in detections:
            self.assertTrue(det.event in cat)
        # Check the short-cut for writing the catalog
        with NamedTemporaryFile() as tf:
            write_catalog(detections, fname=tf.name)
            cat_back = read_events(tf.name)
            self.assertEqual(len(cat), len(cat_back))

    def test_extraction(self):
        """Check the extraction function."""
        detections = match_filter(template_names=self.template_names,
                                  template_list=self.templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plot=False, plotdir='.',
                                  cores=1)
        streams = extract_from_stream(stream=self.st.copy(),
                                      detections=detections)
        self.assertEqual(len(streams), len(detections))
        # Test when a channel is missing in the stream
        streams = extract_from_stream(stream=self.st.copy()[0:-2],
                                      detections=detections)
        self.assertEqual(len(streams), len(detections))

    def test_incorrect_arguments(self):
        with self.assertRaises(MatchFilterError):
            # template_names is not a list
            match_filter(template_names=self.template_names[0],
                         template_list=self.templates, st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # templates is not a list
            match_filter(template_names=self.template_names,
                         template_list=self.templates[0], st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # template and template_names length are not equal
            match_filter(template_names=self.template_names,
                         template_list=[self.templates[0]], st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # templates is not a list of streams
            match_filter(template_names=self.template_names,
                         template_list=['abc'], st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # st is not a Stream
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=np.random.randn(10),
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # threshold_type is wrong
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=self.st,
                         threshold=8.0, threshold_type='albert', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)

    def test_masked_template(self):
        templates = [self.templates[0].copy()]
        tr = templates[0][0].copy()
        tr.stats.starttime += 3600
        templates[0] += tr
        templates[0].merge()
        with self.assertRaises(MatchFilterError):
            match_filter(template_names=[self.template_names[0]],
                         template_list=templates, st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)

    def test_non_equal_template_lengths(self):
        templates = [self.templates[0].copy()]
        templates[0][0].data = np.concatenate([templates[0][0].data,
                                               np.random.randn(10)])
        with self.assertRaises(AssertionError):
            match_filter(template_names=[self.template_names[0]],
                         template_list=templates, st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plot=False, plotdir='.', cores=1)


class TestMatchCopy(unittest.TestCase):
    def test_tribe_copy(self):
        """Test copy method"""
        party = Party().read(
            filename=os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                'test_data', 'test_party.tgz'))
        tribe = Tribe(f.template for f in party.families)
        copied = tribe.copy()
        self.assertEqual(len(tribe), len(copied))
        for t, copy_t in zip(tribe.templates, copied.templates):
            self.assertEqual(t, copy_t)
        self.assertEqual(tribe, copied)


@pytest.mark.network
class TestTribeConstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        client = Client("GEONET")
        cat = client.get_events(
            latitude=-41.7, longitude=174.1, maxradius=0.5, minmagnitude=4.0,
            starttime=UTCDateTime(2013, 7, 30),
            endtime=UTCDateTime(2013, 8, 2))
        cls.cat, cls.client = cat, client

    def test_short_stream(self):
        """ Check that events are correctly associated with templates, #381"""
        process_len = 3600
        length = 6
        swin = "P"

        st = self.client.get_waveforms(
            network="NZ", station="TCW", location="10", channel="EHZ",
            starttime=UTCDateTime(2013, 7, 31, 23),
            endtime=UTCDateTime(2013, 8, 1))
        tribe = Tribe().construct(
            method="from_meta_file", lowcut=2, highcut=10, samp_rate=50,
            filt_order=4, length=length, prepick=0.5, st=st,
            meta_file=self.cat, swin=swin, process_len=process_len)
        self.assertEqual(len(tribe), 2)
        for template in tribe:
            self.assertEqual(process_len, template.process_length)
            for tr in template.st:
                self.assertEqual(tr.stats.npts / tr.stats.sampling_rate,
                                 length)
                matched_pick = [p for p in template.event.picks
                                if p.waveform_id.get_seed_string() == tr.id
                                and p.phase_hint == swin]
                if len(matched_pick) == 0:
                    continue
                self.assertLessEqual(
                    abs(matched_pick[0].time -
                        tr.stats.starttime) - template.prepick, tr.stats.delta)


class TestMatchObjectHeavy(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        client = Client('NCEDC')
        t1 = UTCDateTime(2004, 9, 28, 17)
        process_len = 3600
        t2 = t1 + process_len
        catalog = client.get_events(
            starttime=t1, endtime=t2, minmagnitude=4,
            minlatitude=35.7, maxlatitude=36.1, minlongitude=-120.6,
            maxlongitude=-120.2, includearrivals=True)
        catalog = catalog_utils.filter_picks(
            catalog, channels=['EHZ'], top_n_picks=5)
        template_stachans = []
        for event in catalog:
            for pick in event.picks:
                template_stachans.append(
                    (pick.waveform_id.network_code,
                     pick.waveform_id.station_code,
                     pick.waveform_id.channel_code))
        template_stachans = list(set(template_stachans))
        bulk_info = [(stachan[0], stachan[1], '*', stachan[2],
                      t1 - 5, t2 + 5) for stachan in template_stachans]
        # Just downloading an hour of data
        try:
            st = client.get_waveforms_bulk(bulk_info)
        except FDSNException:
            st = Stream()
            for _bulk in bulk_info:
                st += client.get_waveforms(*_bulk)
        st.merge()
        st.trim(t1, t2)
        for tr in st:
            tr.data = tr.data[0:int(process_len * tr.stats.sampling_rate)]
            assert len(tr.data) == process_len * tr.stats.sampling_rate
            assert tr.stats.starttime - t1 < 0.1
        unproc_st = st.copy()
        tribe = Tribe().construct(
            method='from_meta_file', catalog=catalog, st=st.copy(),
            lowcut=2.0, highcut=9.0, samp_rate=20.0, filt_order=4,
            length=3.0, prepick=0.15, swin='all', process_len=process_len)
        print(tribe)
        onehztribe = Tribe().construct(
            method='from_meta_file', catalog=catalog, st=st.copy(),
            lowcut=0.1, highcut=0.45, samp_rate=1.0, filt_order=4,
            length=20.0, prepick=0.15, swin='all', process_len=process_len)
        st = pre_processing.multi_process(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=20.0,
            num_cores=1, starttime=st[0].stats.starttime,
            endtime=st[0].stats.starttime + process_len)
        party = Party().read(
            filename=os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                'test_data', 'test_party.tgz'))
        cls.family = party.sort()[0].copy()
        cls.t1, cls.t2, cls.template_stachans = (t1, t2, template_stachans)
        cls.unproc_st, cls.tribe, cls.onehztribe, cls.st, cls.party = (
            unproc_st, tribe, onehztribe, st, party)

    @classmethod
    def tearDownClass(cls):
        for f in ['eqcorrscan_temporary_party.pkl']:
            if os.path.isfile(f):
                os.remove(f)

    def test_tribe_detect(self):
        """Test the detect method on Tribe objects"""
        for conc_proc in [True, False]:
            party = self.tribe.detect(
                stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False,
                parallel_process=False, concurrent_processing=conc_proc)
            self.assertEqual(len(party), 4)
            compare_families(
                party=party, party_in=self.party, float_tol=0.05,
                check_event=True)

    def test_tribe_detect_with_empty_streams(self):
        """
        Compare the detect method for a tribe of one vs two templates and check
        that the detection has the same detect time, in a case where the
        continuous data is incomplete. This test should fail in v0.4.2 due to
        a bug.
        """
        for conc_proc in [True, False]:
            # remove trace for station PHA (PHOB, PSR, PCA, PAG remain)
            st = self.unproc_st.copy().remove(
                self.unproc_st.copy().select(station='PHA')[0])
            tribe1 = Tribe([t.copy() for t in self.tribe
                            if (t.name == '2004_09_28t17_19_08' or
                                t.name == '2004_09_28t17_19_25')])
            # run detection with 2 templates in tribe
            party1 = tribe1.detect(
                stream=st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plotvar=False,
                parallel_process=False, concurrent_processing=conc_proc)
            self.assertEqual(len(party1), 2)
            party1 = Party([f for f in party1
                            if f.template.name == '2004_09_28t17_19_25'])
            # run detection with only 1 template in tribe
            tribe2 = Tribe([t.copy() for t in self.tribe
                            if t.name == '2004_09_28t17_19_25'])
            party2 = tribe2.detect(
                stream=st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plotvar=False,
                parallel_process=False, concurrent_processing=conc_proc)
            self.assertEqual(len(party2), 1)
            # This should fail in v0.4.2
            compare_families(
                party=party1, party_in=party2, float_tol=0.05,
                check_event=False)

    def test_tribe_detect_short_data(self):
        """Test the detect method on Tribe objects"""
        for conc_proc in [True, False]:
            short_st = self.unproc_st.copy()
            tribe = self.tribe.copy()
            for template in tribe:
                template.process_length = 2400
            party = tribe.detect(
                stream=short_st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False,
                parallel_process=False, concurrent_processing=conc_proc,
                ignore_bad_data=True)
            self.assertEqual(len(party), 4)

    @pytest.mark.serial
    def test_tribe_detect_parallel_process(self):
        """Test the detect method on Tribe objects"""
        for conc_proc in [True, False]:
            party = self.tribe.detect(
                stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False, parallel_process=True,
                process_cores=2, concurrent_processing=conc_proc)
            self.assertEqual(len(party), 4)
            compare_families(
                party=party, party_in=self.party, float_tol=0.05,
                check_event=False)

    def test_tribe_detect_save_progress(self):
        """Test the detect method on Tribe objects"""
        for conc_proc in [True, False]:
            party = self.tribe.detect(
                stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False, parallel_process=False,
                save_progress=True, concurrent_processing=conc_proc)
            self.assertEqual(len(party), 4)
            # Get all the parties
            party_files = glob.glob(".parties/????/???/*.pkl")
            saved_party = Party()
            for pf in party_files:
                saved_party += Party().read(pf)
            self.assertEqual(party, saved_party)
            shutil.rmtree(".parties")

    @pytest.mark.serial
    def test_tribe_detect_masked_data(self):
        """Test using masked data - possibly raises error at pre-processing.
        Padding may also result in error at correlation stage due to poor
        normalisation."""
        for conc_proc in [True, False]:
            stream = self.unproc_st.copy()
            stream[0] = (stream[0].copy().trim(
                stream[0].stats.starttime, stream[0].stats.starttime + 1800) +
                         stream[0].trim(
                stream[0].stats.starttime + 1900, stream[0].stats.endtime))
            party = self.tribe.detect(
                stream=stream, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False, parallel_process=False,
                xcorr_func='fftw', concurrency='concurrent',
                concurrent_processing=conc_proc)
            self.assertEqual(len(party), 4)

    def test_tribe_detect_no_processing(self):
        """Test that no processing is done when it isn't necessary."""
        for conc_proc in [True, False]:
            tribe = self.tribe.copy()
            for template in tribe:
                template.lowcut = None
                template.highcut = None
            party = tribe.detect(
                stream=self.st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False,
                parallel_process=False, concurrent_processing=conc_proc)
            self.assertEqual(len(party), 4)
            compare_families(
                party=party, party_in=self.party, float_tol=0.05,
                check_event=False)

    @pytest.mark.flaky(reruns=2)
    @pytest.mark.network
    def test_client_detect(self):
        """Test the client_detect method."""
        for conc_proc in [True, False]:
            client = Client('NCEDC')
            party = self.tribe.copy().client_detect(
                client=client, starttime=self.t1 + 2.75, endtime=self.t2,
                threshold=8.0, threshold_type='MAD', trig_int=6.0,
                daylong=False, plot=False, concurrent_processing=conc_proc)
            compare_families(
                party=party, party_in=self.party, float_tol=0.05,
                check_event=False)

    @pytest.mark.flaky(reruns=2)
    @pytest.mark.network
    def test_client_detect_save_progress(self):
        """Test the client_detect method."""
        for conc_proc in [True, False]:
            client = Client('NCEDC')
            party = self.tribe.copy().client_detect(
                client=client, starttime=self.t1 + 2.75, endtime=self.t2,
                threshold=8.0, threshold_type='MAD', trig_int=6.0,
                daylong=False, plot=False, save_progress=True,
                concurrent_processing=conc_proc)
            self.assertTrue(os.path.isdir(".parties"))

            # Get all the parties
            party_files = glob.glob(".parties/????/???/*.pkl")
            saved_party = Party()
            for pf in party_files:
                saved_party += Party().read(pf)
            self.assertEqual(party, saved_party)
            shutil.rmtree(".parties")
            compare_families(
                party=party, party_in=self.party, float_tol=0.05,
                check_event=False)

    @pytest.mark.network
    def test_party_lag_calc(self):
        """Test the lag-calc method on Party objects."""
        # Test the chained method
        for conc_proc in [True, False]:
            chained_cat = self.tribe.detect(
                stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False,
                concurrent_processing=conc_proc).lag_calc(
                stream=self.unproc_st, pre_processed=False)
            catalog = self.party.copy().lag_calc(
                stream=self.unproc_st, pre_processed=False)
            self.assertEqual(len(catalog), 4)
            for ev1, ev2 in zip(catalog, chained_cat):
                ev1.picks.sort(key=lambda p: p.time)
                ev2.picks.sort(key=lambda p: p.time)
            catalog.events.sort(key=lambda e: e.picks[0].time)
            chained_cat.events.sort(key=lambda e: e.picks[0].time)
            for ev, chained_ev in zip(catalog, chained_cat):
                for i in range(len(ev.picks)):
                    for key in ev.picks[i].keys():
                        if key == 'resource_id':
                            continue
                        if key == 'comments':
                            continue
                        if key == 'waveform_id':
                            for _k in ['network_code', 'station_code',
                                       'channel_code']:
                                self.assertEqual(
                                    sorted(ev.picks,
                                           key=lambda p: p.time)[i][key][_k],
                                    sorted(chained_ev.picks,
                                           key=lambda p: p.time)[i][key][_k])
                            continue
                        self.assertEqual(
                            sorted(ev.picks,
                                   key=lambda p: p.time)[i][key],
                            sorted(chained_ev.picks,
                                   key=lambda p: p.time)[i][key])
                    pick_corrs = sorted(ev.picks, key=lambda p: p.time)
                    pick_corrs = [float(p.comments[0].text.split("=")[-1])
                                  for p in pick_corrs]
                    chained_ev_pick_corrs = sorted(ev.picks, key=lambda p: p.time)
                    chained_ev_pick_corrs = [
                        float(p.comments[0].text.split("=")[-1])
                        for p in chained_ev_pick_corrs]
                    assert np.allclose(
                        pick_corrs, chained_ev_pick_corrs, atol=0.001)
                    assert np.allclose(
                        float(ev.comments[0].text.split("=")[-1]),
                        float(chained_ev.comments[0].text.split("=")[-1]),
                        atol=0.001)

    def test_party_lag_calc_preprocessed(self):
        """Test that the lag-calc works on pre-processed data."""
        catalog = self.party.copy().lag_calc(
            stream=self.st, pre_processed=True)
        self.assertEqual(len(catalog), 4)

    def test_party_lag_calc_empty_families(self):
        """ Test that passing some empty families works properly, see #341. """
        party = self.party.copy()
        party[1].detections = []  # Make an empty family
        catalog = party.lag_calc(stream=self.st, pre_processed=True)
        self.assertEqual(len(catalog), 3)

    def test_party_lag_calc_missing_data(self):
        """Check that if data are short, then no events are returned. #406 """
        party = self.party.copy()
        st = self.unproc_st.copy()
        st = st.trim(st[0].stats.starttime,
                     st[0].stats.starttime + (
                         0.75 * party[0].template.process_length))
        catalog = party.lag_calc(stream=st, pre_processed=False)
        self.assertEqual(len(catalog), 0)

    def test_party_lag_calc_short_data(self):
        """ Check that insufficient data with ignore_length are run. #406 """
        party = self.party.copy()
        st = self.unproc_st.copy()
        cut_start = st[0].stats.starttime + (
                0.5 * party[0].template.process_length)
        cut_end = st[0].stats.starttime + (
                0.8 * party[0].template.process_length)
        st = st.cutout(cut_start, cut_end)
        catalog = party.lag_calc(stream=st, pre_processed=False,
                                 ignore_length=True, ignore_bad_data=True)
        self.assertEqual(len(catalog), 4)

    # TODO: tests removed due to method removal - re-instate when relative
    #  magnitudes is viable.
    # def test_party_mag_calc_unpreprocessed(self):
    #     """Test that the lag-calc works on pre-processed data."""
    #     catalog = self.party.copy().lag_calc(
    #         stream=self.unproc_st, pre_processed=False,
    #         relative_magnitudes=True, min_snr=0)
    #     self.assertEqual(len(catalog), 4)
    #     for event in catalog:
    #         self.assertGreater(len(event.station_magnitudes), 0)
    #         template_id = [c.text.split(": ")[-1] for c in event.comments
    #                        if "Template" in c.text][0]
    #         template = [fam.template for fam in self.party
    #                     if fam.template.name == template_id][0]
    #         self.assertAlmostEqual(
    #             template.event.magnitudes[0].mag,
    #             event.magnitudes[0].mag, 1)
    #         self.assertEqual(
    #             template.event.magnitudes[0].magnitude_type,
    #             event.magnitudes[0].magnitude_type)
    #
    # def test_party_mag_calc_preprocessed(self):
    #     """Test that the lag-calc works on pre-processed data."""
    #     catalog = self.party.copy().lag_calc(
    #         stream=self.st, pre_processed=True, relative_magnitudes=True,
    #         min_snr=0)
    #     self.assertEqual(len(catalog), 4)
    #     for event in catalog:
    #         self.assertGreater(len(event.station_magnitudes), 0)
    #         print(event)
    #         print(event.comments)
    #         template_id = [c.text.split(": ")[-1] for c in event.comments
    #                        if "Template" in c.text][0]
    #         template = [fam.template for fam in self.party
    #                     if fam.template.name == template_id][0]
    #         self.assertAlmostEqual(
    #             template.event.magnitudes[0].mag,
    #             event.magnitudes[0].mag, 1)
    #         self.assertEqual(
    #             template.event.magnitudes[0].magnitude_type,
    #             event.magnitudes[0].magnitude_type)

    @pytest.mark.network
    def test_day_long_methods(self):
        """Conduct a test using day-long data."""
        client = Client('NCEDC', timeout=300)
        t1 = UTCDateTime(2004, 9, 28)
        bulk_info = [(stachan[0], stachan[1], '*', stachan[2], t1, t1 + 86400)
                     for stachan in self.template_stachans]
        # Just downloading an hour of data
        print('Downloading continuous day-long data')
        try:
            st = client.get_waveforms_bulk(bulk_info)
        except FDSNException:
            st = Stream()
            for _bulk in bulk_info:
                st += client.get_waveforms(*_bulk)
        st.merge(fill_value='interpolate')
        # Hack day-long templates
        daylong_tribe = self.onehztribe.copy()
        for template in daylong_tribe:
            template.process_length = 86400
        # Aftershock sequence, with 1Hz data, lots of good correlations = high
        # MAD!
        Logger.info(f"Downloaded {len(st)} traces - handing off to detect")
        for conc_proc in [True, False]:
            Logger.info(f"Running conc_proc={conc_proc}")
            day_party = daylong_tribe.detect(
                stream=st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
                daylong=True, plot=False, parallel_process=False,
                concurrent_processing=conc_proc)
            self.assertEqual(len(day_party), 4)
            day_catalog = day_party.lag_calc(stream=st, pre_processed=False,
                                             parallel=False)
            self.assertEqual(len(day_catalog), 4)
            pre_picked_cat = day_party.get_catalog()
            self.assertEqual(len(pre_picked_cat), 4)

    def test_family_lag_calc(self):
        """Test the lag-calc method on family."""
        catalog = self.family.copy().lag_calc(
            stream=self.st, pre_processed=True)
        self.assertEqual(len(catalog), 1)

    def test_template_detect(self):
        """Test detect method on Template objects."""
        test_template = self.family.template.copy()
        for conc_proc in [True, False]:
            party_t = test_template.detect(
                stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
                trig_int=6.0, daylong=False, plot=False, overlap=None,
                concurrent_processing=conc_proc)
            self.assertEqual(len(party_t), 1)

    def test_template_construct_not_implemented(self):
        """Test template construction."""
        with self.assertRaises(NotImplementedError):
            Template().construct(
                method='from_client', client_id='NCEDC', name='bob',
                lowcut=2.0, highcut=9.0, samp_rate=50.0, filt_order=4,
                length=3.0, prepick=0.15, swin='all', process_len=6)


class TestMatchObjectLight(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.party = Party().read(
            filename=os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                'test_data', 'test_party.tgz'))
        cls.tribe = Tribe(templates=[fam.template for fam in cls.party])
        cls.family = cls.party.sort()[0].copy()

    @pytest.mark.mpl_image_compare
    def test_party_plot_individual(self):
        fig = self.party.plot(show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_party_plot_grouped(self):
        fig = self.party.plot(
            plot_grouped=True, show=False, return_figure=True)
        return fig

    @pytest.mark.mpl_image_compare
    def test_party_plot_grouped_rate(self):
        fig = self.party.plot(
            plot_grouped=True, rate=True, show=False, return_figure=True)
        return fig

    def test_party_io_list(self):
        """Test reading and writing party objects."""
        if os.path.isfile('test_party_list.tgz'):
            os.remove('test_party_list.tgz')
        try:
            self.party.write(filename='test_party_list')
            party_back = read_party(fname=['test_party_list.tgz'])
            self.assertEqual(self.party, party_back)
        finally:
            if os.path.isfile('test_party_list.tgz'):
                os.remove('test_party_list.tgz')

    def test_party_io_wildcards(self):
        """Test reading and writing party objects."""
        party_0 = self.party.copy()
        for f in party_0:
            for d in f:
                detect_time = d.detect_time + 3600
                d.detect_time = detect_time
                d.id = (''.join(d.template_name.split(' ')) + '_' +
                        detect_time.strftime('%Y%m%d_%H%M%S%f'))

        for f in ["0", "1"]:
            if os.path.isfile(f'test_party_walrus_{f}.tgz'):
                os.remove(f'test_party_walrus_{f}.tgz')
        try:
            self.party.write(filename='test_party_walrus_0')
            party_0.write(filename='test_party_walrus_1')
            party_back = read_party(fname='test_party_w*.tgz')
            self.assertEqual(self.party + party_0, party_back)
        finally:
            for f in ["0", "1"]:
                if os.path.isfile(f'test_party_walrus_{f}.tgz'):
                    os.remove(f'test_party_walrus_{f}.tgz')

    def test_tribe_internal_methods(self):
        self.assertEqual(len(self.tribe), 4)
        self.assertTrue(self.tribe == self.tribe)
        self.assertFalse(self.tribe != self.tribe)

    def test_tribe_add(self):
        """Test add method"""
        added = self.tribe.copy()
        self.assertEqual(len(added + added[0]), 5)
        self.assertEqual(len(added), 4)
        added += added[-1]
        self.assertEqual(len(added), 5)

    def test_tribe_remove(self):
        """Test remove method"""
        removal = self.tribe.copy()
        self.assertEqual(len(removal.remove(removal[0])), 3)
        for template in self.tribe:
            self.assertTrue(isinstance(template, Template))

    def test_tribe_io_qml(self):
        """Test reading and writing or Tribe objects using tar form."""
        try:
            if os.path.isfile('test_tribe_QML.tgz'):
                os.remove('test_tribe_QML.tgz')
            self.tribe.write(
                filename='test_tribe_QML', catalog_format="QUAKEML")
            tribe_back = read_tribe('test_tribe_QML.tgz')
            self.assertEqual(self.tribe, tribe_back)
        finally:
            if os.path.isfile('test_tribe_QML.tgz'):
                os.remove('test_tribe_QML.tgz')

    def test_tribe_io_sc3ml(self):
        """Test reading and writing or Tribe objects using tar form."""
        try:
            if os.path.isfile('test_tribe_SC3ML.tgz'):
                os.remove('test_tribe_SC3ML.tgz')
            self.tribe.write(
                filename='test_tribe_SC3ML', catalog_format="SC3ML")
            tribe_back = read_tribe('test_tribe_SC3ML.tgz')
            for template_in, template_back in zip(self.tribe, tribe_back):
                assert template_in.__eq__(
                    template_back, verbose=True, shallow_event_check=True)
        finally:
            if os.path.isfile('test_tribe_SC3ML.tgz'):
                os.remove('test_tribe_SC3ML.tgz')

    # Requires bug-fixes in obspy to be deployed.
    # def test_tribe_io_nordic(self):
    #     """Test reading and writing or Tribe objects using tar form."""
    #     try:
    #         if os.path.isfile('test_tribe_nordic.tgz'):
    #             os.remove('test_tribe_nordic.tgz')
    #         self.tribe.write(
    #             filename='test_tribe_nordic', catalog_format="NORDIC")
    #         tribe_back = read_tribe('test_tribe_nordic.tgz')
    #         for template_in, template_back in zip(self.tribe, tribe_back):
    #             assert template_in.__eq__(
    #                 template_back, verbose=True, shallow_event_check=True)
    #     finally:
    #         os.remove('test_tribe_nordic.tgz')

    def test_detection_regenerate_event(self):
        template = self.party[0].template
        test_detection = self.party[0][0]
        test_detection_altered = test_detection.copy()
        # Make sure they are equal going into this.
        self.assertEqual(test_detection, test_detection_altered)
        test_detection_altered._calculate_event(
            template=template, estimate_origin=False)
        self.assertEqual(test_detection, test_detection_altered)
        test_detection_altered._calculate_event(
            template_st=template.st, estimate_origin=False)
        # Picks are adjusted by prepick
        for pick in test_detection_altered.event.picks:
            pick.time += template.prepick
        self.assertEqual(test_detection, test_detection_altered)
        # Check that detection is left alone if wrong template given
        test_detection_altered._calculate_event(
            template=self.party[1].template)
        # Check that the origin gets assigned with a useful time
        test_detection_altered._calculate_event(template=template)
        self.assertAlmostEqual(
            test_detection_altered.event.origins[0].time,
            template.event.origins[0].time, 1)
        for pick in test_detection_altered.event.picks:
            matched_pick = [p for p in template.event.picks
                            if p.waveform_id == pick.waveform_id and
                            p.phase_hint == pick.phase_hint]
            self.assertEqual(len(matched_pick), 1)

    def test_party_basic_methods(self):
        """Test the basic methods on Party objects."""
        self.assertEqual(self.party.__repr__(), 'Party of 4 Families.')
        self.assertFalse(self.party != self.party)

    def test_party_add(self):
        """Test getting items and adding them to party objects, and sorting"""
        test_party = self.party.copy()
        test_family = test_party[0]
        self.assertTrue(isinstance(test_family, Family))
        self.assertEqual(len(test_party + test_family), 5)
        self.assertEqual(len(test_party), 4)
        test_party += test_family
        self.assertEqual(len(test_party), 5)
        test_slice = test_party[0:2]
        # Add new fam with fake det to ensure unique families aren't missed
        new_family = Family(
            template=self.tribe[-1].copy(),
            detections=[Detection(template_name=self.tribe[-1].name,
                                  detect_time=UTCDateTime(), no_chans=5,
                                  detect_val=3.5, threshold=2.0,
                                  typeofdet='corr', threshold_type='MAD',
                                  threshold_input=8.0)],
            catalog=Catalog(events=[Event()]))
        test_slice.families.append(new_family)
        self.assertTrue(isinstance(test_slice, Party))
        self.assertEqual(len(test_slice), 4)
        self.assertEqual(len(test_party + test_slice), 9)
        self.assertEqual(len(test_party), 5)
        test_party += test_slice
        self.assertEqual(len(test_party), 9)
        with self.assertRaises(NotImplementedError):
            test_party += ['bob']
        test_party.sort()
        self.assertEqual(
            [f.template.name for f in test_party.families],
            sorted([f.template.name for f in test_party.families]))

    def test_party_decluster(self):
        """Test the decluster method on party."""
        trig_ints = [40, 15, 3600]
        trig_ints = [3600]
        for trig_int in trig_ints:
            for metric in ['avg_cor', 'cor_sum']:
                declust = self.party.copy().decluster(
                    trig_int=trig_int, metric=metric)
                declustered_dets = [
                    d for family in declust for d in family.detections]
                for det in declustered_dets:
                    for d in declustered_dets:
                        if d == det:
                            continue
                        dif = abs(det.detect_time - d.detect_time)
                        print('Time dif: {0} trig_int: {1}'.format(
                            dif, trig_int))
                        self.assertGreater(dif, trig_int)
                with self.assertRaises(IndexError):
                    self.party.copy().decluster(
                        trig_int=trig_int, timing='origin', metric=metric)

    def test_party_decluster_same_times(self):
        """
        Test that the correct detection is associated with the peak.
        Tests for the case where two detections from different templates are
        made at the same time - the peak-finding finds the best, but decluster
        did not always correctly associate the correct detection with that
        peak.
        """
        # Test insertion before
        test_party = self.party.copy()
        det = test_party[0][0].copy()
        det.detect_time = test_party[1][0].detect_time
        det.detect_val = 4
        test_party[0].detections.append(det)
        test_party.decluster(1)
        assert det not in [d for f in test_party for d in f]
        # Tes insertion after
        test_party = self.party.copy()
        det = test_party[1][0].copy()
        det.detect_time = test_party[0][0].detect_time
        det.detect_val = 4
        test_party[1].detections.append(det)
        test_party.decluster(1)
        assert det not in [d for f in test_party for d in f]

    def test_party_rethreshold(self):
        """Make sure that rethresholding removes the events we want it to."""
        party = self.party.copy()
        # Append a load of detections to the first family
        for i in range(200):
            det = party[0][0].copy()
            det.detect_time += i * 20
            # Include some negative detections!
            fudge = np.random.randint(-1, 2)
            if fudge == 0:
                fudge += 1
            det.detect_val = det.threshold + (i * 1e-1)
            det.detect_val *= fudge
            det.id = str(i)
            party[0].detections.append(det)
        self.assertEqual(len(party), 204)
        negative_count = len([d for f in party for d in f if d.detect_val < 0])
        self.assertGreater(negative_count, 0)
        party1 = party.copy().rethreshold(new_threshold=9)
        for family in party1:
            for d in family:
                self.assertEqual(d.threshold_input, 9.0)
                self.assertGreaterEqual(d.detect_val, d.threshold)
        party2 = party.copy().rethreshold(new_threshold=9, abs_values=True)
        negative_count = 0
        for family in party2:
            for d in family:
                self.assertEqual(d.threshold_input, 9.0)
                self.assertGreaterEqual(abs(d.detect_val), d.threshold)
                if d.detect_val < 0:
                    negative_count += 1
        # Check that there actually are some negative detections...
        self.assertGreater(negative_count, 0)

    def test_family_init(self):
        """Test generating a family with various things."""
        test_family = Family(template=self.family.template.copy())
        self.assertEqual(len(test_family), 0)
        test_family = Family(template=self.family.template.copy(),
                             detections=self.family.detections)
        self.assertEqual(len(test_family), 1)
        test_family = Family(template=self.family.template.copy(),
                             detections=self.family.detections[0])
        self.assertEqual(len(test_family), 1)
        test_family = Family(template=self.family.template.copy(),
                             detections=self.family.detections,
                             catalog=self.family.catalog)
        self.assertEqual(len(test_family), 1)
        test_family = Family(template=self.family.template.copy(),
                             detections=self.family.detections,
                             catalog=self.family.catalog[0])
        self.assertEqual(len(test_family), 1)

    def test_template_printing(self):
        """Check that printing templates works."""
        test_template = Template()
        self.assertEqual(test_template.__repr__(), 'Template()')
        test_template.name = 'bob'
        self.assertTrue('Template bob' in test_template.__repr__())
        test_template = self.family.template.copy()
        self.assertTrue('Template' in test_template.__repr__())

    def test_slicing_by_name(self):
        """Check that slicing by template name works as expected"""
        t_name = self.tribe[2].name
        self.assertEqual(self.tribe[2], self.tribe[t_name])
        t_name = self.party[2].template.name
        self.assertEqual(self.party[2], self.party[t_name])

    def test_template_io(self):
        """Test template read/write."""
        test_template = self.family.template.copy()
        try:
            test_template.write('test_template')
            template_back = Template().read('test_template.tgz')
            self.assertEqual(test_template, template_back)
        finally:
            if os.path.isfile('test_template.tgz'):
                os.remove('test_template.tgz')
        # Make sure we raise a useful error when trying to read from a
        # tribe file.
        try:
            self.tribe.write(filename='test_template')
            with self.assertRaises(IOError):
                Template().read('test_template.tgz')
        finally:
            if os.path.isfile('test_template.tgz'):
                os.remove('test_template.tgz')

    def test_template_grouping(self):
        # PR #524
        # Test that this works directly on the tribe - it should
        tribe_len = len(self.tribe)
        groups = quick_group_templates(self.tribe)
        self.assertEqual(len(groups), 1)
        # Add one copy of a template with a different processing length
        t2 = self.tribe[0].copy()
        t2.process_length -= 100
        templates = [t2]
        templates.extend(self.tribe.templates)
        # Quick check that we haven't changed the tribe
        self.assertEqual(len(self.tribe), tribe_len)
        groups2 = quick_group_templates(templates)
        self.assertEqual(len(groups2), 2)

    def test_party_io(self):
        """Test reading and writing party objects."""
        if os.path.isfile('test_party_out.tgz'):
            os.remove('test_party_out.tgz')
        try:
            self.party.write(filename='test_party_out')
            party_back = read_party(fname='test_party_out.tgz')
            self.assertEqual(self.party, party_back)
        finally:
            if os.path.isfile('test_party_out.tgz'):
                os.remove('test_party_out.tgz')

    def test_party_write_csv(self):
        """ There was an issue (#298) where the header was repeated."""
        party = Party().read()
        party.write("test_party.csv", format="csv")
        with open("test_party.csv", "rb") as f:
            lines = f.read().decode().split("\n")
        self.assertEqual(len(lines), 6)
        os.remove("test_party.csv")

    def test_party_io_no_catalog_writing(self):
        """Test reading and writing party objects."""
        if os.path.isfile('test_party_out_no_cat.tgz'):
            os.remove('test_party_out_no_cat.tgz')
        try:
            self.party.write(
                filename='test_party_out_no_cat',
                write_detection_catalog=False)
            party_back = read_party(fname='test_party_out_no_cat.tgz',
                                    estimate_origin=False)
            self.assertTrue(self.party.__eq__(party_back, verbose=True))
        finally:
            if os.path.isfile('test_party_out_no_cat.tgz'):
                os.remove('test_party_out_no_cat.tgz')

    def test_party_io_no_catalog_reading(self):
        """Test reading and writing party objects."""
        if os.path.isfile('test_party_out_no_cat2.tgz'):
            os.remove('test_party_out_no_cat2.tgz')
        try:
            self.party.write(filename='test_party_out_no_cat2')
            party_back = read_party(
                fname='test_party_out_no_cat2.tgz',
                read_detection_catalog=False, estimate_origin=False)
            # creation times will differ - hack around this to make comparison
            # easier
            for family in self.party:
                family_back = party_back.select(family.template.name)
                for detection_in, detection_back in zip(family, family_back):
                    detection_in.event.creation_info.creation_time = \
                        detection_back.event.creation_info.creation_time
            self.assertTrue(self.party.__eq__(party_back, verbose=True))
        finally:
            if os.path.isfile('test_party_out_no_cat2.tgz'):
                os.remove('test_party_out_no_cat2.tgz')

    def test_family_methods(self):
        """Test basic methods on Family objects."""
        family = self.family.copy()
        self.assertEqual(
            family.__repr__(),
            'Family of 1 detections from template 2004_09_28t17_15_26')

    def test_family_addition(self):
        """Test adding to the family."""
        family = self.family.copy()
        fam_copy = family.copy()
        fam_copy.template.name = 'bob'
        with self.assertRaises(NotImplementedError):
            family += fam_copy
        with self.assertRaises(NotImplementedError):
            family += 'bob'

    def test_family_equality(self):
        """Test that when we check equality all is good."""
        family = self.family.copy()
        fam_copy = family.copy()
        fam_copy.template.name = 'bob'
        self.assertFalse(family == fam_copy)
        self.assertTrue(family != fam_copy)
        fam_copy.template = family.template
        fam_copy.detections = []
        self.assertFalse(family == fam_copy)

    def test_family_slicing(self):
        """Check getting items returns the expected result."""
        family = self.family.copy()
        self.assertTrue(isinstance(family[0], Detection))
        self.assertTrue(isinstance(family[0:], list))

    def test_family_sort(self):
        """Test sorting of family objects."""
        family = self.family.copy()
        for i in np.arange(20):
            randn_det = family.detections[0].copy()
            randn_det.detect_time += np.random.rand() * i * 100
            family.detections.append(randn_det)
        sorted_dets = family.copy().detections
        sorted_dets.sort(key=lambda x: x.detect_time)
        self.assertEqual(family.sort().detections, sorted_dets)

    def test_family_io(self):
        """Test the write method of family."""
        family = self.family.copy()
        try:
            family.write('test_family', overwrite=True)
            party_back = read_party('test_family.tgz')
            self.assertEqual(len(party_back), 1)
            self.assertEqual(party_back[0], family)
        finally:
            if os.path.isfile('test_family.tgz'):
                os.remove('test_family.tgz')

    def test_family_catalogs(self):
        """Check that the catalog always represents the detections"""
        family = self.family.copy()
        self.assertEqual(family.catalog, get_catalog(family.detections))
        additional_detection = family.detections[0].copy()
        additional_detection.detect_time += 3600
        for pick in additional_detection.event.picks:
            pick.time += 3600
        added_family = family + additional_detection
        self.assertEqual(added_family.catalog,
                         get_catalog(added_family.detections))
        family.detections.append(additional_detection)
        self.assertEqual(family.catalog, get_catalog(family.detections))


class TestTemplateGrouping(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        templates = []
        station_names = [
            'ALPH', 'BETA', 'GAMM', 'KAPP', 'ZETA', 'BOB', 'MAGG',
            'ALF', 'WALR', 'ALBA', 'PENG', 'BANA', 'WIGG', 'SAUS',
            'MALC']
        for i in range(20):
            template = Stream()
            for j in range(10):
                for c in ["EHZ", "EHN", "EHE"]:
                    tr = Trace(
                        data=np.random.randn(600),
                        header=dict(station=station_names[j],
                                    channel=c, network="NZ", location="10",
                                    sampling_rate=100.,
                                    starttime=UTCDateTime() + (i * 1000)))
                    template += tr
            templates.append(template)
        st = templates[0].copy()
        for tr in st:
            tr.data = np.random.randn(3600 * 100)
        templates = [Template(name=str(i), st=t)
                     for i, t in enumerate(templates)]
        cls.templates = templates
        cls.st = st
        cls.st_seed_ids = {tr.id for tr in st}

    def test_all_grouped(self):
        groups = group_templates_by_seedid(
            templates=self.templates, st_seed_ids=self.st_seed_ids,
            group_size=100)
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]), len(self.templates))

    def test_group_size_respected(self):
        groups = group_templates_by_seedid(
            templates=self.templates, st_seed_ids=self.st_seed_ids,
            group_size=10)
        self.assertEqual(len(groups), 2)
        self.assertEqual(len(groups[0]) + len(groups[1]), len(self.templates))

    def test_all_different_seeds(self):
        edited_templates = []
        rng = np.random.default_rng()
        for template in self.templates:
            template = template.copy()
            choices = rng.choice(len(template.st), 5)
            template.st = [tr for i, tr in enumerate(template.st)
                           if i in choices]
            edited_templates.append(template)
        groups = group_templates_by_seedid(
            templates=edited_templates, st_seed_ids=self.st_seed_ids,
            group_size=10)
        self.assertEqual(len(groups), 2)
        self.assertEqual(len(groups[0]) + len(groups[1]), len(self.templates))

    def test_unmatched_dropped(self):
        edited_templates = [t.copy() for t in self.templates]
        for tr in edited_templates[0].st:
            tr.stats.channel = "ABC"

        groups = group_templates_by_seedid(
             templates=edited_templates, st_seed_ids=self.st_seed_ids,
             group_size=10)
        self.assertEqual(len(groups), 2)
        self.assertEqual(len(groups[0]) + len(groups[1]),
                         len(self.templates) - 1)


def compare_families(party, party_in, float_tol=0.001, check_event=True):
    party.sort()
    party_in.sort()
    for fam, check_fam in zip(party, party_in):
        assert fam.template.name == check_fam.template.name
        fam.detections.sort(key=lambda d: d.detect_time)
        check_fam.detections.sort(key=lambda d: d.detect_time)
        for det, check_det in zip(fam.detections, check_fam.detections):
            for key in det.__dict__.keys():
                if key == 'event':
                    if not check_event:
                        continue
                    assert (
                        len(det.__dict__[key].picks) ==
                        len(check_det.__dict__[key].picks))
                    # Check that the number of picks equals the number of
                    # traces in the template
                    assert len(det.__dict__[key].picks) == len(fam.template.st)
                    min_template_time = min(
                        [tr.stats.starttime for tr in fam.template.st])
                    min_pick_time = min(
                        [pick.time for pick in det.event.picks])
                    for pick in det.event.picks:
                        traces = fam.template.st.select(
                            id=pick.waveform_id.get_seed_string())
                        lags = []
                        for tr in traces:
                            lags.append(
                                tr.stats.starttime - min_template_time)
                        assert pick.time - min_pick_time in lags
                    continue
                if isinstance(det.__dict__[key], float):
                    if not np.allclose(
                            det.__dict__[key], check_det.__dict__[key],
                            atol=float_tol):
                        print(key)
                    assert np.allclose(
                        det.__dict__[key], check_det.__dict__[key],
                        atol=float_tol)
                elif isinstance(det.__dict__[key], np.float32):
                    if not np.allclose(
                            det.__dict__[key], check_det.__dict__[key],
                            atol=float_tol):
                        print("{0}: new: {1}\tcheck-against: {2}".format(
                            key, det.__dict__[key], check_det.__dict__[key]))
                        print(det)
                    assert np.allclose(
                        det.__dict__[key], check_det.__dict__[key],
                        atol=float_tol)
                elif isinstance(det.__dict__[key], UTCDateTime):
                    if not det.__dict__[key] == check_det.__dict__[key]:
                        print("{0}: new: {1}\tcheck-against: {2}".format(
                            key, det.__dict__[key], check_det.__dict__[key]))
                        print(det)
                        print(check_det)
                    assert (abs(
                        det.__dict__[key] - check_det.__dict__[key]) <= 0.1)
                elif key in ['template_name', 'id']:
                    continue
                    # Name relies on creation-time, which is checked elsewhere,
                    # ignore it.
                else:
                    if not det.__dict__[key] == check_det.__dict__[key]:
                        print(key)
                    assert det.__dict__[key] == check_det.__dict__[key]


def test_match_filter(plot=False, extract_detections=False,
                      threshold_type='MAD', threshold=10,
                      template_excess=False, stream_excess=False,
                      conc_proc=False):
    """
    Function to test the capabilities of match_filter and just check that \
    it is working!  Uses synthetic templates and seeded, randomised data.

    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz to use
    """
    from eqcorrscan.utils import pre_processing
    from obspy import UTCDateTime
    import string
    import inspect
    # Read in the synthetic dataset
    templates = []
    testing_path = os.path.join(os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe()))), 'test_data',
        'synthetic_data')
    templates.append(read(os.path.join(testing_path, 'synth_template_0.ms')))
    templates.append(read(os.path.join(testing_path, 'synth_template_1.ms')))
    data = read(os.path.join(testing_path, 'synth_data.ms'))
    seeds = np.load(os.path.join(testing_path, 'seeds.npz'))
    seeds = [{'SNR': seeds['SNR_0'], 'time': seeds['time_0']},
             {'SNR': seeds['SNR_1'], 'time': seeds['time_1']}]
    if stream_excess:
        i = 0
        for template in templates:
            template = template.remove(template[i])
            if i < len(template):
                i += 1
            else:
                i = 0
    if template_excess:
        data = data[0:-1]
    # Filter the data and the templates
    for template in templates:
        for tr in template:
            tr.data += 1  # Make the synthetic data not be all zeros
        pre_processing.multi_process(
            st=template, lowcut=1.0, highcut=4.0, filt_order=3, samp_rate=10.0,
            seisan_chan_names=True)
    template_names = list(string.ascii_lowercase)[0:len(templates)]
    detections = match_filter(
        template_names=template_names, template_list=templates, st=data,
        threshold=threshold, threshold_type=threshold_type, trig_int=6.0,
        plot=plot, plotdir='.', cores=1, output_cat=False,
        extract_detections=extract_detections, concurrent_processing=conc_proc)
    if extract_detections:
        detection_streams = detections[1]
        detections = detections[0]
    # Compare the detections to the seeds
    print('This test made ' + str(len(detections)) + ' detections')
    ktrue = 0
    kfalse = 0
    for detection in detections:
        print(detection)
        i = template_names.index(detection.template_name)
        t_seeds = seeds[i]
        dtime_samples = int((detection.detect_time - UTCDateTime(0)) *
                            10.0)
        if dtime_samples in t_seeds['time']:
            j = list(t_seeds['time']).index(dtime_samples)
            print('Detection at SNR of: ' + str(t_seeds['SNR'][j]))
            ktrue += 1
        else:
            min_diff = min(abs(t_seeds['time'] - dtime_samples))
            if min_diff < 10:
                # If there is a match within ten samples then it is
                # good enough
                j = list(abs(t_seeds['time'] -
                             dtime_samples)).index(min_diff)
                print('Detection at SNR of: ' + str(t_seeds['SNR'][j]))
                ktrue += 1
            else:
                print('Detection at sample: ' + str(dtime_samples) +
                      ' does not match anything in seed times:')
                kfalse += 1
            print('Minimum difference in samples is: ' + str(min_diff))
    # print('Catalog created is of length: ' + str(len(out_cat)))
    # Set an 'acceptable' ratio of positive to false detections
    print(str(ktrue) + ' true detections and ' + str(kfalse) +
          ' false detections')
    if not extract_detections:
        return kfalse, ktrue
    else:
        return kfalse, ktrue, detection_streams


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
