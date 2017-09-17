"""
A series of test functions for the core functions in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import os
import unittest
import pytest

import numpy as np
from obspy import read, UTCDateTime, read_events, Catalog, Stream, Trace
from obspy.clients.fdsn import Client
from obspy.core.event import Pick, Event
from obspy.core.util.base import NamedTemporaryFile

from eqcorrscan.core.match_filter import MatchFilterError
from eqcorrscan.core.match_filter import match_filter, normxcorr2, Detection
from eqcorrscan.core.match_filter import read_detections, get_catalog
from eqcorrscan.core.match_filter import write_catalog, extract_from_stream
from eqcorrscan.core.match_filter import Tribe, Template, Party, Family
from eqcorrscan.core.match_filter import read_party, read_tribe, _spike_test
from eqcorrscan.tutorials.get_geonet_events import get_geonet_events
from eqcorrscan.utils import pre_processing, catalog_utils
from eqcorrscan.utils.correlate import fftw_normxcorr, numpy_normxcorr


class TestCoreMethods(unittest.TestCase):
    """
    Tests for internal _template_loop and normxcorr2 functions.
    """
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
        test_match_filter(stream_excess=True)

    def test_extra_templates(self):
        # Test case where there are non-matching streams in the data
        test_match_filter(template_excess=True)

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
        match_filter(template_names=['1'], template_list=templates, st=stream,
                     threshold=8, threshold_type='MAD', trig_int=1,
                     plotvar=False)


@pytest.mark.network
class TestGeoNetCase(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        client = Client('GEONET')
        cls.t1 = UTCDateTime(2016, 9, 4)
        cls.t2 = cls.t1 + 86400
        catalog = get_geonet_events(
            startdate=cls.t1, enddate=cls.t2, minmag=4, minlat=-49, maxlat=-35,
            minlon=175.0, maxlon=185.0)
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
        cls.st = pre_processing.shortproc(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=50.0,
            debug=0, num_cores=1)
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
        detections = match_filter(template_names=self.template_names,
                                  template_list=templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plotvar=False, plotdir='.',
                                  cores=1)
        self.assertEqual(len(detections), 1)
        self.assertEqual(detections[0].no_chans, 6)

    def test_duplicate_cont_data(self):
        """ Check that error is raised if duplicate channels are present in
        the continuous data."""
        tr = self.st[0].copy()
        tr.data = np.random.randn(100)
        st = self.st.copy() + tr
        with self.assertRaises(MatchFilterError):
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=st, threshold=8.0,
                         threshold_type='MAD', trig_int=6.0, plotvar=False,
                         plotdir='.', cores=1)

    def test_missing_cont_channel(self):
        """ Remove one channel from continuous data and check that everything
        still works. """
        st = self.st.copy()
        st.remove(st[-1])
        detections, det_cat = match_filter(
            template_names=self.template_names, template_list=self.templates,
            st=st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
            plotvar=False, plotdir='.', cores=1, output_cat=True)
        self.assertEqual(len(detections), 1)
        self.assertEqual(detections[0].no_chans, 5)
        self.assertEqual(len(detections), len(det_cat))

    def test_no_matching_data(self):
        """ No matching data between continuous and templates."""
        st = self.st.copy()
        for tr, staname in zip(st, ['a', 'b', 'c', 'd', 'e']):
            tr.stats.station = staname
        with self.assertRaises(IndexError):
            match_filter(
                template_names=self.template_names,
                template_list=self.templates, st=st, threshold=8.0,
                threshold_type='MAD', trig_int=6.0, plotvar=False,
                plotdir='.', cores=1)

    @pytest.mark.flaky(reruns=2)
    def test_geonet_tribe_detect(self):
        client = Client('GEONET')
        # Try to force issues with starting samples on wrong day for geonet
        # data
        tribe = self.tribe.copy()
        for template in tribe.templates:
            template.process_length = 86400
            template.st = Stream(template.st[0])
            # Only run one channel templates
        party = self.tribe.copy().client_detect(
            client=client, starttime=self.t1, endtime=self.t2,
            threshold=8.0, threshold_type='MAD', trig_int=6.0,
            daylong=False, plotvar=False)
        self.assertEqual(len(party), 16)


@pytest.mark.network
class TestNCEDCCases(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        print('\t\t\t Downloading data')
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
        print('Constructing tribe')
        cls.tribe.construct(
            method='from_client', catalog=catalog, client_id='NCEDC',
            lowcut=2.0, highcut=9.0, samp_rate=50.0, filt_order=4,
            length=3.0, prepick=0.15, swin='all', process_len=process_len)
        print('Constructing templates')
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
        print('Downloading continuous data')
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        cls.unproc_st = st.copy()
        print('Processing data')
        cls.st = pre_processing.shortproc(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=50.0,
            debug=0, num_cores=1, starttime=st[0].stats.starttime,
            endtime=st[0].stats.starttime + process_len)
        cls.template_names = [str(template[0].stats.starttime)
                              for template in cls.templates]
        print('Set Up finished')

    def test_detection_extraction(self):
        # Test outputting the streams works
        detections, detection_streams = \
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=self.st,
                         threshold=8.0, threshold_type='MAD',
                         trig_int=6.0, plotvar=False, plotdir='.',
                         cores=1, extract_detections=True)
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
        self.assertTrue(np.allclose(ccc, ccc_numpy, atol=0.03))

    def test_catalog_extraction(self):
        detections, det_cat, detection_streams = \
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=self.st,
                         threshold=8.0, threshold_type='MAD',
                         trig_int=6.0, plotvar=False, plotdir='.',
                         cores=1, extract_detections=True, output_cat=True)
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
                trig_int=6.0, plotvar=False, plotdir='.', cores=1)
        individual_dict = []
        for detection in individual_detections:
            individual_dict.append({'template_name': detection.template_name,
                                    'time': detection.detect_time,
                                    'cccsum': detection.detect_val.round(6)})
        detections = match_filter(template_names=self.template_names,
                                  template_list=self.templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plotvar=False, plotdir='.',
                                  cores=1)
        self.assertEqual(len(individual_detections), len(detections))
        for detection in detections:
            detection_dict = {'template_name': detection.template_name,
                              'time': detection.detect_time,
                              'cccsum': detection.detect_val.round(6)}
            self.assertTrue(detection_dict in individual_dict)

    def test_read_write_detections(self):
        """Check that we can read and write detections accurately."""
        detections = match_filter(template_names=self.template_names,
                                  template_list=self.templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plotvar=False, plotdir='.',
                                  cores=1)
        detection_dict = []
        for detection in detections:
            detection_dict.append(
                {'template_name': detection.template_name,
                 'time': detection.detect_time,
                 'cccsum': float(str(detection.detect_val)),
                 'no_chans': detection.no_chans,
                 'chans': detection.chans,
                 'threshold': float(str(detection.threshold)),
                 'typeofdet': detection.typeofdet})
            detection.write(fname='dets_out.txt')
        detections_back = read_detections('dets_out.txt')
        print(detection_dict)
        self.assertEqual(len(detections), len(detections_back))
        for detection in detections_back:
            d = {'template_name': detection.template_name,
                 'time': detection.detect_time,
                 'cccsum': detection.detect_val,
                 'no_chans': detection.no_chans,
                 'chans': detection.chans,
                 'threshold': detection.threshold,
                 'typeofdet': detection.typeofdet}
            print(d)
            self.assertTrue(d in detection_dict)
        os.remove('dets_out.txt')

    def test_get_catalog(self):
        """Check that catalog objects are created properly."""
        detections = match_filter(
            template_names=self.template_names, template_list=self.templates,
            st=self.st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
            plotvar=False, plotdir='.', cores=1)
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
                                  trig_int=6.0, plotvar=False, plotdir='.',
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
                         plotvar=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # templates is not a list
            match_filter(template_names=self.template_names,
                         template_list=self.templates[0], st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plotvar=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # template and template_names length are not equal
            match_filter(template_names=self.template_names,
                         template_list=[self.templates[0]], st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plotvar=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # templates is not a list of streams
            match_filter(template_names=self.template_names,
                         template_list=['abc'], st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plotvar=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # st is not a Stream
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=np.random.randn(10),
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plotvar=False, plotdir='.', cores=1)
        with self.assertRaises(MatchFilterError):
            # threshold_type is wrong
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=self.st,
                         threshold=8.0, threshold_type='albert', trig_int=6.0,
                         plotvar=False, plotdir='.', cores=1)

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
                         plotvar=False, plotdir='.', cores=1)

    def test_non_equal_template_lengths(self):
        templates = [self.templates[0].copy()]
        templates[0][0].data = np.concatenate([templates[0][0].data,
                                               np.random.randn(10)])
        with self.assertRaises(MatchFilterError):
            match_filter(template_names=[self.template_names[0]],
                         template_list=templates, st=self.st,
                         threshold=8.0, threshold_type='MAD', trig_int=6.0,
                         plotvar=False, plotdir='.', cores=1)


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
class TestMatchObjects(unittest.TestCase):
    @classmethod
    @pytest.mark.flaky(reruns=2)
    def setUpClass(cls):
        print('\t\t\t Downloading data')
        client = Client('NCEDC')
        cls.t1 = UTCDateTime(2004, 9, 28, 17)
        cls.t2 = cls.t1 + 3600
        process_len = 3600
        # t1 = UTCDateTime(2004, 9, 28)
        # t2 = t1 + 80000
        # process_len = 80000
        catalog = client.get_events(
            starttime=cls.t1, endtime=cls.t2, minmagnitude=4,
            minlatitude=35.7, maxlatitude=36.1, minlongitude=-120.6,
            maxlongitude=-120.2, includearrivals=True)
        catalog = catalog_utils.filter_picks(
            catalog, channels=['EHZ'], top_n_picks=5)
        cls.tribe = Tribe()
        print('Constructing tribe')
        cls.tribe.construct(
            method='from_client', catalog=catalog, client_id='NCEDC',
            lowcut=2.0, highcut=9.0, samp_rate=20.0, filt_order=4,
            length=3.0, prepick=0.15, swin='all', process_len=process_len)
        cls.onehztribe = Tribe().construct(
            method='from_client', catalog=catalog, client_id='NCEDC',
            lowcut=0.1, highcut=0.45, samp_rate=1.0, filt_order=4,
            length=20.0, prepick=0.15, swin='all', process_len=process_len)
        # Download and process the day-long data
        template_stachans = []
        for template in cls.tribe.templates:
            for tr in template.st:
                template_stachans.append(
                    (tr.stats.network, tr.stats.station, tr.stats.channel))
        cls.template_stachans = list(set(template_stachans))
        bulk_info = [(stachan[0], stachan[1], '*', stachan[2],
                      cls.t1, cls.t1 + process_len)
                     for stachan in template_stachans]
        # Just downloading an hour of data
        print('Downloading continuous data')
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        cls.unproc_st = st.copy()
        print('Processing data')
        cls.st = pre_processing.shortproc(
            st, lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=20.0,
            debug=0, num_cores=1, starttime=st[0].stats.starttime,
            endtime=st[0].stats.starttime + process_len)
        print('Reading party')
        cls.party = Party().read(
            filename=os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                'test_data', 'test_party.tgz'))
        cls.family = cls.party.sort()[0].copy()
        print('Set Up finished')

    def test_tribe_internal_methods(self):
        self.assertEqual(len(self.tribe), 4)
        self.assertTrue(self.tribe == self.tribe)
        self.assertFalse(self.tribe != self.tribe)

    def test_tribe_add(self):
        """Test add method"""
        added = self.tribe.copy()
        self.assertEqual(len(added + added[0]), 5)
        added += added[-1]
        self.assertEqual(len(added), 6)

    def test_tribe_remove(self):
        """Test remove method"""
        removal = self.tribe.copy()
        self.assertEqual(len(removal.remove(removal[0])), 3)
        for template in self.tribe:
            self.assertTrue(isinstance(template, Template))

    def test_tribe_io(self):
        """Test reading and writing or Tribe objects using tar form."""
        try:
            if os.path.isfile('test_tribe.tgz'):
                os.remove('test_tribe.tgz')
            self.tribe.write(filename='test_tribe')
            tribe_back = read_tribe('test_tribe.tgz')
            self.assertEqual(self.tribe, tribe_back)
        finally:
            os.remove('test_tribe.tgz')

    def test_tribe_detect(self):
        """Test the detect method on Tribe objects"""
        party = self.tribe.detect(
            stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
            trig_int=6.0, daylong=False, plotvar=False, parallel_process=False)
        self.assertEqual(len(party), 4)
        for fam, check_fam in zip(party, self.party):
            for det, check_det in zip(fam.detections, check_fam.detections):
                for key in det.__dict__.keys():
                    if key == 'event':
                        continue
                    if isinstance(det.__dict__[key], float):
                        if not np.allclose(det.__dict__[key],
                                           check_det.__dict__[key], atol=0.1):
                            print(key)
                        self.assertTrue(np.allclose(
                            det.__dict__[key], check_det.__dict__[key],
                            atol=0.2))
                    else:
                        self.assertEqual(
                            det.__dict__[key], check_det.__dict__[key])
            # self.assertEqual(fam.template, check_fam.template)

    @pytest.mark.flaky(reruns=2)
    def test_client_detect(self):
        """Test the client_detect method."""
        client = Client('NCEDC')
        party = self.tribe.copy().client_detect(
            client=client, starttime=self.t1, endtime=self.t2,
            threshold=8.0, threshold_type='MAD', trig_int=6.0,
            daylong=False, plotvar=False)
        self.assertEqual(len(party), 4)

    def test_party_io(self):
        """Test reading and writing party objects."""
        if os.path.isfile('test_party_out.tgz'):
            os.remove('test_party_out.tgz')
        try:
            self.party.write(filename='test_party_out')
            party_back = read_party(fname='test_party_out.tgz')
            self.assertEqual(self.party, party_back)
        finally:
            os.remove('test_party_out.tgz')

    def test_party_basic_methods(self):
        """Test the basic methods on Party objects."""
        self.assertEqual(self.party.__repr__(), 'Party of 4 Families.')
        self.assertFalse(self.party != self.party)

    def test_party_add(self):
        """Test getting items and adding them to party objects, and sorting"""
        test_party = self.party.copy()
        test_family = test_party[0]
        self.assertTrue(isinstance(test_family, Family))
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
        for trig_int in [40, 15, 3600]:
            for metric in ['avg_cor', 'cor_sum']:
                declust = self.party.copy().decluster(trig_int=trig_int,
                                                      metric=metric)
                declustered_dets = [d for family in declust
                                    for d in family.detections]
                for det in declustered_dets:
                    time_difs = [abs(det.detect_time - d.detect_time)
                                 for d in declustered_dets]
                    for dif in time_difs:
                        if dif != 0:
                            self.assertTrue(dif > trig_int)
                with self.assertRaises(MatchFilterError):
                    self.party.copy().decluster(trig_int=trig_int,
                                                timing='origin',
                                                metric=metric)

    def test_party_lag_calc(self):
        """Test the lag-calc method on Party objects."""
        # Test the chained method
        chained_cat = self.tribe.detect(
            stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
            trig_int=6.0, daylong=False, plotvar=False).lag_calc(
            stream=self.unproc_st, pre_processed=False)
        catalog = self.party.lag_calc(stream=self.unproc_st,
                                      pre_processed=False)
        self.assertEqual(len(catalog), 3)
        # Check that the party is unaltered
        self.assertEqual(self.party, read_party(
            fname=os.path.join(os.path.abspath(os.path.dirname(__file__)),
                               'test_data', 'test_party.tgz')))
        for ev, chained_ev in zip(catalog, chained_cat):
            for i in range(len(ev.picks)):
                for key in ev.picks[i].keys():
                    if key == 'resource_id':
                        continue
                    if key == 'comments':
                        self.assertEqual(
                            sorted(ev.picks,
                                   key=lambda p: p.time)[i][key][0].text,
                            sorted(chained_ev.picks,
                                   key=lambda p: p.time)[i][key][0].text)
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
                self.assertEqual(ev.resource_id, chained_ev.resource_id)
                self.assertEqual(ev.comments[0].text,
                                 chained_ev.comments[0].text)

    def test_party_lag_calc_preprocessed(self):
        """Test that the lag-calc works on pre-processed data."""
        catalog = self.party.lag_calc(stream=self.st, pre_processed=True)
        self.assertEqual(len(catalog), 3)

    def test_day_long_methods(self):
        """Conduct a test using day-long data."""
        client = Client('NCEDC')
        t1 = UTCDateTime(2004, 9, 28)
        bulk_info = [(stachan[0], stachan[1], '*', stachan[2], t1, t1 + 86400)
                     for stachan in self.template_stachans]
        # Just downloading an hour of data
        print('Downloading continuous day-long data')
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        # Hack day-long templates
        daylong_tribe = self.onehztribe.copy()
        for template in daylong_tribe:
            template.process_length = 86400
        # Aftershock sequence, with 1Hz data, lots of good correlations = high
        # MAD!
        day_party = daylong_tribe.detect(
            stream=st, threshold=8.0, threshold_type='MAD', trig_int=6.0,
            daylong=True, plotvar=False, parallel_process=False)
        self.assertEqual(len(day_party), 4)
        day_catalog = day_party.lag_calc(stream=st, pre_processed=False,
                                         parallel=False)
        self.assertEqual(len(day_catalog), 3)
        pre_picked_cat = day_party.get_catalog()
        self.assertEqual(len(pre_picked_cat), 4)

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
            family.write('test_family')
            party_back = read_party('test_family.tgz')
            self.assertEqual(len(party_back), 1)
            self.assertEqual(party_back[0], family)
        finally:
            os.remove('test_family.tgz')

    def test_family_lag_calc(self):
        """Test the lag-calc method on family."""
        catalog = self.family.lag_calc(stream=self.st, pre_processed=True)
        self.assertEqual(len(catalog), 1)

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

    def test_template_io(self):
        """Test template read/write."""
        test_template = self.family.template.copy()
        try:
            test_template.write('test_template')
            template_back = Template().read('test_template.tgz')
            self.assertEqual(test_template, template_back)
        finally:
            os.remove('test_template.tgz')
        # Make sure we raise a useful error when trying to read from a
        # tribe file.
        try:
            self.tribe.write(filename='test_template')
            with self.assertRaises(IOError):
                Template().read('test_template.tgz')
        finally:
            os.remove('test_template.tgz')

    def test_template_detect(self):
        """Test detect method on Template objects."""
        test_template = self.family.template.copy()
        party_t = test_template.detect(
            stream=self.unproc_st, threshold=8.0, threshold_type='MAD',
            trig_int=6.0, daylong=False, plotvar=False, overlap=None)
        self.assertEqual(len(party_t), 1)

    def test_template_construct(self):
        """Test template construction."""
        test_template = Template()
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_',
                                    '15-0931-08L.S201309')
        test_template.construct(
            method='from_sfile', name='tester', lowcut=2, highcut=8,
            samp_rate=20, filt_order=4, length=10, swin='all', prepick=0.2,
            debug=3, sfile=testing_path)
        self.assertTrue(isinstance(test_template, Template))
        self.assertEqual(test_template.name, 'tester')
        with self.assertRaises(NotImplementedError):
            test_template.construct(
                method='from_client', client_id='NCEDC', name='bob',
                lowcut=2.0, highcut=9.0, samp_rate=50.0, filt_order=4,
                length=3.0, prepick=0.15, swin='all', process_len=6)

    def test_slicing_by_name(self):
        """Check that slicing by template name works as expected"""
        t_name = self.tribe[2].name
        self.assertTrue(self.tribe[2] == self.tribe[t_name])
        self.assertTrue(self.party[2] == self.party[t_name])


def test_match_filter(
        debug=0, plotvar=False, extract_detections=False, threshold_type='MAD',
        threshold=10, template_excess=False, stream_excess=False):
    """
    Function to test the capabilities of match_filter and just check that \
    it is working!  Uses synthetic templates and seeded, randomised data.

    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz to use
    :type debug: int
    :param debug: Debug level, higher the number the more output.
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
        pre_processing.shortproc(
            st=template, lowcut=1.0, highcut=4.0, filt_order=3, samp_rate=10.0,
            seisan_chan_names=True)
    template_names = list(string.ascii_lowercase)[0:len(templates)]
    detections = match_filter(
        template_names=template_names, template_list=templates, st=data,
        threshold=threshold, threshold_type=threshold_type, trig_int=6.0,
        plotvar=plotvar, plotdir='.', cores=1, debug=debug, output_cat=False,
        extract_detections=extract_detections)
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
    # Plot the detections
    if debug > 3:
        for i, template in enumerate(templates):
            times = [d.detect_time.datetime for d in detections
                     if d.template_name == template_names[i]]
            print(times)
            # plotting.detection_multiplot(data, template, times)
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
