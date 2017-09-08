"""
Functions to test generating templates from SAC data.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import pytest
import glob
import os
import numpy as np
import warnings
import shutil
import inspect
import copy

from obspy import read, UTCDateTime, read_events, Stream
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog, Event, Origin, Pick, WaveformStreamID

from eqcorrscan.core.template_gen import from_sac, _group_events, from_seishub
from eqcorrscan.core.template_gen import from_meta_file, from_client
from eqcorrscan.core.template_gen import multi_template_gen, from_contbase
from eqcorrscan.core.template_gen import template_gen, extract_from_stack
from eqcorrscan.core.template_gen import from_sfile, TemplateGenError
from eqcorrscan.tutorials.template_creation import mktemplates
from eqcorrscan.tutorials.get_geonet_events import get_geonet_events
from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.utils.sfile_util import eventtosfile, read_event


class TestTemplateGeneration(unittest.TestCase):
    """Test the reading a writing of pick info."""
    def test_sac_template_gen(self):
        """Test template generation."""
        samp_rate = 20
        length = 8

        for event in ['2014p611252', 'No_head']:
            test_files = os.path.join(os.path.abspath(os.path.
                                                      dirname(__file__)),
                                      'test_data', 'SAC', event, '*')
            # Test with various input types
            filelist = glob.glob(test_files)
            streamlist = [read(f) for f in glob.glob(test_files)]
            stream = read(test_files)
            for sac_files in [filelist, streamlist, stream]:
                template = from_sac(sac_files, lowcut=2.0, highcut=8.0,
                                    samp_rate=samp_rate, filt_order=4,
                                    length=length, swin='all', prepick=0.1,
                                    debug=0, plot=False)
                for tr in template:
                    self.assertEqual(len(tr.data), length * samp_rate)

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_tutorial_template_gen(self):
        """Test template generation from tutorial, uses from_client method.

        Checks that the tutorial generates the templates we expect it to!
        """
        testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        mktemplates(plot=False)
        for template_no in range(4):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            expected_template = read(os.path.join(
                testing_path, 'tutorial_template_' + str(template_no) + '.ms'))
            for tr in template:
                expected_tr = expected_template.select(
                    station=tr.stats.station, channel=tr.stats.channel)[0]
                self.assertTrue((expected_tr.data.astype(np.float32) ==
                                 tr.data.astype(np.float32)).all())
            del(template)
            os.remove('tutorial_template_' + str(template_no) + '.ms')

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_not_delayed(self):
        """Test the method of template_gen without applying delays to
        channels."""
        cat = get_geonet_events(
            minlat=-40.98, maxlat=-40.85, minlon=175.4, maxlon=175.5,
            startdate=UTCDateTime(2016, 5, 1), enddate=UTCDateTime(2016, 5, 2))
        cat = filter_picks(catalog=cat, top_n_picks=5)
        template = from_client(
            catalog=cat, client_id='GEONET', lowcut=None, highcut=None,
            samp_rate=100.0, filt_order=4, length=10.0, prepick=0.5,
            swin='all', process_len=3600, debug=0, plot=False,
            delayed=False)[0]
        for tr in template:
            tr.stats.starttime.precision = 6
        starttime = template[0].stats.starttime
        length = template[0].stats.npts
        print(template)
        for tr in template:
            self.assertTrue(abs((tr.stats.starttime - starttime)) <=
                            tr.stats.delta)
            self.assertEqual(tr.stats.npts, length)

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_download_various_methods(self):
        """
        Will download data from server and store in various databases,
        then create templates using the various methods.
        """
        client = Client('GEONET')
        # get the events
        catalog = Catalog()
        data_stream = client._download(
            'http://quakeml.geonet.org.nz/quakeml/1.2/2016p008194')
        data_stream.seek(0, 0)
        catalog += read_events(data_stream, format="quakeml")
        data_stream.close()
        # Select 3 channels to use and download
        sta_chans = [(pick.waveform_id.station_code,
                      pick.waveform_id.channel_code)
                     for pick in catalog[0].picks[0:3]]
        t1 = UTCDateTime(catalog[0].origins[0].time.date)
        t2 = t1 + 86400
        bulk = [('NZ', sta_chan[0], '*', sta_chan[1], t1, t2)
                for sta_chan in sta_chans]
        continuous_st = client.get_waveforms_bulk(bulk)
        continuous_st.merge(fill_value=0)
        # Test multi_template_gen
        templates = multi_template_gen(catalog, continuous_st, length=3)
        self.assertEqual(len(templates), 1)
        # Test without an event
        templates = multi_template_gen(Catalog(), continuous_st, length=3)
        self.assertEqual(len(templates), 0)
        # Test from contbase method
        sfile = eventtosfile(catalog[0], 'TEST', 'L', '.', 'None',
                             overwrite=True)
        os.makedirs(catalog[0].origins[0].time.date.strftime('Y%Y'))
        os.makedirs(catalog[0].origins[0].time.date.
                    strftime('Y%Y' + os.sep + 'R%j.01'))
        for tr in continuous_st:
            tr.write(catalog[0].origins[0].time.date.
                     strftime('Y%Y' + os.sep + 'R%j.01') + os.sep +
                     tr.stats.station + '.' + tr.stats.network + '.' +
                     tr.stats.location + '.' + tr.stats.channel +
                     tr.stats.starttime.strftime('%Y.%j'), format='MSEED')
        template = from_contbase(sfile,
                                 contbase_list=[('.', 'Yyyyy/Rjjj.01', 'NZ')],
                                 lowcut=1.0, highcut=5.0, samp_rate=20,
                                 filt_order=4, length=3, prepick=0.5,
                                 swin='all')
        self.assertTrue(isinstance(template, Stream))
        shutil.rmtree(continuous_st[0].stats.starttime.strftime('Y%Y'))

    def test_seishub(self):
        """Test the seishub method, use obspy default seishub client."""
        import sys
        if sys.version_info.major == 2:
            from future.backports.urllib.request import URLError
        else:
            from urllib.request import URLError
        t = UTCDateTime(2009, 9, 3)
        test_cat = Catalog()
        test_cat.append(Event())
        test_cat[0].origins.append(Origin())
        test_cat[0].origins[0].time = t
        test_cat[0].origins[0].latitude = 45
        test_cat[0].origins[0].longitude = 45
        test_cat[0].origins[0].depth = 5000
        test_cat[0].picks.append(Pick(
            waveform_id=WaveformStreamID(
                station_code='MANZ', channel_code='EHZ', network_code='BW'),
            phase_hint='PG', time=t + 2000))
        test_cat[0].picks.append(Pick(
            waveform_id=WaveformStreamID(
                station_code='MANZ', channel_code='EHN', network_code='BW'),
            phase_hint='SG', time=t + 2005))
        test_cat[0].picks.append(Pick(
            waveform_id=WaveformStreamID(
                station_code='MANZ', channel_code='EHE', network_code='BW'),
            phase_hint='SG', time=t + 2005.5))

        test_url = "http://teide.geophysik.uni-muenchen.de:8080"

        if sys.version_info.major == 3:
            try:
                template = from_seishub(
                    test_cat, url=test_url, lowcut=1.0, highcut=5.0,
                    samp_rate=20, filt_order=4, length=3, prepick=0.5,
                    swin='all', process_len=300)
            except URLError:
                warnings.warn('Timed out connection to seishub')
        else:
            warnings.warn('URLError would not be caught on py2.')
        if 'template' in locals():
            self.assertEqual(len(template), 3)

    def test_catalog_grouping(self):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_', '*.S??????')
        catalog = Catalog()
        sfiles = glob.glob(testing_path)
        for sfile in sfiles:
            catalog.append(read_event(sfile=sfile))
        for process_len, pads in [(60, [5]),
                                  (300, [5, 60]),
                                  (3600, [5, 60, 300]),
                                  (86400, [5, 60, 300])]:
            for data_pad in pads:
                sub_catalogs = _group_events(catalog=catalog,
                                             process_len=process_len,
                                             data_pad=data_pad)
                k_events = 0
                for sub_catalog in sub_catalogs:
                    min_time = min([event.origins[0].time
                                    for event in sub_catalog])
                    min_time -= data_pad
                    for event in sub_catalog:
                        self.assertTrue((event.origins[0].time +
                                         data_pad) - min_time < process_len)
                        k_events += 1
                self.assertEqual(k_events, len(catalog))

    def test_missing_waveform_id(self):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        quakeml = os.path.join(testing_path,
                               '20130901T041115_missingwavid.xml')
        st = read(os.path.join(testing_path, 'WAV', 'TEST_',
                               '2013-09-01-0410-35.DFDPC_024_00'))
        templates = from_meta_file(meta_file=quakeml, st=st, lowcut=2.0,
                                   highcut=9.0, samp_rate=20.0, filt_order=3,
                                   length=2, prepick=0.1, swin='S')
        self.assertEqual(len(templates), 1)

    def test_from_sfile(self):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_',
                                    '15-0931-08L.S201309')
        event = read_event(sfile=testing_path)
        template = from_sfile(sfile=testing_path, lowcut=2, highcut=8,
                              samp_rate=20, filt_order=4, length=10,
                              swin='all', prepick=0.2, debug=3)
        self.assertEqual(len(template), len(event.picks))

    def test_upsample_error(self):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_',
                                    '15-0931-08L.S201309')
        with self.assertRaises(TemplateGenError):
            from_sfile(sfile=testing_path, lowcut=2, highcut=8, samp_rate=200,
                       filt_order=4, length=10, swin='all', prepick=0.2)


class TestEdgeGen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testing_path = os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe())))
        cls.st = read(os.path.join(
            cls.testing_path, 'test_data', 'WAV', 'TEST_',
            '2013-09-15-0930-28.DFDPC_027_00'))
        for tr in cls.st:
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        event = read_event(os.path.join(
            cls.testing_path, 'test_data', 'REA', 'TEST_',
            '15-0931-08L.S201309'))
        cls.picks = event.picks

    def test_undefined_phase_type(self):
        with self.assertRaises(IOError):
            template_gen(
                picks=self.picks, st=self.st.copy(), length=2, swin='bob')

    def test_warn_zeros(self):
        st = self.st.copy()
        template = template_gen(self.picks, st.copy(), 10)
        self.assertTrue('LABE' in [tr.stats.station for tr in template])
        st.select(station='LABE', channel='SN')[0].data = np.zeros(10000)
        template = template_gen(self.picks, st, 10)
        self.assertFalse('LABE' in [tr.stats.station for tr in template])

    def test_missing_data(self):
        picks = copy.deepcopy(self.picks)
        picks.append(picks[-1])
        picks[-1].waveform_id.station_code = 'DUMMY'
        template = template_gen(picks, self.st.copy(), 10)
        self.assertFalse('DUMMY' in [tr.stats.station for tr in template])

    def test_no_matched_picks(self):
        picks = [copy.deepcopy(self.picks[0])]
        picks[0].waveform_id.station_code = 'DUMMY'
        template = template_gen(picks, self.st.copy(), 10)
        self.assertFalse(template)

    def test_debug_levels(self):
        template = template_gen(self.picks, self.st.copy(), 10, debug=3)
        self.assertEqual(len(template), len(self.picks))

    def test_extract_from_stack(self):
        length = 3
        stack = self.st.copy()
        template = template_gen(self.picks, self.st.copy(), 2)
        extracted = extract_from_stack(stack, template, length=length,
                                       pre_pick=0.3, pre_pad=45)
        self.assertEqual(len(template), len(extracted))
        for tr in extracted:
            self.assertEqual(tr.stats.endtime - tr.stats.starttime,
                             length)

    def test_extract_from_stack_and_process(self):
        length = 3
        stack = self.st.copy()
        template = template_gen(self.picks, self.st.copy(), 2)
        extracted = extract_from_stack(
            stack, template, length=length, pre_pick=0.3, pre_pad=45,
            pre_processed=False, samp_rate=20, lowcut=2, highcut=8)
        self.assertEqual(len(template), len(extracted))
        for tr in extracted:
            self.assertEqual(tr.stats.endtime - tr.stats.starttime,
                             length)

    def test_extract_from_stack_including_z(self):
        length = 3
        stack = self.st.copy()
        template = template_gen(self.picks, self.st.copy(), 2)
        extracted = extract_from_stack(
            stack, template, length=length, pre_pick=0.3, pre_pad=45,
            Z_include=True)
        self.assertLess(len(template), len(extracted))
        for tr in extracted:
            self.assertEqual(tr.stats.endtime - tr.stats.starttime,
                             length)

if __name__ == '__main__':
    unittest.main()
