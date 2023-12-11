"""
Functions to test generating templates from SAC data.
"""
import unittest
import pytest
import glob
import os
import numpy as np
import logging
import inspect
import copy
import shutil

from obspy import read, UTCDateTime, read_events
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog, Event, Origin, Pick, WaveformStreamID

from eqcorrscan.core import template_gen as template_gen_module
from eqcorrscan.core.template_gen import (
    _group_events, extract_from_stack, _template_gen, template_gen,
    TemplateGenError)
from eqcorrscan.tutorials.template_creation import mktemplates
from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.utils.sac_util import sactoevent
from eqcorrscan.helpers.mock_logger import MockLoggingHandler


class _StreamTestClient:
    """A simple waveform client using a stream."""

    def __init__(self, st=None):
        self.st = st or read()

    def get_waveforms(self, network, station, location, channel,
                      starttime, endtime):
        """Get waveforms contained by stream."""
        out = self.st.copy()
        st_trimed = out.trim(starttime=starttime, endtime=endtime)
        st_filtered = st_trimed.select(network=network, station=station,
                                       location=location, channel=channel)
        return st_filtered

    def get_default_catalog(self):
        """Get a catalog with picks from the default stream."""
        pick1 = Pick(
            time=UTCDateTime(2009, 8, 24, 0, 20, 7, 696381),
            waveform_id=WaveformStreamID(seed_string='BW.RJOB..EHZ'),
            phase_hint='P'
        )
        origin = Origin(
            time=UTCDateTime(2009, 8, 24, 0, 20, 6, 410034),
            longitude=0,
            latitude=0,
            depth=0,
        )
        event = Event(picks=[pick1], origins=[origin])
        return Catalog(events=[event])


class TestTemplateGeneration(unittest.TestCase):
    """Test the reading a writing of pick info."""
    def test_sac_template_gen(self):
        """Test template generation."""
        samp_rate = 20
        length = 8

        for event in ['2014p611252', 'No_head']:
            test_files = os.path.join(os.path.abspath(
                os.path.dirname(__file__)), 'test_data', 'SAC', event, '*')
            # Test with various input types
            filelist = glob.glob(test_files)
            streamlist = [read(f) for f in glob.glob(test_files)]
            stream = read(test_files)
            for sac_files in [filelist, streamlist, stream]:
                templates = template_gen(
                    method="from_sac", sac_files=sac_files, lowcut=2.0,
                    highcut=8.0, samp_rate=samp_rate,
                    filt_order=4, length=length, swin='all', prepick=0.1,
                    plot=False)
                self.assertEqual(len(templates), 1)
                template = templates[0]
                self.assertEqual(len(template), len(sactoevent(stream).picks))
                for tr in template:
                    self.assertEqual(len(tr.data), length * samp_rate)

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_tutorial_template_gen(self):
        """
        Test template generation from tutorial, uses from_client method.

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
                self.assertTrue(np.allclose(
                    expected_tr.data.astype(np.float32),
                    tr.data.astype(np.float32), rtol=0.0001))
            del template
            os.remove('tutorial_template_' + str(template_no) + '.ms')

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_not_delayed(self):
        """Test the method of template_gen without applying delays to
        channels."""
        client = Client("GEONET")
        cat = client.get_events(
            minlatitude=-40.98, maxlatitude=-40.85, minlongitude=175.4,
            maxlongitude=175.5, starttime=UTCDateTime(2016, 5, 1),
            endtime=UTCDateTime(2016, 5, 2))
        cat = filter_picks(catalog=cat, top_n_picks=5)
        template = template_gen(
            method="from_client", catalog=cat, client_id='GEONET',
            lowcut=None, highcut=None, samp_rate=100.0, filt_order=4,
            length=10.0, prepick=0.5, swin='all', process_len=3600,
            plot=False, delayed=False)[0]
        for tr in template:
            tr.stats.starttime.precision = 6
        starttime = template[0].stats.starttime
        length = template[0].stats.npts
        self.assertEqual(len(template), 5)
        for tr in template:
            self.assertTrue(abs((tr.stats.starttime - starttime)) <=
                            tr.stats.delta)
            self.assertEqual(tr.stats.npts, length)
        template = template_gen(
            method="from_client", catalog=cat, client_id='GEONET',
            lowcut=None, highcut=None, samp_rate=100.0, filt_order=4,
            length=10.0, prepick=0.5, swin='P_all', process_len=3600,
            plot=False, delayed=False)[0]
        for tr in template:
            tr.stats.starttime.precision = 6
        starttime = template[0].stats.starttime
        length = template[0].stats.npts
        self.assertEqual(len(template), 15)
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
        catalog = client.get_events(eventid='2016p008194')
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
        kwargs = {"process": False, "lowcut": None, "highcut": None,
                  "filt_order": None, "swin": "all", "prepick": 0.05,
                  "all_horiz": False, "delayed": True, "plot": False,
                  "return_event": False, "min_snr": None,
                  "samp_rate": continuous_st[0].stats.sampling_rate}
        templates = template_gen(
            method="from_meta_file", meta_file=catalog, st=continuous_st,
            length=3, **kwargs)
        self.assertEqual(len(templates), 1)
        # Test without an event
        templates = template_gen(
            method="from_meta_file", meta_file=Catalog(),
            st=continuous_st, length=3, **kwargs)
        self.assertEqual(len(templates), 0)

    @pytest.mark.network
    def test_save_progress(self):
        """ Test template creation with progress saving """
        client = Client('GEONET')
        catalog = client.get_events(
            starttime=UTCDateTime(2016, 1, 4, 0, 50),
            endtime=UTCDateTime(2016, 1, 4, 1, 20))
        # Gets a catalog of 2 events separated by 127s
        # Need a bigger gap to allow moveouts
        catalog[0].origins[0].time -= 600
        for pick in catalog[0].picks:
            pick.time -= 600
        catalog = filter_picks(catalog=catalog, top_n_picks=5)
        templates = template_gen(
            method="from_client", catalog=catalog, client_id="GEONET",
            lowcut=2, highcut=5, samp_rate=20, filt_order=4, prepick=0.4,
            process_len=600, swin="P", save_progress=True, length=2)
        assert(os.path.isdir("eqcorrscan_temporary_templates"))
        saved_templates = [
            read(f) for f in sorted(
                glob.glob("eqcorrscan_temporary_templates/*.ms"))]
        # Writing to miniseed adds miniseed stats dict
        for saved_template, template in zip(saved_templates, templates):
            for saved_tr in saved_template:
                tr = template.select(id=saved_tr.id)[0]
                assert(np.allclose(saved_tr.data, tr.data, atol=0.01))
        shutil.rmtree("eqcorrscan_temporary_templates")

    def test_all_phase_methods(self):
        sfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'test_data', 'REA', 'TEST_',
                             '01-0411-15L.S201309')
        catalog = read_events(sfile)
        p_stations = list(set(
            [pick.waveform_id.station_code
             for pick in catalog[0].picks if pick.phase_hint == 'P']))
        s_stations = list(set(
            [pick.waveform_id.station_code
             for pick in catalog[0].picks if pick.phase_hint == 'S']))
        st = read(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'test_data', 'WAV', 'TEST_',
                               '2013-09-01-0410-35.DFDPC_024_00'))
        templates = template_gen(
            method="from_meta_file", meta_file=sfile, st=st, lowcut=2,
            highcut=20, samp_rate=100, filt_order=4, length=6, prepick=0.2,
            swin='P_all')
        self.assertEqual(len(templates), 1)
        self.assertEqual(len(templates[0]), len(p_stations) * 3)
        for tr in templates[0]:
            pick = [
                p for p in catalog[0].picks
                if p.waveform_id.station_code == tr.stats.station and
                p.phase_hint.upper() == 'P'][0]
            print(tr)
            print(pick)
            self.assertLess(abs(tr.stats.starttime - (pick.time - 0.2)),
                            tr.stats.delta)
        templates = template_gen(
            method="from_meta_file", meta_file=sfile, st=st, lowcut=2,
            highcut=20, samp_rate=100, filt_order=4, length=6, prepick=0.2,
            swin='S_all')
        self.assertEqual(len(templates), 1)
        self.assertEqual(len(templates[0]), len(s_stations) * 3)
        for tr in templates[0]:
            pick = [
                p for p in catalog[0].picks
                if p.waveform_id.station_code == tr.stats.station and
                p.phase_hint.upper() == 'S'][0]
            print(tr)
            print(pick)
            self.assertLess(abs(tr.stats.starttime - (pick.time - 0.2)),
                            tr.stats.delta)

    def test_catalog_grouping(self):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_', '*.S??????')
        catalog = Catalog()
        sfiles = glob.glob(testing_path)
        for sfile in sfiles:
            catalog += read_events(sfile)
        for process_len, pads in [
           (60, [5]), (300, [5, 60]), (3600, [5, 60, 300]),
           (86400, [5, 60, 300])]:
            for data_pad in pads:
                sub_catalogs = _group_events(
                    catalog=catalog, process_len=process_len,
                    template_length=10, data_pad=data_pad)
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
        templates = template_gen(
            method="from_meta_file", meta_file=quakeml, st=st, lowcut=2.0,
            highcut=9.0, samp_rate=20.0, filt_order=3, length=2, prepick=0.1,
            swin='S')
        self.assertEqual(len(templates), 1)

    def test_from_client(self):
        """Test for using a waveform client not related to obspy's Clients."""
        client = _StreamTestClient()
        cat = client.get_default_catalog()
        temps = template_gen('from_client', client_id=client, catalog=cat,
                             highcut=None, lowcut=None, filt_order=4,
                             samp_rate=100, prepick=0.1, length=10,
                             process_len=20, data_pad=5)
        self.assertEqual(len(temps), 1)  # there should be one template stream
        self.assertEqual(len(temps[0]), 1)  # with one trace
        self.assertGreater(len(temps[0][0].data), 1)  # with some data

    def test_bad_client(self):
        """Ensure passing a non-client raises."""
        client = _StreamTestClient()
        cat = client.get_default_catalog()
        with self.assertRaises(NotImplementedError):
            template_gen('from_client', client_id=cat, catalog=cat,
                         highcut=None, lowcut=None, filt_order=4,
                         samp_rate=100, prepick=0.1, length=10,
                         process_len=20, data_pad=5)


class TestEdgeGen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import eqcorrscan
        cls.testing_path = os.path.dirname(eqcorrscan.__file__) + '/tests'
        log = logging.getLogger(template_gen_module.__name__)
        cls._log_handler = MockLoggingHandler(level='DEBUG')
        log.addHandler(cls._log_handler)
        cls.log_messages = cls._log_handler.messages
        cls.st = read(os.path.join(
            cls.testing_path, 'test_data', 'WAV', 'TEST_',
            '2013-09-15-0930-28.DFDPC_027_00'))
        for tr in cls.st:
            tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        event = read_events(os.path.join(
            cls.testing_path, 'test_data', 'REA', 'TEST_',
            '15-0931-08L.S201309'))[0]
        cls.picks = event.picks

    def setUp(self):
        self._log_handler.reset()

    def test_undefined_phase_type(self):
        with self.assertRaises(AssertionError):
            _template_gen(
                picks=self.picks, st=self.st.copy(), length=2, swin='bob')

    def test_warn_zeros(self):
        st = self.st.copy()
        template = _template_gen(self.picks, st.copy(), 10)
        self.assertTrue('LABE' in [tr.stats.station for tr in template])
        st.select(station='LABE', channel='SN')[0].data = np.zeros(10000)
        template = _template_gen(self.picks, st, 10)
        self.assertFalse('LABE' in [tr.stats.station for tr in template])

    def test_missing_data(self):
        picks = copy.deepcopy(self.picks)
        picks.append(picks[-1])
        picks[-1].waveform_id.station_code = 'DUMMY'
        template = _template_gen(picks, self.st.copy(), 10)
        self.assertFalse('DUMMY' in [tr.stats.station for tr in template])

    def test_no_matched_picks(self):
        picks = [copy.deepcopy(self.picks[0])]
        picks[0].waveform_id.station_code = 'DUMMY'
        template = _template_gen(picks, self.st.copy(), 10)
        self.assertFalse(template)

    def test_missing_phase_hints(self):
        picks = copy.deepcopy(self.picks)
        template = _template_gen(picks, self.st.copy(), 10, swin="P_all")
        self.assertIn("WV03", {tr.stats.station for tr in template})
        # Remove phase_hint for WV03 P pick
        for pick in picks:
            if pick.waveform_id.station_code == "WV03":
                pick.phase_hint = None
        template_2 = _template_gen(picks, self.st.copy(), 10, swin="P_all")
        self.assertNotIn("WV03", {tr.stats.station for tr in template_2})
        for tr in template:
            if tr.stats.station == "WV03":
                continue
            # check that other than WV03 the template is the same
            tr2 = template_2.select(id=tr.id)[0]
            self.assertEqual(tr, tr2)

    def test_misc(self):
        template = _template_gen(self.picks, self.st.copy(), 10)
        self.assertEqual(len(template), len(self.picks))

    def test_extract_from_stack(self):
        length = 3
        stack = self.st.copy()
        template = _template_gen(self.picks, self.st.copy(), 2)
        extracted = extract_from_stack(stack, template, length=length,
                                       pre_pick=0.3, pre_pad=45)
        self.assertEqual(len(template), len(extracted))
        for tr in extracted:
            self.assertEqual(tr.stats.endtime - tr.stats.starttime,
                             length)

    def test_extract_from_stack_and_process(self):
        length = 3
        stack = self.st.copy()
        template = _template_gen(self.picks, self.st.copy(), 2)
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
        template = _template_gen(self.picks, self.st.copy(), 2)
        extracted = extract_from_stack(
            stack, template, length=length, pre_pick=0.3, pre_pad=45,
            Z_include=True)
        self.assertLess(len(template), len(extracted))
        for tr in extracted:
            self.assertEqual(tr.stats.endtime - tr.stats.starttime,
                             length)

    def test_warn_no_phase_hint(self):
        picks = copy.deepcopy(self.picks)
        for pick in picks:
            setattr(pick, 'phase_hint', None)
        template = _template_gen(picks, self.st.copy(), 10)
        w = self.log_messages['warning']
        self.assertGreater(len(w), 0)
        self.assertEqual(len(template), 11)
        _w = ' '.join([warning for warning in w])
        self.assertTrue("no phase hint given" in _w)

    def test_no_station_code(self):
        picks = copy.deepcopy(self.picks)
        for pick in picks:
            setattr(pick.waveform_id, 'station_code', None)
        template = _template_gen(picks, self.st.copy(), 10)
        self.assertEqual(len(template), 0)

    def test_no_channel_code(self):
        picks = copy.deepcopy(self.picks)
        for pick in picks:
            setattr(pick.waveform_id, 'channel_code', None)
        template = _template_gen(picks, self.st.copy(), 10)
        self.assertEqual(len(template), 0)

    def test_swin_P(self):
        p_picks = []
        for pick in self.picks:
            if pick.phase_hint == 'P':
                p_picks.append(pick)
        template = _template_gen(self.picks, self.st.copy(), 10, swin='P')
        self.assertEqual(len(template), len(p_picks))
        used_stations = [tr.stats.station for tr in template]
        for pick in p_picks:
            self.assertTrue(pick.waveform_id.station_code in used_stations)

    def test_swin_all_and_all_horiz(self):
        template = _template_gen(self.picks, self.st.copy(), 10, swin='all',
                                 all_horiz=True,
                                 horizontal_chans=['E', 'N', '1', '2', '3'])
        for pick in self.picks:
            if pick.phase_hint == 'S':
                self.assertGreaterEqual(
                    len(template.select(
                        station=pick.waveform_id.station_code)), 2)

    def test_snr_cutoff(self):
        template = _template_gen(self.picks, self.st.copy(), 10, min_snr=100)
        self.assertEqual(len(template), 1)

    def test_no_data_for_channel(self):
        st = self.st.copy()
        st.select(station='LABE', channel='SN')[0].stats.starttime += 2000
        template = _template_gen(self.picks, st, 10)
        self.assertEqual(len(template), 10)

    @pytest.mark.network
    def test_triggered_data(self):
        client = Client("GEONET")
        catalog = client.get_events(eventid="1481730")
        templates = template_gen(
            "from_client", lowcut=2., highcut=15., samp_rate=40., swin="all",
            filt_order=4, prepick=0.2, catalog=catalog, length=3.0,
            client_id="GEONET", all_horiz=True, process_len=600,
            min_snr=5., skip_short_chans=True)
        self.assertEqual(len(templates), 0)


class TestEdgeGenObs(unittest.TestCase):
    @classmethod
    # Extra test case with OBS data with hydrophone channels (HDH) and T-phases
    def setUpClass(cls):
        import eqcorrscan
        cls.testing_path = os.path.dirname(eqcorrscan.__file__) + '/tests'
        log = logging.getLogger(template_gen_module.__name__)
        cls._log_handler = MockLoggingHandler(level='DEBUG')
        log.addHandler(cls._log_handler)
        cls.log_messages = cls._log_handler.messages
        cls.st = read(os.path.join(
            cls.testing_path, 'test_data', 'WAV', 'TEST_',
            '2019-08-09-1558-47M.NNSN__038'))
        # for tr in cls.st:
        #     tr.stats.channel = tr.stats.channel[0] + tr.stats.channel[-1]
        # Sfile in New Nordic format
        event = read_events(os.path.join(
            cls.testing_path, 'test_data', 'REA', 'TEST_',
            '09-1558-48R.S201908'))[0]
        cat = filter_picks(
            Catalog([event]), stations=['KBS', 'OBIN1', 'OBIN2', 'SPA0',
                                        'NOR', 'DAG', 'HOPEN', 'HSPB'])
        cls.picks = cat[0].picks

    def setUp(self):
        self._log_handler.reset()

    def test_swin_all_and_all_vert_and_all_horiz(self):
        # Test that the hydrophone channel on an OBS is included in the
        # creation of the vertical channel (P-arrival) template.
        template = _template_gen(self.picks, self.st.copy(), 20, swin='all',
                                 all_horiz=True, all_vert=True,
                                 vertical_chans=['Z', 'H'])
        for pick in self.picks:
            if pick.phase_hint and pick.phase_hint[0] == 'P':
                self.assertGreaterEqual(
                    len(template.select(
                        station=pick.waveform_id.station_code,
                        channel='??[ZH]')), 1)
            if pick.phase_hint and pick.phase_hint[0] == 'S':
                self.assertGreaterEqual(
                    len(template.select(
                        station=pick.waveform_id.station_code,
                        channel='??[NE12]')), 2)
            if pick.phase_hint and pick.phase_hint[0] == 'T':
                self.assertGreaterEqual(
                    len(template.select(
                        station=pick.waveform_id.station_code,
                        channel='??[ZH]')), 2)


class TestDayLong(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testing_path = os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe())))
        cls.st = read(os.path.join(cls.testing_path, 'test_data',
                                   'day_vols', 'Y2012', 'R086.01', '*'))
        event = read_events(os.path.join(
            cls.testing_path, 'test_data', 'REA', 'TEST_',
            '15-0931-08L.S201309'))[0]
        event.picks = event.picks[0:3]
        for pick, station in zip(event.picks, ['GOVA', 'EORO', 'WHYM']):
            setattr(pick, 'time', UTCDateTime(2012, 3, 26, 2, 0, 0))
            setattr(pick.waveform_id, 'channel_code', 'SHZ')
            setattr(pick.waveform_id, 'network_code', 'AF')
            setattr(pick.waveform_id, 'station_code', station)
            setattr(pick, 'phase_hint', 'P')
        cls.cat = Catalog([event])

    @pytest.mark.slow
    def test_day_long_processing(self):
        templates = template_gen(
            method='from_meta_file', meta_file=self.cat, st=self.st,
            lowcut=2.0, highcut=9.0, samp_rate=20.0, filt_order=3, length=2,
            prepick=0.1, swin='P', all_horiz=True)
        self.assertEqual(len(templates), 1)
        self.assertEqual(len(templates[0]), 3)


if __name__ == '__main__':
    unittest.main()
