"""
Functions to test the mag_calc functions within EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import pytest
import numpy as np
import os
import datetime as dt
import glob
import logging

from obspy.core.util import NamedTemporaryFile
from obspy import UTCDateTime, read, read_events
from obspy.clients.fdsn import Client
from obspy.clients.iris import Client as OldIris_Client
from obspy.io.nordic.core import readwavename

from eqcorrscan.utils import mag_calc
from eqcorrscan.utils.mag_calc import dist_calc, _sim_WA, _max_p2t
from eqcorrscan.utils.mag_calc import _GSE2_PAZ_read, _find_resp, _pairwise
from eqcorrscan.utils.mag_calc import svd_moments, amp_pick_event
from eqcorrscan.utils.clustering import svd
from eqcorrscan.helpers.mock_logger import MockLoggingHandler


class TestMagCalcMethods(unittest.TestCase):
    """Test all mag_calc functions."""
    def test_dist_calc(self):
        """
        Test the distance calculation that comes with mag_calc.
        """
        self.assertEqual(dist_calc((0, 0, 0), (0, 0, 10)), 10)
        self.assertEqual(round(dist_calc((0, 0, 0), (0, 1, 0))), 111)
        self.assertEqual(round(dist_calc((0, 0, 0), (1, 0, 0))), 111)
        self.assertEqual(round(dist_calc((45, 45, 0), (45, 45, 10))), 10)
        self.assertEqual(round(dist_calc((45, 45, 0), (45, 46, 0))), 79)
        self.assertEqual(round(dist_calc((45, 45, 0), (46, 45, 0))), 111)
        self.assertEqual(round(dist_calc((90, 90, 0), (90, 90, 10))), 10)
        self.assertEqual(round(dist_calc((90, 90, 0), (90, 89, 0))), 0)
        self.assertEqual(round(dist_calc((90, 90, 0), (89, 90, 0))), 111)
        self.assertEqual(round(dist_calc((0, 180, 0), (0, -179, 0))), 111)
        self.assertEqual(round(dist_calc((0, -180, 0), (0, 179, 0))), 111)

    @pytest.mark.network
    @pytest.mark.flaky(reruns=2)
    def test_sim_WA(self):
        """Test feeding both PAZ and seedresp."""
        t1 = UTCDateTime("2010-09-3T16:30:00.000")
        t2 = UTCDateTime("2010-09-3T17:00:00.000")
        fdsn_client = Client('IRIS')
        st = fdsn_client.get_waveforms(
            network='NZ', station='BFZ', location='10', channel='HHZ',
            starttime=t1, endtime=t2, attach_response=True)
        tr = st[0]
        PAZ = {'poles': [-4.440 + 4.440j, -4.440 - 4.440j, -1.083 + 0.0j],
               'zeros': [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0],
               'sensitivity': 0.4,
               'gain': 60077000.0}
        tr_safe = tr.copy()
        # Test with PAZ
        _sim_WA(trace=tr, PAZ=PAZ, seedresp=None, water_level=10)
        tr = tr_safe.copy()
        # Test without PAZ or seedresp
        _sim_WA(trace=tr, PAZ=None, seedresp=None, water_level=10)
        tr = tr_safe.copy()
        with open("Temp_resp", "w") as tf:
            respf = tf.name
            old_iris_client = OldIris_Client()
            # fetch RESP information from "old" IRIS web service, see
            # obspy.fdsn for accessing the new IRIS FDSN web services
            old_iris_client.resp('NZ', 'BFZ', '10', 'HHZ', t1, t2,
                                 filename=respf)
            # Hack around unit issues
        with open("Temp_resp", "r") as tf:
            resp_contents = [line for line in tf]
        corrected_contents = []
        for line in resp_contents:
            if "COUNT" in line:
                line = line.replace("COUNT", "COUNTS")
            corrected_contents.append(line)
        with open("Temp_resp", "w") as tf:
            for line in corrected_contents:
                tf.write(line)
        date = t1

        seedresp = {
            'filename': respf, 'date': date, 'network': tr.stats.network,
            'station': tr.stats.station, 'channel': tr.stats.channel,
            'location': tr.stats.location, 'units': 'DIS'}
        _sim_WA(trace=tr, PAZ=None, seedresp=seedresp, water_level=10)
        os.remove(respf)

    def test_max_p2t(self):
        """Test the minding of maximum peak-to-trough."""
        data = np.random.randn(1000)
        data[500] = 20
        delta = 0.01
        amplitude, period, time = _max_p2t(data=data, delta=delta)
        self.assertTrue(amplitude > 15)
        self.assertTrue(amplitude < 25)
        self.assertEqual(round(time), 5)

    def test_GSE_read(self):
        """Test reading GSE PAZ."""
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        GSEfile = os.path.join(testing_path, 'POCR2SH_1.2008-01-01-0000_GSE')
        PAZ, date, station, channel, sensor = _GSE2_PAZ_read(gsefile=GSEfile)
        # Check that all the elements are there
        self.assertEqual(date, dt.datetime(2008, 11, 6, 0, 0))
        self.assertEqual(station, 'POCR2')
        self.assertEqual(channel, 'SH1')
        self.assertEqual(sensor, 'HS1-2')
        self.assertTrue('gain' in PAZ)
        self.assertTrue('poles' in PAZ)
        self.assertTrue('sensitivity' in PAZ)
        self.assertTrue('zeros' in PAZ)
        # Check that we only cope with CAL2 files
        with NamedTemporaryFile() as tf:
            f = open(GSEfile, 'r')
            corrf = open(tf.name, 'w')
            corrf.write(f.readline().replace('CAL2', 'BOB2'))
            for line in f:
                corrf.write(line)
            corrf.close()
            f.close()
            with self.assertRaises(IOError):
                _GSE2_PAZ_read(gsefile=tf.name)
        # Check that we only cope with PAZ2 files
        with NamedTemporaryFile() as tf:
            f = open(GSEfile, 'r')
            corrf = open(tf.name, 'w')
            corrf.write(f.readline())
            corrf.write(f.readline().replace('PAZ2', 'BOB2'))
            for line in f:
                corrf.write(line)
            corrf.close()
            f.close()
            with self.assertRaises(IOError):
                _GSE2_PAZ_read(gsefile=tf.name)

    def test_find_resp(self):
        """Test the ability to find response info"""
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        station = 'POCR2'
        channel = 'SH1'
        network = 'AF'
        time = dt.datetime(2008, 11, 9, 0, 0)

        PAZ = _find_resp(station=station, channel=channel, network=network,
                         time=time, delta=0.005, directory=testing_path)
        # Will find the POCR2 respfile
        self.assertTrue('gain' in PAZ)
        self.assertTrue('poles' in PAZ)
        self.assertTrue('sensitivity' in PAZ)
        self.assertTrue('zeros' in PAZ)

        station = 'GCSZ'
        channel = 'EHZ'
        network = 'NZ'
        time = dt.datetime(2013, 1, 1, 0, 0)

        resp = _find_resp(station=station, channel=channel, network=network,
                          time=time, delta=0.005, directory=testing_path)
        for key in ['channel', 'date', 'filename', 'location',
                    'network', 'station', 'units']:
            self.assertTrue(key in resp)

        station = 'WZ01'
        channel = 'ELE'
        network = 'ZT'
        time = dt.datetime(2013, 1, 1, 0, 0)

        resp = _find_resp(station=station, channel=channel, network=network,
                          time=time, delta=0.005, directory=testing_path)
        for key in ['channel', 'date', 'filename', 'location',
                    'network', 'station', 'units']:
            self.assertTrue(key in resp)

    def test_pairwise(self):
        """Test the itertools wrapper"""
        pairs = _pairwise(range(20))
        for i, pair in enumerate(pairs):
            self.assertEqual(pair[0] + 1, pair[1])
        self.assertEqual(i, 18)

    def test_SVD_mag(self):
        """Test the SVD magnitude calculator."""
        # Do the set-up
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events_processed')
        stream_files = glob.glob(os.path.join(testing_path, '*DFDPC*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        event_list = []
        for i, stream in enumerate(stream_list):
            st_list = []
            for tr in stream:
                if (tr.stats.station, tr.stats.channel) not in\
                        [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]:
                    stream.remove(tr)
                    continue
                st_list.append(i)
            event_list.append(st_list)
        event_list = np.asarray(event_list).T.tolist()
        SVectors, SValues, Uvectors, stachans = svd(stream_list=stream_list)
        M, events_out = svd_moments(u=Uvectors, s=SValues, v=SVectors,
                                    stachans=stachans, event_list=event_list)
        self.assertEqual(len(M), len(stream_list))
        self.assertEqual(len(events_out), len(stream_list))


class TestAmpPickEvent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        log = logging.getLogger(mag_calc.__name__)
        cls._log_handler = MockLoggingHandler(level='DEBUG')
        log.addHandler(cls._log_handler)
        cls.log_messages = cls._log_handler.messages
        cls.testing_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'test_data')
        sfile = os.path.join(cls.testing_path, 'REA', 'TEST_',
                             '01-0411-15L.S201309')
        cls.event = read_events(sfile)[0]
        cls.wavfiles = readwavename(sfile)
        cls.datapath = os.path.join(cls.testing_path, 'WAV', 'TEST_')
        cls.st = read(os.path.join(cls.datapath, cls.wavfiles[0]))
        cls.respdir = cls.testing_path

    def setUp(self):
        self._log_handler.reset()

    def test_amp_pick_event(self):
        """Test the main amplitude picker."""
        picked_event = amp_pick_event(event=self.event.copy(),
                                      st=self.st.copy(),
                                      respdir=self.respdir)
        self.assertEqual(len(picked_event.picks), len(self.event.picks) + 1)

    def test_amp_pick_remove_old_picks(self):
        picked_event = amp_pick_event(event=self.event.copy(),
                                      st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_missing_channel(self):
        picked_event = amp_pick_event(
            event=self.event.copy(), st=self.st.copy()[0:-4],
            respdir=self.respdir, remove_old=True)
        missed = False
        for warning in self.log_messages['warning']:
            if 'no station and channel match' in warning:
                missed = True
        self.assertTrue(missed)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_not_varwin(self):
        picked_event = amp_pick_event(event=self.event.copy(),
                                      st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True,
                                      var_wintype=False)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_not_varwin_no_S(self):
        event = self.event.copy()
        for pick in event.picks:
            if pick.phase_hint.upper() == 'S':
                event.picks.remove(pick)
        picked_event = amp_pick_event(event=event, st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True,
                                      var_wintype=False)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_varwin_no_S(self):
        event = self.event.copy()
        for pick in event.picks:
            if pick.phase_hint.upper() == 'S':
                event.picks.remove(pick)
        picked_event = amp_pick_event(event=event, st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True,
                                      var_wintype=True)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_varwin_no_P(self):
        event = self.event.copy()
        for pick in event.picks:
            if pick.phase_hint.upper() == 'P':
                event.picks.remove(pick)
        picked_event = amp_pick_event(event=event, st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True,
                                      var_wintype=True)
        self.assertEqual(len(picked_event.amplitudes), 1)

    def test_amp_pick_high_min_snr(self):
        picked_event = amp_pick_event(event=self.event.copy(),
                                      st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True,
                                      var_wintype=False, min_snr=15)
        self.assertEqual(len(picked_event.amplitudes), 0)

    def test_amp_pick_no_prefilt(self):
        picked_event = amp_pick_event(event=self.event.copy(),
                                      st=self.st.copy(),
                                      respdir=self.respdir, remove_old=True,
                                      var_wintype=False, pre_filt=False)
        self.assertEqual(len(picked_event.amplitudes), 1)


if __name__ == '__main__':
    unittest.main()
