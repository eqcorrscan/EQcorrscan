"""
Functions to test the mag_calc functions within EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestMagCalcMethods(unittest.TestCase):
    """Test all mag_calc functions."""
    def test_dist_calc(self):
        """
        Test the distance calculation that comes with mag_calc.
        """
        from eqcorrscan.utils.mag_calc import dist_calc
        self.assertEqual(dist_calc((0, 0, 0), (0, 0, 10)), 10)
        self.assertEqual(round(dist_calc((0, 0, 0), (0, 1, 0))), 111)
        self.assertEqual(round(dist_calc((0, 0, 0), (1, 0, 0))), 111)
        self.assertEqual(round(dist_calc((45, 45, 0), (45, 45, 10))), 10)
        self.assertEqual(round(dist_calc((45, 45, 0), (45, 46, 0))), 79)
        self.assertEqual(round(dist_calc((45, 45, 0), (46, 45, 0))), 111)
        self.assertEqual(round(dist_calc((90, 90, 0), (90, 90, 10))), 10)
        self.assertEqual(round(dist_calc((90, 90, 0), (90, 89, 0))), 0)
        self.assertEqual(round(dist_calc((90, 90, 0), (89, 90, 0))), 111)

    def test_sim_WA(self):
        """Test feeding both PAZ and seedresp."""
        from eqcorrscan.utils.mag_calc import _sim_WA
        from obspy.core.util import NamedTemporaryFile
        from obspy import UTCDateTime
        from obspy.clients.fdsn import Client
        from obspy.clients.iris import Client as OldIris_Client

        t1 = UTCDateTime("2010-09-3T16:30:00.000")
        t2 = UTCDateTime("2010-09-3T17:00:00.000")
        fdsn_client = Client('IRIS')
        st = fdsn_client.get_waveforms(network='NZ', station='BFZ',
                                       location='10',
                                       channel='HHZ',
                                       starttime=t1, endtime=t2,
                                       attach_response=True)
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
        with NamedTemporaryFile() as tf:
            respf = tf.name
            old_iris_client = OldIris_Client()
            # fetch RESP information from "old" IRIS web service, see
            # obspy.fdsn for accessing the new IRIS FDSN web services
            old_iris_client.resp('NZ', 'BFZ', '10', 'HHZ', t1, t2,
                                 filename=respf)
            date = t1

            seedresp = {'filename': respf,  # RESP filename
                        'date': date,
                        'network': tr.stats.network,
                        'station': tr.stats.station,
                        'channel': tr.stats.channel,
                        'location': tr.stats.location,
                        # Units to return response in ('DIS', 'VEL' or ACC)
                        'units': 'DIS'
                        }
            _sim_WA(trace=tr, PAZ=None, seedresp=seedresp, water_level=10)

    def test_max_p2t(self):
        """Test the minding of maximum peak-to-trough."""
        import numpy as np
        from eqcorrscan.utils.mag_calc import _max_p2t
        data = np.random.randn(1000)
        data[500] = 20
        delta = 0.01
        amplitude, period, time = _max_p2t(data=data, delta=delta)
        self.assertTrue(amplitude > 15)
        self.assertTrue(amplitude < 25)
        self.assertEqual(round(time), 5)

    def test_GSE_read(self):
        """Test reading GSE PAZ."""
        import os
        from eqcorrscan.utils.mag_calc import _GSE2_PAZ_read
        import datetime as dt
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

    def test_find_resp(self):
        """Test the ability to find response info"""
        from eqcorrscan.utils.mag_calc import _find_resp
        import datetime as dt
        import os
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
        from eqcorrscan.utils.mag_calc import _pairwise
        pairs = _pairwise(range(20))
        for i, pair in enumerate(pairs):
            self.assertEqual(pair[0] + 1, pair[1])
        self.assertEqual(i, 18)

    def test_amp_pick_sfile(self):
        """Test the amplitude picker wrapper."""
        from eqcorrscan.utils.mag_calc import amp_pick_sfile
        from obspy.core.event import Event
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        sfile = os.path.join(testing_path, 'REA', 'TEST_',
                             '01-0411-15L.S201309')
        datapath = os.path.join(testing_path, 'WAV', 'TEST_')
        respdir = testing_path
        event = amp_pick_sfile(sfile=sfile, datapath=datapath,
                               respdir=respdir, chans=['Z'], var_wintype=True,
                               winlen=0.9, pre_pick=0.2, pre_filt=True,
                               lowcut=1.0, highcut=20.0, corners=4)
        self.assertTrue(isinstance(event, Event))
        self.assertTrue(os.path.isfile('mag_calc.out'))
        os.remove('mag_calc.out')

    def test_SVD_mag(self):
        """Test the SVD magnitude calcualtor."""
        from eqcorrscan.utils.mag_calc import SVD_moments
        from obspy import read
        import glob
        import os
        from eqcorrscan.utils.clustering import SVD
        import numpy as np
        # Do the set-up
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'similar_events')
        stream_files = glob.glob(os.path.join(testing_path, '*'))
        stream_list = [read(stream_file) for stream_file in stream_files]
        event_list = []
        for i, stream in enumerate(stream_list):
            st_list = []
            for tr in stream:
                if (tr.stats.station, tr.stats.channel) not in\
                        [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]:
                    stream.remove(tr)
                    continue
                tr.detrend('simple')
                tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
                tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
                st_list.append(i)
            event_list.append(st_list)
        event_list = np.asarray(event_list).T.tolist()
        SVectors, SValues, Uvectors, stachans = SVD(stream_list=stream_list)
        M, events_out = SVD_moments(U=Uvectors, s=SValues, V=SVectors,
                                    stachans=stachans, event_list=event_list)
        self.assertEqual(len(M), len(stream_list))
        self.assertEqual(len(events_out), len(stream_list))

if __name__ == '__main__':
    unittest.main()
