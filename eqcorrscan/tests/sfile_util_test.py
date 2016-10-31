"""
Functions for testing the utils.sfile_util functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.utils.sfile_util import eventtosfile, readwavename, readpicks
from eqcorrscan.utils.sfile_util import _nortoevmag, _evmagtonor, nordpick
from eqcorrscan.utils.sfile_util import _int_conv, _float_conv, _str_conv
from eqcorrscan.utils.sfile_util import read_event, read_select, blanksfile
import unittest


class TestSfileMethods(unittest.TestCase):
    def test_download_write(self):
        """
        Function to download quakeML files from a range of datacenters and \
        attempt to write miniseed files
        """
        import os
        from eqcorrscan.utils import sfile_util
        import obspy
        if int(obspy.__version__.split('.')[0]) >= 1:
            from obspy.clients.fdsn import Client
            from obspy import read_events
            from obspy.clients.fdsn.header import FDSNException
        else:
            from obspy.fdsn import Client
            from obspy import readEvents as read_events
            from obspy.fdsn.header import FDSNException
        import warnings

        event_list = [('GEONET', '2016p008122'),
                      ('NCEDC', '72572665'),
                      ('USGS', 'nc72597260')]
        for event_info in event_list:
            try:
                client = Client(event_info[0])
                if event_info[0] == 'GEONET':
                    data_stream = client.\
                        _download('http://quakeml.geonet.org.nz/' +
                                  'quakeml/1.2/' + event_info[1])
                    data_stream.seek(0, 0)
                    event = read_events(data_stream, format="quakeml")
                    data_stream.close()
                else:
                    event = client.get_events(eventid=event_info[1],
                                              includearrivals=True)
            except FDSNException:
                warnings.warn('FDSNException')
                continue
            test_Sfile_name = sfile_util.eventtosfile(event, 'test', 'L', '.',
                                                      'null', overwrite=True)
            os.remove(test_Sfile_name)
        return True

    def test_read_write(self):
        """
        Function to test the read and write capabilities of sfile_util.
        """
        import os
        from obspy.core.event import Catalog
        import obspy
        if int(obspy.__version__.split('.')[0]) >= 1:
            from obspy.core.event import read_events
        else:
            from obspy.core.event import readEvents as read_events

        # Set-up a test event
        test_event = full_test_event()
        # Add the event to a catalogue which can be used for QuakeML testing
        test_cat = Catalog()
        test_cat += test_event
        # Write the catalog
        test_cat.write("Test_catalog.xml", format='QUAKEML')
        # Read and check
        read_cat = read_events("Test_catalog.xml")
        os.remove("Test_catalog.xml")
        self.assertEqual(read_cat[0].resource_id, test_cat[0].resource_id)
        for i in range(len(read_cat[0].picks)):
            for key in read_cat[0].picks[i].keys():
                # Ignore backazimuth errors and horizontal_slowness_errors
                if key in ['backazimuth_errors', 'horizontal_slowness_errors']:
                    continue
                self.assertEqual(read_cat[0].picks[i][key],
                                 test_cat[0].picks[i][key])
        self.assertEqual(read_cat[0].origins[0].resource_id,
                         test_cat[0].origins[0].resource_id)
        self.assertEqual(read_cat[0].origins[0].time,
                         test_cat[0].origins[0].time)
        # Note that time_residual_RMS is not a quakeML format
        self.assertEqual(read_cat[0].origins[0].longitude,
                         test_cat[0].origins[0].longitude)
        self.assertEqual(read_cat[0].origins[0].latitude,
                         test_cat[0].origins[0].latitude)
        self.assertEqual(read_cat[0].origins[0].depth,
                         test_cat[0].origins[0].depth)
        # Check magnitudes
        self.assertEqual(read_cat[0].magnitudes, test_cat[0].magnitudes)
        self.assertEqual(read_cat[0].event_descriptions,
                         test_cat[0].event_descriptions)
        # Check local magnitude amplitude
        self.assertEqual(read_cat[0].amplitudes[0].resource_id,
                         test_cat[0].amplitudes[0].resource_id)
        self.assertEqual(read_cat[0].amplitudes[0].period,
                         test_cat[0].amplitudes[0].period)
        self.assertEqual(read_cat[0].amplitudes[0].unit,
                         test_cat[0].amplitudes[0].unit)
        self.assertEqual(read_cat[0].amplitudes[0].generic_amplitude,
                         test_cat[0].amplitudes[0].generic_amplitude)
        self.assertEqual(read_cat[0].amplitudes[0].pick_id,
                         test_cat[0].amplitudes[0].pick_id)
        self.assertEqual(read_cat[0].amplitudes[0].waveform_id,
                         test_cat[0].amplitudes[0].waveform_id)
        # Check coda magnitude pick
        self.assertEqual(read_cat[0].amplitudes[1].resource_id,
                         test_cat[0].amplitudes[1].resource_id)
        self.assertEqual(read_cat[0].amplitudes[1].type,
                         test_cat[0].amplitudes[1].type)
        self.assertEqual(read_cat[0].amplitudes[1].unit,
                         test_cat[0].amplitudes[1].unit)
        self.assertEqual(read_cat[0].amplitudes[1].generic_amplitude,
                         test_cat[0].amplitudes[1].generic_amplitude)
        self.assertEqual(read_cat[0].amplitudes[1].pick_id,
                         test_cat[0].amplitudes[1].pick_id)
        self.assertEqual(read_cat[0].amplitudes[1].waveform_id,
                         test_cat[0].amplitudes[1].waveform_id)
        self.assertEqual(read_cat[0].amplitudes[1].magnitude_hint,
                         test_cat[0].amplitudes[1].magnitude_hint)
        self.assertEqual(read_cat[0].amplitudes[1].snr,
                         test_cat[0].amplitudes[1].snr)
        self.assertEqual(read_cat[0].amplitudes[1].category,
                         test_cat[0].amplitudes[1].category)

        # Check the read-write s-file functionality
        sfile = eventtosfile(test_cat[0], userID='TEST',
                             evtype='L', outdir='.',
                             wavefiles='test', explosion=True, overwrite=True)
        del read_cat
        self.assertEqual(readwavename(sfile), ['test'])
        read_cat = Catalog()
        read_cat += readpicks(sfile)
        os.remove(sfile)
        for i in range(len(read_cat[0].picks)):
            self.assertEqual(read_cat[0].picks[i].time,
                             test_cat[0].picks[i].time)
            self.assertEqual(read_cat[0].picks[i].backazimuth,
                             test_cat[0].picks[i].backazimuth)
            self.assertEqual(read_cat[0].picks[i].onset,
                             test_cat[0].picks[i].onset)
            self.assertEqual(read_cat[0].picks[i].phase_hint,
                             test_cat[0].picks[i].phase_hint)
            self.assertEqual(read_cat[0].picks[i].polarity,
                             test_cat[0].picks[i].polarity)
            self.assertEqual(read_cat[0].picks[i].waveform_id.station_code,
                             test_cat[0].picks[i].waveform_id.station_code)
            self.assertEqual(read_cat[0].picks[i].waveform_id.channel_code[-1],
                             test_cat[0].picks[i].waveform_id.channel_code[-1])
        # assert read_cat[0].origins[0].resource_id ==\
        #     test_cat[0].origins[0].resource_id
        self.assertEqual(read_cat[0].origins[0].time,
                         test_cat[0].origins[0].time)
        # Note that time_residual_RMS is not a quakeML format
        self.assertEqual(read_cat[0].origins[0].longitude,
                         test_cat[0].origins[0].longitude)
        self.assertEqual(read_cat[0].origins[0].latitude,
                         test_cat[0].origins[0].latitude)
        self.assertEqual(read_cat[0].origins[0].depth,
                         test_cat[0].origins[0].depth)
        self.assertEqual(read_cat[0].magnitudes[0].mag,
                         test_cat[0].magnitudes[0].mag)
        self.assertEqual(read_cat[0].magnitudes[1].mag,
                         test_cat[0].magnitudes[1].mag)
        self.assertEqual(read_cat[0].magnitudes[2].mag,
                         test_cat[0].magnitudes[2].mag)
        self.assertEqual(read_cat[0].magnitudes[0].creation_info,
                         test_cat[0].magnitudes[0].creation_info)
        self.assertEqual(read_cat[0].magnitudes[1].creation_info,
                         test_cat[0].magnitudes[1].creation_info)
        self.assertEqual(read_cat[0].magnitudes[2].creation_info,
                         test_cat[0].magnitudes[2].creation_info)
        self.assertEqual(read_cat[0].magnitudes[0].magnitude_type,
                         test_cat[0].magnitudes[0].magnitude_type)
        self.assertEqual(read_cat[0].magnitudes[1].magnitude_type,
                         test_cat[0].magnitudes[1].magnitude_type)
        self.assertEqual(read_cat[0].magnitudes[2].magnitude_type,
                         test_cat[0].magnitudes[2].magnitude_type)
        self.assertEqual(read_cat[0].event_descriptions,
                         test_cat[0].event_descriptions)
        # assert read_cat[0].amplitudes[0].resource_id ==\
        #     test_cat[0].amplitudes[0].resource_id
        self.assertEqual(read_cat[0].amplitudes[0].period,
                         test_cat[0].amplitudes[0].period)
        self.assertEqual(read_cat[0].amplitudes[0].snr,
                         test_cat[0].amplitudes[0].snr)
        # Check coda magnitude pick
        # Resource ids get overwritten because you can't have two the same in
        # memory
        # self.assertEqual(read_cat[0].amplitudes[1].resource_id,
        #                  test_cat[0].amplitudes[1].resource_id)
        self.assertEqual(read_cat[0].amplitudes[1].type,
                         test_cat[0].amplitudes[1].type)
        self.assertEqual(read_cat[0].amplitudes[1].unit,
                         test_cat[0].amplitudes[1].unit)
        self.assertEqual(read_cat[0].amplitudes[1].generic_amplitude,
                         test_cat[0].amplitudes[1].generic_amplitude)
        # Resource ids get overwritten because you can't have two the same in
        # memory
        # self.assertEqual(read_cat[0].amplitudes[1].pick_id,
        #                  test_cat[0].amplitudes[1].pick_id)
        self.assertEqual(read_cat[0].amplitudes[1].waveform_id.station_code,
                         test_cat[0].amplitudes[1].waveform_id.station_code)
        self.assertEqual(read_cat[0].amplitudes[1].waveform_id.channel_code,
                         test_cat[0].amplitudes[1].
                         waveform_id.channel_code[0] +
                         test_cat[0].amplitudes[1].
                         waveform_id.channel_code[-1])
        self.assertEqual(read_cat[0].amplitudes[1].magnitude_hint,
                         test_cat[0].amplitudes[1].magnitude_hint)
        # snr is not supported in s-file
        # self.assertEqual(read_cat[0].amplitudes[1].snr,
        #                  test_cat[0].amplitudes[1].snr)
        self.assertEqual(read_cat[0].amplitudes[1].category,
                         test_cat[0].amplitudes[1].category)
        del read_cat

        # Test a deliberate fail
        test_cat.append(full_test_event())
        with self.assertRaises(IOError):
            # Raises error due to multiple events in catalog
            sfile = eventtosfile(test_cat, userID='TEST',
                                 evtype='L', outdir='.',
                                 wavefiles='test', explosion=True,
                                 overwrite=True)
            # Raises error due to too long userID
            sfile = eventtosfile(test_cat[0], userID='TESTICLE',
                                 evtype='L', outdir='.',
                                 wavefiles='test', explosion=True,
                                 overwrite=True)
            # Raises error due to unrecognised event type
            sfile = eventtosfile(test_cat[0], userID='TEST',
                                 evtype='U', outdir='.',
                                 wavefiles='test', explosion=True,
                                 overwrite=True)
            # Raises error due to no output directory
            sfile = eventtosfile(test_cat[0], userID='TEST',
                                 evtype='L', outdir='albatross',
                                 wavefiles='test', explosion=True,
                                 overwrite=True)
            # Raises error due to incorrect wavefil formatting
            sfile = eventtosfile(test_cat[0], userID='TEST',
                                 evtype='L', outdir='.',
                                 wavefiles=1234, explosion=True,
                                 overwrite=True)
        with self.assertRaises(IndexError):
            invalid_origin = test_cat[0].copy()
            invalid_origin.origins = []
            sfile = eventtosfile(invalid_origin, userID='TEST',
                                 evtype='L', outdir='.',
                                 wavefiles='test', explosion=True,
                                 overwrite=True)
        with self.assertRaises(ValueError):
            invalid_origin = test_cat[0].copy()
            invalid_origin.origins[0].time = None
            sfile = eventtosfile(invalid_origin, userID='TEST',
                                 evtype='L', outdir='.',
                                 wavefiles='test', explosion=True,
                                 overwrite=True)
        # Write a near empty origin
        valid_origin = test_cat[0].copy()
        valid_origin.origins[0].latitude = None
        valid_origin.origins[0].longitude = None
        valid_origin.origins[0].depth = None
        sfile = eventtosfile(valid_origin, userID='TEST',
                             evtype='L', outdir='.',
                             wavefiles='test', explosion=True,
                             overwrite=True)
        self.assertTrue(os.path.isfile(sfile))
        os.remove(sfile)

    def test_blanksfile(self):
        import os
        from obspy import UTCDateTime
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_',
                                    '2013-09-01-0410-35.DFDPC_024_00')
        sfile = blanksfile(testing_path, 'L', 'TEST', '.', overwrite=True)
        self.assertTrue(os.path.isfile(sfile))
        os.remove(sfile)
        sfile = blanksfile(testing_path, 'L', 'TEST', '.', overwrite=True,
                           evtime=UTCDateTime())
        self.assertTrue(os.path.isfile(sfile))
        os.remove(sfile)
        with self.assertRaises(IOError):
            # No wavefile
            blanksfile('albert', 'L', 'TEST', '.', overwrite=True)
            # No outdir
            blanksfile(testing_path, 'L', 'TEST', 'albert', overwrite=True)
            # USER ID too long
            blanksfile(testing_path, 'L', 'TESTICLE', '.', overwrite=True)
            # Unknown event type
            blanksfile(testing_path, 'U', 'TEST', '.', overwrite=True)

    def test_write_empty(self):
        """
        Function to check that writing a blank event works as it should.
        """
        from obspy.core.event import Event, Origin
        from obspy import UTCDateTime
        import os
        test_event = Event()
        with self.assertRaises(IndexError):
            eventtosfile(test_event, 'TEST', 'L', '.', 'test')
        test_event.origins.append(Origin())
        with self.assertRaises(ValueError):
            eventtosfile(test_event, 'TEST', 'L', '.', 'test')
        test_event.origins[0].time = UTCDateTime()
        test_sfile = eventtosfile(test_event, 'TEST', 'L', '.', 'test')
        self.assertTrue(os.path.isfile(test_sfile))
        os.remove(test_sfile)

    def test_read_empty_header(self):
        """
        Function to check a known issue, empty header info S-file: Bug found \
        by Dominic Evanzia.
        """
        import os
        import numpy as np
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        test_event = readpicks(os.path.join(testing_path, 'Sfile_no_header'))
        self.assertTrue(np.isnan(test_event.origins[0].latitude))
        self.assertTrue(np.isnan(test_event.origins[0].longitude))
        self.assertTrue(np.isnan(test_event.origins[0].depth))

    def test_read_extra_header(self):
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'Sfile_extra_header')
        not_extra_header = os.path.join(os.path.abspath(os.path.
                                                        dirname(__file__)),
                                        'test_data', 'REA', 'TEST_',
                                        '01-0411-15L.S201309')
        test_event = readpicks(testing_path)
        header_event = readpicks(not_extra_header)
        self.assertEqual(test_event.origins[0].time,
                         header_event.origins[0].time)
        self.assertEqual(test_event.origins[0].latitude,
                         header_event.origins[0].latitude)
        self.assertEqual(test_event.origins[0].longitude,
                         header_event.origins[0].longitude)
        self.assertEqual(test_event.origins[0].depth,
                         header_event.origins[0].depth)

    def test_mag_conv(self):
        """Check that we convert magnitudes as we should!"""
        magnitude_map = [('L', 'ML'),
                         ('b', 'mB'),
                         ('s', 'Ms'),
                         ('S', 'MS'),
                         ('W', 'MW'),
                         ('G', 'MbLg'),
                         ('C', 'Mc'),
                         ]
        for magnitude in magnitude_map:
            self.assertEqual(magnitude[0], _evmagtonor(magnitude[1]))
            self.assertEqual(_nortoevmag(magnitude[0]), magnitude[1])

    def test_str_conv(self):
        """Test the simple string conversions."""
        self.assertEqual(_int_conv('albert'), 999)
        self.assertEqual(_float_conv('albert'), 999.0)
        self.assertEqual(_str_conv('albert'), 'albert')
        self.assertEqual(_int_conv('1'), 1)
        self.assertEqual(_float_conv('1'), 1.0)
        self.assertEqual(_str_conv(1), '1')
        self.assertEqual(_int_conv('1.0256'), 999)
        self.assertEqual(_float_conv('1.0256'), 1.0256)
        self.assertEqual(_str_conv(1.0256), '1.0256')

    def test_read_wavename(self):
        from eqcorrscan.utils.sfile_util import readwavename
        import os

        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_',
                                    '19-0926-59L.S201309')
        wavefiles = readwavename(testing_path)
        self.assertEqual(len(wavefiles), 1)

    def test_station_to_seisan(self):
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from eqcorrscan.utils.sfile_util import stationtoseisan

        t1 = UTCDateTime(2012, 3, 26)
        t2 = UTCDateTime(2012, 4, 26)
        client = Client('GEONET')
        bulk = [('NZ', 'FOZ', '*', '*', t1, t2),
                ('NZ', 'JCZ', '*', '*', t1, t2),
                ('NZ', 'WVZ', '*', '*', t1, t2)]
        inventory = client.get_stations_bulk(bulk, level="channel")
        for station in inventory[0]:
            sta_str = stationtoseisan(station)
            self.assertEqual(len(sta_str), 27)

        for station in inventory[0]:
            station.latitude = abs(station.latitude)
            station.longitude = abs(station.longitude)
            sta_str = stationtoseisan(station)
            self.assertEqual(len(sta_str), 27)

        with self.assertRaises(IOError):
            inventory = client.get_stations_bulk(bulk)
            for station in inventory[0]:
                sta_str = stationtoseisan(station)

    def test_read_event(self):
        """Test the wrapper."""
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_',
                                    '01-0411-15L.S201309')
        event = read_event(testing_path)
        self.assertEqual(len(event.origins), 1)

    def test_read_many_events(self):
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'select.out')
        catalog = read_select(testing_path)
        self.assertEqual(len(catalog), 50)

    def test_inaccurate_picks(self):
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'bad_picks.sfile')
        event = readpicks(testing_path)
        pick_string = nordpick(event)
        for pick in pick_string:
            self.assertEqual(len(pick), 80)

    def test_round_len(self):
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'round_len_undef.sfile')
        event = readpicks(testing_path)
        pick_string = nordpick(event)
        for pick in pick_string:
            self.assertEqual(len(pick), 80)


def full_test_event():
    """
    Function to generate a basic, full test event
    """
    from obspy.core.event import Pick, WaveformStreamID, Arrival, Amplitude
    from obspy.core.event import Event, Origin, Magnitude
    from obspy.core.event import EventDescription, CreationInfo
    from obspy import UTCDateTime

    test_event = Event()
    test_event.origins.append(Origin())
    test_event.origins[0].time = UTCDateTime("2012-03-26") + 1
    test_event.event_descriptions.append(EventDescription())
    test_event.event_descriptions[0].text = 'LE'
    test_event.origins[0].latitude = 45.0
    test_event.origins[0].longitude = 25.0
    test_event.origins[0].depth = 15000
    test_event.creation_info = CreationInfo(agency_id='TES')
    test_event.origins[0].time_errors['Time_Residual_RMS'] = 0.01
    test_event.magnitudes.append(Magnitude())
    test_event.magnitudes[0].mag = 0.1
    test_event.magnitudes[0].magnitude_type = 'ML'
    test_event.magnitudes[0].creation_info = CreationInfo('TES')
    test_event.magnitudes[0].origin_id = test_event.origins[0].resource_id
    test_event.magnitudes.append(Magnitude())
    test_event.magnitudes[1].mag = 0.5
    test_event.magnitudes[1].magnitude_type = 'Mc'
    test_event.magnitudes[1].creation_info = CreationInfo('TES')
    test_event.magnitudes[1].origin_id = test_event.origins[0].resource_id
    test_event.magnitudes.append(Magnitude())
    test_event.magnitudes[2].mag = 1.3
    test_event.magnitudes[2].magnitude_type = 'Ms'
    test_event.magnitudes[2].creation_info = CreationInfo('TES')
    test_event.magnitudes[2].origin_id = test_event.origins[0].resource_id

    # Define the test pick
    _waveform_id_1 = WaveformStreamID(station_code='FOZ', channel_code='SHZ',
                                      network_code='NZ')
    _waveform_id_2 = WaveformStreamID(station_code='WTSZ', channel_code='BH1',
                                      network_code=' ')
    # Pick to associate with amplitude
    test_event.picks.append(Pick(waveform_id=_waveform_id_1,
                                 phase_hint='IAML',
                                 polarity='undecidable',
                                 time=UTCDateTime("2012-03-26") + 1.68))
    # Need a second pick for coda
    test_event.picks.append(Pick(waveform_id=_waveform_id_1,
                                 onset='impulsive', phase_hint='PN',
                                 polarity='positive',
                                 time=UTCDateTime("2012-03-26") + 1.68))
    # Unassociated pick
    test_event.picks.append(Pick(waveform_id=_waveform_id_2,
                                 onset='impulsive', phase_hint='SG',
                                 polarity='undecidable',
                                 time=UTCDateTime("2012-03-26") + 1.72))
    # Unassociated pick
    test_event.picks.append(Pick(waveform_id=_waveform_id_2,
                                 onset='impulsive', phase_hint='PN',
                                 polarity='undecidable',
                                 time=UTCDateTime("2012-03-26") + 1.62))
    # Test a generic local magnitude amplitude pick
    test_event.amplitudes.append(Amplitude(generic_amplitude=2.0,
                                           period=0.4,
                                           pick_id=test_event.picks[0].
                                           resource_id,
                                           waveform_id=test_event.picks[0].
                                           waveform_id,
                                           unit='m',
                                           magnitude_hint='Ml'))
    # Test a coda magnitude pick
    test_event.amplitudes.append(Amplitude(generic_amplitude=10,
                                           pick_id=test_event.picks[1].
                                           resource_id,
                                           waveform_id=test_event.picks[1].
                                           waveform_id,
                                           type='END',
                                           category='duration',
                                           unit='s',
                                           magnitude_hint='Mc',
                                           snr=2.3))
    test_event.origins[0].arrivals.append(Arrival(time_weight=2,
                                                  phase=test_event.
                                                  picks[2].
                                                  phase_hint,
                                                  pick_id=test_event.
                                                  picks[2].
                                                  resource_id,
                                                  backazimuth_residual=5,
                                                  time_residual=0.2,
                                                  distance=15,
                                                  azimuth=25))
    test_event.origins[0].arrivals.append(Arrival(time_weight=2,
                                                  phase=test_event.
                                                  picks[3].
                                                  phase_hint,
                                                  pick_id=test_event.
                                                  picks[3].
                                                  resource_id,
                                                  backazimuth_residual=5,
                                                  time_residual=0.2,
                                                  distance=15,
                                                  azimuth=25))
    return test_event

if __name__ == '__main__':
    unittest.main()
