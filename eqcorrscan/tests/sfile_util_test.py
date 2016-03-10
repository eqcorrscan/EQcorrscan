"""
Functions for testing the utils.sfile_util functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.utils.sfile_util import eventtosfile, readwavename, readpicks
from eqcorrscan.utils.sfile_util import eventtopick, picktoevent
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
            client = Client(event_info[0])
            if event_info[0] == 'GEONET':
                try:
                    data_stream = client.\
                        _download('http://quakeml.geonet.org.nz/' +
                                  'quakeml/1.2/' + event_info[1])
                    data_stream.seek(0, 0)
                    event = read_events(data_stream, format="quakeml")
                    data_stream.close()
                except FDSNException:
                    warnings.warn('FDSNException')
                    continue
            else:
                try:
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
        test_event = basic_test_event()
        # Add the event to a catalogue which can be used for QuakeML testing
        test_cat = Catalog()
        test_cat += test_event
        # Write the catalog
        test_cat.write("Test_catalog.xml", format='QUAKEML')
        # Read and check
        read_cat = read_events("Test_catalog.xml")
        os.remove("Test_catalog.xml")
        self.assertEqual(read_cat[0].resource_id, test_cat[0].resource_id)
        self.assertEqual(read_cat[0].picks, test_cat[0].picks)
        self.assertEqual(read_cat[0].origins[0].resource_id,
                         test_cat[0].origins[0].resource_id)
        self.assertEqual(read_cat[0].origins[0].time,
                         test_cat[0].origins[0].time)
        # Note that time_residuel_RMS is not a quakeML format
        self.assertEqual(read_cat[0].origins[0].longitude,
                         test_cat[0].origins[0].longitude)
        self.assertEqual(read_cat[0].origins[0].latitude,
                         test_cat[0].origins[0].latitude)
        self.assertEqual(read_cat[0].origins[0].depth,
                         test_cat[0].origins[0].depth)
        self.assertEqual(read_cat[0].magnitudes, test_cat[0].magnitudes)
        self.assertEqual(read_cat[0].event_descriptions,
                         test_cat[0].event_descriptions)
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

        # Check the read-write s-file functionality
        sfile = eventtosfile(test_cat[0], userID='TEST',
                             evtype='L', outdir='.',
                             wavefiles='test', explosion=True, overwrite=True)
        del read_cat
        self.assertEqual(readwavename(sfile), ['test'])
        read_cat = Catalog()
        read_cat += readpicks(sfile)
        os.remove(sfile)
        self.assertEqual(read_cat[0].picks[0].time,
                         test_cat[0].picks[0].time)
        self.assertEqual(read_cat[0].picks[0].backazimuth,
                         test_cat[0].picks[0].backazimuth)
        self.assertEqual(read_cat[0].picks[0].onset,
                         test_cat[0].picks[0].onset)
        self.assertEqual(read_cat[0].picks[0].phase_hint,
                         test_cat[0].picks[0].phase_hint)
        self.assertEqual(read_cat[0].picks[0].polarity,
                         test_cat[0].picks[0].polarity)
        self.assertEqual(read_cat[0].picks[0].waveform_id.station_code,
                         test_cat[0].picks[0].waveform_id.station_code)
        self.assertEqual(read_cat[0].picks[0].waveform_id.channel_code[-1],
                         test_cat[0].picks[0].waveform_id.channel_code[-1])
        # assert read_cat[0].origins[0].resource_id ==\
        #     test_cat[0].origins[0].resource_id
        self.assertEqual(read_cat[0].origins[0].time,
                         test_cat[0].origins[0].time)
        # Note that time_residuel_RMS is not a quakeML format
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
        del read_cat
        # assert read_cat[0].amplitudes[0].pick_id ==\
        #     test_cat[0].amplitudes[0].pick_id
        # assert read_cat[0].amplitudes[0].waveform_id ==\
        #     test_cat[0].amplitudes[0].waveform_id

        # Test the wrappers for PICK and EVENTINFO classes
        picks, evinfo = eventtopick(test_cat)
        # Test the conversion back
        conv_cat = Catalog()
        conv_cat.append(picktoevent(evinfo, picks))
        self.assertEqual(conv_cat[0].picks[0].time, test_cat[0].picks[0].time)
        self.assertEqual(conv_cat[0].picks[0].backazimuth,
                         test_cat[0].picks[0].backazimuth)
        self.assertEqual(conv_cat[0].picks[0].onset,
                         test_cat[0].picks[0].onset)
        self.assertEqual(conv_cat[0].picks[0].phase_hint,
                         test_cat[0].picks[0].phase_hint)
        self.assertEqual(conv_cat[0].picks[0].polarity,
                         test_cat[0].picks[0].polarity)
        self.assertEqual(conv_cat[0].picks[0].waveform_id.station_code,
                         test_cat[0].picks[0].waveform_id.station_code)
        self.assertEqual(conv_cat[0].picks[0].waveform_id.channel_code[-1],
                         test_cat[0].picks[0].waveform_id.channel_code[-1])
        # self.assertEqual(read_cat[0].origins[0].resource_id,
        #                  test_cat[0].origins[0].resource_id)
        self.assertEqual(conv_cat[0].origins[0].time,
                         test_cat[0].origins[0].time)
        # Note that time_residuel_RMS is not a quakeML format
        self.assertEqual(conv_cat[0].origins[0].longitude,
                         test_cat[0].origins[0].longitude)
        self.assertEqual(conv_cat[0].origins[0].latitude,
                         test_cat[0].origins[0].latitude)
        self.assertEqual(conv_cat[0].origins[0].depth,
                         test_cat[0].origins[0].depth)
        self.assertEqual(conv_cat[0].magnitudes[0].mag,
                         test_cat[0].magnitudes[0].mag)
        self.assertEqual(conv_cat[0].magnitudes[1].mag,
                         test_cat[0].magnitudes[1].mag)
        self.assertEqual(conv_cat[0].magnitudes[2].mag,
                         test_cat[0].magnitudes[2].mag)
        self.assertEqual(conv_cat[0].magnitudes[0].creation_info,
                         test_cat[0].magnitudes[0].creation_info)
        self.assertEqual(conv_cat[0].magnitudes[1].creation_info,
                         test_cat[0].magnitudes[1].creation_info)
        self.assertEqual(conv_cat[0].magnitudes[2].creation_info,
                         test_cat[0].magnitudes[2].creation_info)
        self.assertEqual(conv_cat[0].magnitudes[0].magnitude_type,
                         test_cat[0].magnitudes[0].magnitude_type)
        self.assertEqual(conv_cat[0].magnitudes[1].magnitude_type,
                         test_cat[0].magnitudes[1].magnitude_type)
        self.assertEqual(conv_cat[0].magnitudes[2].magnitude_type,
                         test_cat[0].magnitudes[2].magnitude_type)
        self.assertEqual(conv_cat[0].event_descriptions,
                         test_cat[0].event_descriptions)
        # self.assertEqual(read_cat[0].amplitudes[0].resource_id,
        #                  test_cat[0].amplitudes[0].resource_id)
        self.assertEqual(conv_cat[0].amplitudes[0].period,
                         test_cat[0].amplitudes[0].period)
        self.assertEqual(conv_cat[0].amplitudes[0].snr,
                         test_cat[0].amplitudes[0].snr)

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


def basic_test_event():
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
    _waveform_id = WaveformStreamID(station_code='FOZ', channel_code='SHZ',
                                    network_code='NZ')
    test_event.picks.append(Pick(waveform_id=_waveform_id,
                                 onset='impulsive', phase_hint='PN',
                                 polarity='positive',
                                 time=UTCDateTime("2012-03-26") + 1.68,
                                 horizontal_slowness=12, backazimuth=20))
    test_event.amplitudes.append(Amplitude(generic_amplitude=2.0,
                                           period=0.4,
                                           pick_id=test_event.picks[0].
                                           resource_id,
                                           waveform_id=test_event.picks[0].
                                           waveform_id,
                                           unit='m'))
    test_event.origins[0].arrivals.append(Arrival(time_weight=2,
                                                  phase=test_event.
                                                  picks[0].
                                                  phase_hint,
                                                  pick_id=test_event.
                                                  picks[0].
                                                  resource_id,
                                                  backazimuth_residual=5,
                                                  time_residual=0.2,
                                                  distance=15,
                                                  azimuth=25))
    return test_event

if __name__ == '__main__':
    unittest.main()
