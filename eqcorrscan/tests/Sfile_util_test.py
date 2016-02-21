"""
Functions for testing the utils.Sfile_util functions
"""

from eqcorrscan.utils.Sfile_util import eventtoSfile, readwavename, readpicks
from eqcorrscan.utils.Sfile_util import eventtopick, picktoevent


def test_read_write():
    """
    Function to test the read and write capabilities of Sfile_util.
    """
    import os
    from obspy.core.event import Pick, WaveformStreamID, Arrival, Amplitude
    from obspy.core.event import Catalog, Event, Origin, Magnitude
    from obspy.core.event import EventDescription, CreationInfo
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.core.event import read_events
    else:
        from obspy.core.event import readEvents as read_events
    from obspy import UTCDateTime

    # Set-up a test event
    test_event = Event()
    test_event.origins.append(Origin())
    test_event.origins[0].time = UTCDateTime("2012-03-26") + 1
    test_event.event_descriptions.append(EventDescription())
    test_event.event_descriptions[0].text = 'LE'
    test_event.origins[0].latitude = 45.0
    test_event.origins[0].longitude = 25.0
    test_event.origins[0].depth = 15.0
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
    test_event.amplitudes.append(Amplitude(generic_amplitude=2.0, period=0.4,
                                           pick_id=test_event.picks[0].
                                           resource_id,
                                           waveform_id=test_event.picks[0].
                                           waveform_id,
                                           unit='m'))
    test_event.origins[0].arrivals.append(Arrival(time_weight=2,
                                                  phase=test_event.picks[0].
                                                  phase_hint,
                                                  pick_id=test_event.picks[0].
                                                  resource_id,
                                                  backazimuth_residual=5,
                                                  time_residual=0.2,
                                                  distance=15,
                                                  azimuth=25))
    # Add the event to a catalogue which can be used for QuakeML testing
    test_cat = Catalog()
    test_cat += test_event
    # Write the catalog
    test_cat.write("Test_catalog.xml", format='QUAKEML')
    # Read and check
    read_cat = read_events("Test_catalog.xml")
    os.remove("Test_catalog.xml")
    assert read_cat[0].resource_id == test_cat[0].resource_id
    assert read_cat[0].picks == test_cat[0].picks
    assert read_cat[0].origins[0].resource_id ==\
        test_cat[0].origins[0].resource_id
    assert read_cat[0].origins[0].time == test_cat[0].origins[0].time
    # Note that time_residuel_RMS is not a quakeML format
    assert read_cat[0].origins[0].longitude == test_cat[0].origins[0].longitude
    assert read_cat[0].origins[0].latitude == test_cat[0].origins[0].latitude
    assert read_cat[0].origins[0].depth == test_cat[0].origins[0].depth
    assert read_cat[0].magnitudes == test_cat[0].magnitudes
    assert read_cat[0].event_descriptions == test_cat[0].event_descriptions
    assert read_cat[0].amplitudes[0].resource_id ==\
        test_cat[0].amplitudes[0].resource_id
    assert read_cat[0].amplitudes[0].period == test_cat[0].amplitudes[0].period
    assert read_cat[0].amplitudes[0].unit == test_cat[0].amplitudes[0].unit
    assert read_cat[0].amplitudes[0].generic_amplitude ==\
        test_cat[0].amplitudes[0].generic_amplitude
    assert read_cat[0].amplitudes[0].pick_id ==\
        test_cat[0].amplitudes[0].pick_id
    assert read_cat[0].amplitudes[0].waveform_id ==\
        test_cat[0].amplitudes[0].waveform_id

    # Check the read-write s-file functionality
    sfile = eventtoSfile(test_cat[0], userID='TEST', evtype='L', outdir='.',
                         wavefiles='test', explosion=True, overwrite=True)
    del read_cat
    assert readwavename(sfile) == ['test']
    read_cat = Catalog()
    read_cat += readpicks(sfile)
    os.remove(sfile)
    assert read_cat[0].picks[0].time == test_cat[0].picks[0].time
    assert read_cat[0].picks[0].backazimuth == test_cat[0].picks[0].backazimuth
    assert read_cat[0].picks[0].onset == test_cat[0].picks[0].onset
    assert read_cat[0].picks[0].phase_hint == test_cat[0].picks[0].phase_hint
    assert read_cat[0].picks[0].polarity == test_cat[0].picks[0].polarity
    assert read_cat[0].picks[0].waveform_id.station_code ==\
        test_cat[0].picks[0].waveform_id.station_code
    assert read_cat[0].picks[0].waveform_id.channel_code[-1] ==\
        test_cat[0].picks[0].waveform_id.channel_code[-1]
    # assert read_cat[0].origins[0].resource_id ==\
    #     test_cat[0].origins[0].resource_id
    assert read_cat[0].origins[0].time == test_cat[0].origins[0].time
    # Note that time_residuel_RMS is not a quakeML format
    assert read_cat[0].origins[0].longitude == test_cat[0].origins[0].longitude
    assert read_cat[0].origins[0].latitude == test_cat[0].origins[0].latitude
    assert read_cat[0].origins[0].depth == test_cat[0].origins[0].depth
    assert read_cat[0].magnitudes[0].mag == test_cat[0].magnitudes[0].mag
    assert read_cat[0].magnitudes[1].mag == test_cat[0].magnitudes[1].mag
    assert read_cat[0].magnitudes[2].mag == test_cat[0].magnitudes[2].mag
    assert read_cat[0].magnitudes[0].creation_info ==\
        test_cat[0].magnitudes[0].creation_info
    assert read_cat[0].magnitudes[1].creation_info ==\
        test_cat[0].magnitudes[1].creation_info
    assert read_cat[0].magnitudes[2].creation_info ==\
        test_cat[0].magnitudes[2].creation_info
    assert read_cat[0].magnitudes[0].magnitude_type ==\
        test_cat[0].magnitudes[0].magnitude_type
    assert read_cat[0].magnitudes[1].magnitude_type ==\
        test_cat[0].magnitudes[1].magnitude_type
    assert read_cat[0].magnitudes[2].magnitude_type ==\
        test_cat[0].magnitudes[2].magnitude_type
    assert read_cat[0].event_descriptions == test_cat[0].event_descriptions
    # assert read_cat[0].amplitudes[0].resource_id ==\
    #     test_cat[0].amplitudes[0].resource_id
    assert read_cat[0].amplitudes[0].period == test_cat[0].amplitudes[0].period
    assert read_cat[0].amplitudes[0].snr == test_cat[0].amplitudes[0].snr
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
    assert conv_cat[0].picks[0].time == test_cat[0].picks[0].time
    assert conv_cat[0].picks[0].backazimuth == test_cat[0].picks[0].backazimuth
    assert conv_cat[0].picks[0].onset == test_cat[0].picks[0].onset
    assert conv_cat[0].picks[0].phase_hint == test_cat[0].picks[0].phase_hint
    assert conv_cat[0].picks[0].polarity == test_cat[0].picks[0].polarity
    assert conv_cat[0].picks[0].waveform_id.station_code ==\
        test_cat[0].picks[0].waveform_id.station_code
    assert conv_cat[0].picks[0].waveform_id.channel_code[-1] ==\
        test_cat[0].picks[0].waveform_id.channel_code[-1]
    # assert read_cat[0].origins[0].resource_id ==\
    #     test_cat[0].origins[0].resource_id
    assert conv_cat[0].origins[0].time == test_cat[0].origins[0].time
    # Note that time_residuel_RMS is not a quakeML format
    assert conv_cat[0].origins[0].longitude == test_cat[0].origins[0].longitude
    assert conv_cat[0].origins[0].latitude == test_cat[0].origins[0].latitude
    assert conv_cat[0].origins[0].depth == test_cat[0].origins[0].depth
    assert conv_cat[0].magnitudes[0].mag == test_cat[0].magnitudes[0].mag
    assert conv_cat[0].magnitudes[1].mag == test_cat[0].magnitudes[1].mag
    assert conv_cat[0].magnitudes[2].mag == test_cat[0].magnitudes[2].mag
    assert conv_cat[0].magnitudes[0].creation_info ==\
        test_cat[0].magnitudes[0].creation_info
    assert conv_cat[0].magnitudes[1].creation_info ==\
        test_cat[0].magnitudes[1].creation_info
    assert conv_cat[0].magnitudes[2].creation_info ==\
        test_cat[0].magnitudes[2].creation_info
    assert conv_cat[0].magnitudes[0].magnitude_type ==\
        test_cat[0].magnitudes[0].magnitude_type
    assert conv_cat[0].magnitudes[1].magnitude_type ==\
        test_cat[0].magnitudes[1].magnitude_type
    assert conv_cat[0].magnitudes[2].magnitude_type ==\
        test_cat[0].magnitudes[2].magnitude_type
    assert conv_cat[0].event_descriptions == test_cat[0].event_descriptions
    # assert read_cat[0].amplitudes[0].resource_id ==\
    #     test_cat[0].amplitudes[0].resource_id
    assert conv_cat[0].amplitudes[0].period == test_cat[0].amplitudes[0].period
    assert conv_cat[0].amplitudes[0].snr == test_cat[0].amplitudes[0].snr
    return True

if __name__ == '__main__':
    test_read_write()
