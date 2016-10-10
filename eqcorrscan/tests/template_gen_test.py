"""
Functions to test generating templates from SAC data.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import glob
import os
import obspy
import numpy as np
import warnings
import shutil

from obspy import read, UTCDateTime, read_events, Stream
from obspy.clients.fdsn.header import FDSNException
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog, Event, Origin, Pick, WaveformStreamID

from eqcorrscan.core.template_gen import from_sac, _group_events, from_seishub
from eqcorrscan.core.template_gen import from_meta_file, from_client
from eqcorrscan.core.template_gen import multi_template_gen, from_contbase
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
            sac_files = glob.glob(test_files)

            # We currently do not support SAC template generation below version
            # 1.0.0 as before this, SACIO did not fill the reference time,
            # which is needed for defining pick times.  This is usually the
            # trace start time, but isn't always...
            if int(obspy.__version__.split('.')[0]) >= 1:
                template = from_sac(sac_files, lowcut=2.0, highcut=8.0,
                                    samp_rate=samp_rate, filt_order=4,
                                    length=length,
                                    swin='all', prepick=0.1, debug=0,
                                    plot=False)
                for tr in template:
                    self.assertEqual(len(tr.data), length * samp_rate)
            else:
                with self.assertRaises(NotImplementedError):
                    template = from_sac(sac_files, lowcut=2.0, highcut=8.0,
                                        samp_rate=samp_rate, filt_order=4,
                                        length=length,
                                        swin='all', prepick=0.1, debug=0,
                                        plot=False)

    def test_tutorial_template_gen(self):
        """Test template generation from tutorial, uses from_client method.

        Checks that the tutorial generates the templates we expect it to!
        """
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        try:
            mktemplates(plot=False)
        except FDSNException:
            warnings.warn('FDSN error, is server down?')
            return
        for template_no in range(4):
            template = read('tutorial_template_' + str(template_no) + '.ms')
            expected_template = read(os.path.join(testing_path,
                                                  'tutorial_template_' +
                                                  str(template_no) + '.ms'))
            for tr in template:
                expected_tr = expected_template.select(station=tr.stats.
                                                       station,
                                                       channel=tr.stats.
                                                       channel)[0]
                self.assertTrue((expected_tr.data.astype(np.float32) ==
                                 tr.data.astype(np.float32)).all())
            del(template)
            os.remove('tutorial_template_' + str(template_no) + '.ms')

    def test_not_delayed(self):
        """Test the method of template_gen without applying delays to
        channels."""
        cat = get_geonet_events(minlat=-40.98, maxlat=-40.85, minlon=175.4,
                                maxlon=175.5,
                                startdate=UTCDateTime(2016, 5, 1),
                                enddate=UTCDateTime(2016, 5, 2))
        cat = filter_picks(catalog=cat, top_n_picks=5)
        template = from_client(catalog=cat, client_id='GEONET',
                               lowcut=None, highcut=None, samp_rate=100.0,
                               filt_order=4, length=10.0, prepick=0.5,
                               swin='all', process_len=3600,
                               debug=0, plot=False, delayed=False)[0]
        for tr in template:
            tr.stats.starttime.precision = 6
        starttime = template[0].stats.starttime
        length = template[0].stats.npts
        print(template)
        for tr in template:
            self.assertTrue(abs((tr.stats.starttime - starttime)) <=
                            tr.stats.delta)
            self.assertEqual(tr.stats.npts, length)

    def test_download_various_methods(self):
        """
        Will download data from server and store in various databases,
        then create templates using the various methods.
        """
        client = Client('GEONET')
        # get the events
        catalog = Catalog()
        data_stream = client._download('http://quakeml.geonet.org.nz/' +
                                       'quakeml/1.2/2016p008194')
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
        from future import standard_library
        with standard_library.hooks():
            from urllib.request import URLError
        t = UTCDateTime(2009, 9, 3)
        test_cat = Catalog()
        test_cat.append(Event())
        test_cat[0].origins.append(Origin())
        test_cat[0].origins[0].time = t
        test_cat[0].origins[0].latitude = 45
        test_cat[0].origins[0].longitude = 45
        test_cat[0].origins[0].depth = 5000
        test_cat[0].\
            picks.append(Pick(waveform_id=WaveformStreamID(station_code='MANZ',
                                                           channel_code='EHZ',
                                                           network_code='BW'),
                              phase_hint='PG', time=t + 2000))
        test_cat[0].\
            picks.append(Pick(waveform_id=WaveformStreamID(station_code='MANZ',
                                                           channel_code='EHN',
                                                           network_code='BW'),
                              phase_hint='SG', time=t + 2005))
        test_cat[0].\
            picks.append(Pick(waveform_id=WaveformStreamID(station_code='MANZ',
                                                           channel_code='EHE',
                                                           network_code='BW'),
                              phase_hint='SG', time=t + 2005.5))

        test_url = 'http://teide.geophysik.uni-muenchen.de:8080'

        try:
            template = from_seishub(test_cat, url=test_url, lowcut=1.0,
                                    highcut=5.0, samp_rate=20, filt_order=4,
                                    length=3, prepick=0.5, swin='all',
                                    process_len=300)
        except URLError:
            warnings.warn('Timed out connection to seishub')
        if 'template' in locals():
            self.assertEqual(len(template), 3)

    def test_catalog_grouping(self):
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_', '*')
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


if __name__ == '__main__':
    unittest.main()
