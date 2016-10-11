"""
A series of test functions for the core functions in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import unittest
import os
import warnings

from obspy import read, Stream, Trace, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Pick

from eqcorrscan.core import template_gen
from eqcorrscan.utils import pre_processing, catalog_utils
from eqcorrscan.core.match_filter import match_filter, normxcorr2
from eqcorrscan.core.match_filter import _template_loop
from eqcorrscan.tutorials.get_geonet_events import get_geonet_events


class TestCoreMethods(unittest.TestCase):
    def test_perfect_normxcorr2(self):
        """Simple test of normxcorr2 to ensure data are detected
        """
        template = np.random.randn(100).astype(np.float32)
        image = np.zeros(1000).astype(np.float32)
        image[200] = 1.0
        image = np.convolve(template, image)
        ccc = normxcorr2(template, image).astype(np.float16)
        self.assertEqual(ccc.max(), 1.0)

    def test_fail_normxcorr2(self):
        """Ensure if template is nan then return is nan
        """
        template = np.array([np.nan] * 100)
        image = np.zeros(1000)
        image[200] = 1.0
        image = np.convolve(template, image)
        ccc = normxcorr2(template, image)
        self.assertTrue(np.all(ccc == 0.0))

    def test_normal_normxcorr2(self):
        """Check that if match is not perfect correlation max isn't unity
        """
        template = np.random.randn(100) * 10.0
        image = np.zeros(1000)
        image[200] = 1.0
        image = np.convolve(template, image)
        image += np.random.randn(1099)  # Add random noise
        ccc = normxcorr2(template, image)
        self.assertNotEqual(ccc.max(), 1.0)

    def test_set_normxcorr2(self):
        """Check that correlations output are the same irrespective of version.
        """
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        template = read(os.path.join(testing_path, 'test_template.ms'))
        template = template[0].data.astype(np.float32)
        image = read(os.path.join(testing_path, 'test_image.ms'))
        image = image[0].data.astype(np.float32)
        ccc = normxcorr2(template, image)[0]
        expected_ccc = np.load(os.path.join(testing_path, 'test_ccc.npy'))
        # We know that conda installs of openCV give a different results
        # to source built - allow this and allow it to pass.
        self.assertTrue((np.gradient(expected_ccc).round(2) ==
                         np.gradient(ccc).round(2)).all())
        if not (ccc == expected_ccc).all():
            warnings.warn('The expected result was not achieved, ' +
                          'but it has the same shape')

    def test_perfect_template_loop(self):
        """Check that perfect correlations are carried through.
        """
        template = Stream(Trace(np.random.randn(100).astype(np.float32)))
        template[0].stats.station = 'test'
        template[0].stats.channel = 'SZ'
        image = np.zeros(1000).astype(np.float32)
        image[200] = 1.0
        image = np.convolve(image, template[0].data)
        chan = image
        i, ccc = _template_loop(template=template, chan=chan, stream_ind=0)
        self.assertEqual(ccc.astype(np.float16).max(), 1.0)

    def test_false_template_loop(self):
        """Check that perfect correlations are carried through.
        """
        template = Stream(Trace(np.array([np.nan] * 100)))
        template[0].stats.station = 'test'
        template[0].stats.channel = 'SZ'
        image = np.zeros(1000)
        image[200] = 1.0
        image = np.convolve(image, template[0].data)
        chan = image
        i, ccc = _template_loop(template=template, chan=chan, stream_ind=0)
        self.assertTrue(np.all(ccc == 0))

    def test_normal_template_loop(self):
        """Check that perfect correlations are carried through.
        """
        template = Stream(Trace(np.random.randn(100) * 10.0))
        template[0].stats.station = 'test'
        template[0].stats.channel = 'SZ'
        image = np.zeros(1000)
        image[200] = 1.0
        image = np.convolve(image, template[0].data)
        image += np.random.randn(1099)  # Add random noise
        chan = image
        i, ccc = _template_loop(template=template, chan=chan, stream_ind=0)
        self.assertNotEqual(ccc.max(), 1.0)


class TestSynthData(unittest.TestCase):
    def test_debug_range(self):
        """Test range of debug outputs"""
        # debug == 3 fails on travis due to plotting restrictions.
        for debug in range(0, 3):
            print('Testing for debug level=%s' % debug)
            kfalse, ktrue = test_match_filter(debug=debug)
            if ktrue > 0:
                self.assertTrue(kfalse / ktrue < 0.25)
            else:
                # Randomised data occasionally yields 0 detections
                kfalse, ktrue = test_match_filter(debug=debug)
                self.assertTrue(kfalse / ktrue < 0.25)

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


class TestGeoNetCase(unittest.TestCase):
    def setUp(self):
        client = Client('GEONET')
        t1 = UTCDateTime(2016, 9, 4)
        t2 = t1 + 86400
        catalog = get_geonet_events(startdate=t1, enddate=t2, minmag=4,
                                    minlat=-49, maxlat=-35, minlon=175.0,
                                    maxlon=185.0)
        catalog = catalog_utils.filter_picks(catalog, channels=['EHZ'],
                                             top_n_picks=5)
        for event in catalog:
            extra_pick = Pick()
            extra_pick.phase_hint = 'S'
            extra_pick.time = event.picks[0].time + 10
            extra_pick.waveform_id = event.picks[0].waveform_id
            event.picks.append(extra_pick)
        self.templates = template_gen.from_client(catalog=catalog,
                                                  client_id='GEONET',
                                                  lowcut=2.0, highcut=9.0,
                                                  samp_rate=50.0, filt_order=4,
                                                  length=3.0, prepick=0.15,
                                                  swin='all', process_len=3600)
        # Download and process the day-long data
        bulk_info = [(tr.stats.network, tr.stats.station, '*',
                      tr.stats.channel[0] + 'H' + tr.stats.channel[1],
                      t1 + (4 * 3600), t1 + (5 * 3600))
                     for tr in self.templates[0]]
        # Just downloading an hour of data
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        self.st = pre_processing.shortproc(st, lowcut=2.0, highcut=9.0,
                                           filt_order=4, samp_rate=50.0,
                                           debug=0, num_cores=1)
        st.trim(t1 + (4 * 3600), t1 + (5 * 3600))
        self.template_names = [str(template[0].stats.starttime)
                               for template in self.templates]

    def test_duplicate_channels_in_template(self):
        """
        Test using a template with duplicate channels.
        """
        detections = match_filter(template_names=self.template_names,
                                  template_list=self.templates, st=self.st,
                                  threshold=8.0, threshold_type='MAD',
                                  trig_int=6.0, plotvar=False, plotdir='.',
                                  cores=4)
        self.assertEqual(len(detections), 1)
        self.assertEqual(detections[0].no_chans, 6)


class TestNCEDCCases(unittest.TestCase):
    def setUp(self):
        print('\t\t\t Downloading data')
        client = Client('NCEDC')
        t1 = UTCDateTime(2004, 9, 28, 17)
        t2 = t1 + 3600
        catalog = client.get_events(starttime=t1, endtime=t2,
                                    minmagnitude=4,
                                    minlatitude=35.7, maxlatitude=36.1,
                                    minlongitude=-120.6,
                                    maxlongitude=-120.2,
                                    includearrivals=True)
        catalog = catalog_utils.filter_picks(catalog, channels=['EHZ'],
                                             top_n_picks=5)
        self.templates = template_gen.from_client(catalog=catalog,
                                                  client_id='NCEDC',
                                                  lowcut=2.0, highcut=9.0,
                                                  samp_rate=50.0, filt_order=4,
                                                  length=3.0, prepick=0.15,
                                                  swin='all', process_len=3600)
        for template in self.templates:
            template.sort()
        # Download and process the day-long data
        bulk_info = [(tr.stats.network, tr.stats.station, '*',
                      tr.stats.channel[0] + 'H' + tr.stats.channel[1],
                      t1, t1 + 3600)
                     for tr in self.templates[0]]
        # Just downloading an hour of data
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        self.st = pre_processing.shortproc(st, lowcut=2.0, highcut=9.0,
                                           filt_order=4, samp_rate=50.0,
                                           debug=0, num_cores=4)
        self.template_names = [str(template[0].stats.starttime)
                               for template in self.templates]

    def test_detection_extraction(self):
        # Test outputting the streams works
        detections, detection_streams = \
            match_filter(template_names=self.template_names,
                         template_list=self.templates, st=self.st,
                         threshold=8.0, threshold_type='MAD',
                         trig_int=6.0, plotvar=False, plotdir='.',
                         cores=4, extract_detections=True)
        self.assertEqual(len(detections), 4)
        self.assertEqual(len(detection_streams), len(detections))


def test_match_filter(debug=0, plotvar=False,
                      extract_detections=False, threshold_type='MAD',
                      threshold=10, template_excess=False,
                      stream_excess=False):
    """
    Function to test the capabilities of match_filter and just check that \
    it is working!  Uses synthetic templates and seeded, randomised data.

    :type samp_rate: float
    :param samp_rate: Sampling rate in Hz to use
    :type debug: int
    :param debug: Debug level, higher the number the more output.
    """
    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils import plotting
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
        pre_processing.shortproc(st=template, lowcut=1.0, highcut=4.0,
                                 filt_order=3, samp_rate=10.0)
    template_names = list(string.ascii_lowercase)[0:len(templates)]
    detections =\
        match_filter(template_names=template_names,
                     template_list=templates, st=data, threshold=threshold,
                     threshold_type=threshold_type, trig_int=6.0,
                     plotvar=plotvar, plotdir='.', cores=1, debug=debug,
                     output_cat=False, extract_detections=extract_detections)
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
            plotting.detection_multiplot(data, template, times)
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
