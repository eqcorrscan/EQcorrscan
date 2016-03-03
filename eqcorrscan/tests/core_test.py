"""
A series of test functions for the core functions in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestCoreMethods(unittest.TestCase):
    # def test_lag_calc():
    #     """
    #     Function to test the capabilites of lag_calc.  When bugs are found \
    #     in the code and fixed the cause of the bug should be emulated here \
    #     to ensure that it is truely fixed.
    #     """
    #     # We need to develop a synthetic dataset...
    #     # This needs to have a list of detection objects, a day of synth
    #     # data, and a series of templates.
    #     return True

    def test_match_filter(self, samp_rate=20.0, debug=0):
        """
        Function to test the capabilities of match_filter and just check that \
        it is working!  Uses synthetic templates and seeded, randomised data.

        :type debug: int
        :param debug: Debug level, higher the number the more output.
        """
        from eqcorrscan.utils import pre_processing
        from eqcorrscan.utils import plotting
        from eqcorrscan.core import match_filter
        from eqcorrscan.utils.synth_seis import generate_synth_data
        from obspy import UTCDateTime
        import string
        # Generate a random dataset
        templates, data, seeds = generate_synth_data(nsta=5, ntemplates=2,
                                                     nseeds=50,
                                                     samp_rate=samp_rate,
                                                     t_length=6.0, max_amp=5.0,
                                                     debug=debug)
        # Notes to the user: If you use more templates you should ensure they
        # are more different, e.g. set the data to have larger moveouts,
        # otherwise similar templates will detect events seeded by another
        # template.
        # Test the pre_processing functions
        data = pre_processing.dayproc(st=data, lowcut=2.0, highcut=8.0,
                                      filt_order=3, samp_rate=samp_rate,
                                      debug=0, starttime=UTCDateTime(0))
        if debug > 0:
            data.plot()
        # Filter the data and the templates
        for template in templates:
            pre_processing.shortproc(st=template, lowcut=2.0, highcut=8.0,
                                     filt_order=3, samp_rate=samp_rate)
            if debug > 0:
                template.plot()
        template_names = list(string.ascii_lowercase)[0:len(templates)]
        detections = match_filter.match_filter(template_names=template_names,
                                               template_list=templates,
                                               st=data, threshold=10.0,
                                               threshold_type='MAD',
                                               trig_int=6.0,
                                               plotvar=False,
                                               plotdir='.',
                                               cores=1,
                                               debug=0)
        # Compare the detections to the seeds
        print('This test made ' + str(len(detections)) + ' detections')
        ktrue = 0
        kfalse = 0
        for detection in detections:
            print(detection.template_name)
            i = template_names.index(detection.template_name)
            t_seeds = seeds[i]
            dtime_samples = int((detection.detect_time - UTCDateTime(0)) *
                                samp_rate)
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
        self.assertTrue(kfalse / ktrue < 0.25)

if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
