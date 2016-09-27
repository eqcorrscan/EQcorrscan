"""
A series of test functions for the core functions in EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest
import os

from eqcorrscan.core import lag_calc
from eqcorrscan.core.template_gen import from_sfile


class TestMethods(unittest.TestCase):
    def test_channel_loop(self):
        """Test the main lag_calc function"""
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_')
        template = from_sfile(sfile=os.path.join(testing_path,
                                                 '21-1412-02L.S201309'),
                              lowcut=5, highcut=15, samp_rate=40,
                              filt_order=4, length=3, swin='all', prepick=0.05)
        detection = from_sfile(sfile=os.path.join(testing_path,
                                                  '21-1759-04L.S201309'),
                               lowcut=5, highcut=15, samp_rate=40,
                               filt_order=4, length=4, swin='all',
                               prepick=0.55)

        i, event = lag_calc._channel_loop(detection=detection,
                                          template=template,
                                          min_cc=0.4, i=0,
                                          detection_id='Tester_01',
                                          interpolate=False)
        matched_traces = []
        detection_stachans = [(tr.stats.station, tr.stats.channel)
                              for tr in detection]
        picked_stachans = [(pick.waveform_id.station_code,
                            pick.waveform_id.channel_code)
                           for pick in event.picks]
        for master_tr in template:
            stachan = (master_tr.stats.station, master_tr.stats.channel)
            if stachan in detection_stachans:
                matched_traces.append(stachan)

        for picked_stachan in picked_stachans:
            self.assertTrue(picked_stachan in matched_traces)

    def test_interpolate(self):
        """Test channel loop with interpolation."""
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'REA', 'TEST_')
        template = from_sfile(sfile=os.path.join(testing_path,
                                                 '21-1412-02L.S201309'),
                              lowcut=5, highcut=15, samp_rate=40,
                              filt_order=4, length=3, swin='all', prepick=0.05)
        detection = from_sfile(sfile=os.path.join(testing_path,
                                                  '21-1759-04L.S201309'),
                               lowcut=5, highcut=15, samp_rate=40,
                               filt_order=4, length=4, swin='all',
                               prepick=0.55)

        i, event = lag_calc._channel_loop(detection=detection,
                                          template=template,
                                          min_cc=0.4, i=0,
                                          detection_id='Tester_01',
                                          interpolate=True)
        matched_traces = []
        detection_stachans = [(tr.stats.station, tr.stats.channel)
                              for tr in detection]
        picked_stachans = [(pick.waveform_id.station_code,
                            pick.waveform_id.channel_code)
                           for pick in event.picks]
        for master_tr in template:
            stachan = (master_tr.stats.station, master_tr.stats.channel)
            if stachan in detection_stachans:
                matched_traces.append(stachan)

        for picked_stachan in picked_stachans:
            self.assertTrue(picked_stachan in matched_traces)

if __name__ == '__main__':
    unittest.main()
