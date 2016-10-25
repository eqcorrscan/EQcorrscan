"""
Test functions for the eqcorrscan.utils.trigger module
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestTriggerMethods(unittest.TestCase):
    def test_parameter_read_write(self):
        """Test parameter writing and reading."""
        from eqcorrscan.utils.trigger import read_trigger_parameters
        from eqcorrscan.utils.trigger import TriggerParameters
        import os

        test_parameters = TriggerParameters()
        test_parameters.write('test_parameters.par')
        check_parameters = read_trigger_parameters('test_parameters.par')
        self.assertEqual(test_parameters, check_parameters[0])
        os.remove('test_parameters.par')

    def test_channel_loop(self):
        """Test trigger generation in internal loop."""
        import numpy as np
        from eqcorrscan.utils.trigger import _channel_loop
        from eqcorrscan.utils.trigger import TriggerParameters
        from obspy import Trace

        parameters = [TriggerParameters({'station': 'TEST',
                                         'channel': 'SHZ',
                                         'sta_len': 0.3,
                                         'lta_len': 10.0,
                                         'thr_on': 10,
                                         'thr_off': 3,
                                         'lowcut': 2,
                                         'highcut': 20})]
        tr = Trace()
        tr.data = np.random.randn(2000)
        tr.data[1000:1010] = [100, -80, 70, -65, 60, -52, 45, -30, 15, 5]
        tr.stats.sampling_rate = 100
        tr.stats.station = parameters[0]['station']
        tr.stats.channel = parameters[0]['channel']
        # Test without despike
        triggers = _channel_loop(tr=tr, parameters=parameters,
                                 max_trigger_length=100,
                                 despike=False, debug=0)
        self.assertEqual(len(triggers), 1)
        # Test with despike
        triggers = _channel_loop(tr=tr, parameters=parameters,
                                 max_trigger_length=100,
                                 despike=True, debug=0)
        self.assertEqual(len(triggers), 1)
        # Test with no filter
        parameters[0]['lowcut'] = None
        parameters[0]['highcut'] = None
        triggers = _channel_loop(tr=tr, parameters=parameters,
                                 max_trigger_length=100,
                                 despike=False, debug=0)
        self.assertEqual(len(triggers), 1)
        # Test with lowpass
        parameters[0]['highcut'] = 20
        triggers = _channel_loop(tr=tr, parameters=parameters,
                                 max_trigger_length=100,
                                 despike=False, debug=0)
        self.assertEqual(len(triggers), 1)
        # Test with highpass
        parameters[0]['highcut'] = None
        parameters[0]['lowcut'] = 2
        triggers = _channel_loop(tr=tr, parameters=parameters,
                                 max_trigger_length=100,
                                 despike=False, debug=0)
        self.assertEqual(len(triggers), 1)

    def test_main_trigger_routine(self):
        """Test the network_trigger function."""
        from eqcorrscan.utils.trigger import network_trigger, TriggerParameters
        from obspy import read
        import os
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_')
        st = read(os.path.join(testing_path,
                               '2013-09-01-0410-36.DFDPC_027_00'))
        parameters = []
        for tr in st:
            parameters.append(TriggerParameters({'station': tr.stats.station,
                                                 'channel': tr.stats.channel,
                                                 'sta_len': 0.5,
                                                 'lta_len': 10.0,
                                                 'thr_on': 10.0,
                                                 'thr_off': 3.0,
                                                 'lowcut': 2.0,
                                                 'highcut': 15.0}))
        triggers = network_trigger(st=st, parameters=parameters,
                                   thr_coincidence_sum=5, moveout=30,
                                   max_trigger_length=60, despike=False)
        self.assertEqual(len(triggers), 1)


if __name__ == '__main__':
    unittest.main()
