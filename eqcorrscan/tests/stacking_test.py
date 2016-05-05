"""
Functions for testing the utils.stacking functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from eqcorrscan.utils.stacking import linstack, PWS_stack, align_traces
import unittest


class TestStackingMethods(unittest.TestCase):
    def test_linstack(self):
        """Test the utils.stacking.linstack function."""
        # Generate synth data
        import numpy as np
        from obspy import Stream, Trace

        synth = Stream(Trace())
        synth[0].data = np.zeros(200)
        synth[0].data[100] = 1.0
        sine_x = np.arange(0, 10.0, 0.5)
        damped_sine = np.exp(-sine_x) * np.sin(2 * np.pi * sine_x)
        synth[0].data = np.convolve(synth[0].data, damped_sine)
        # Normalize:
        synth[0].data = synth[0].data / synth[0].data.max()
        maximum_synth = synth[0].data.max()
        RMS_max = np.sqrt(np.mean(np.square(synth[0].data)))
        streams = [synth.copy() for i in range(10)]

        stack = linstack(streams, normalize=True)
        # Check normalized amplitude is correct
        self.assertEqual(np.float32(stack[0].data.max()),
                         np.float32(10 * maximum_synth / RMS_max))
        stack = linstack(streams, normalize=False)
        # Check amplitude is preserved
        self.assertEqual(stack[0].data.max(), 10 * maximum_synth)
        # Check length is preserved
        self.assertEqual(len(synth[0].data), len(stack[0].data))

    def test_phase_weighted_stack(self):
        """Test the utils.stacking.PWS_stack."""
        # Generate synth data
        import numpy as np
        from obspy import Stream, Trace

        synth = Stream(Trace())
        synth[0].data = np.zeros(200)
        synth[0].data[100] = 1.0
        sine_x = np.arange(0, 10.0, 0.5)
        damped_sine = np.exp(-sine_x) * np.sin(2 * np.pi * sine_x)
        synth[0].data = np.convolve(synth[0].data, damped_sine)
        # Normalize:
        synth[0].data = synth[0].data / synth[0].data.max()
        # maximum_synth = synth[0].data.max()
        # RMS_max = np.sqrt(np.mean(np.square(synth[0].data)))
        streams = [synth.copy() for i in range(10)]

        stack = PWS_stack(streams, weight=2, normalize=True)
        # Check length is preserved
        self.assertEqual(len(synth[0].data), len(stack[0].data))

    def test_align_traces(self):
        """Test the utils.stacking.align_traces fucntion."""
        # Generate synth data
        import numpy as np
        from obspy import Trace

        synth = Trace()
        synth.data = np.zeros(200)
        synth.data[100] = 1.0
        sine_x = np.arange(0, 10.0, 0.5)
        damped_sine = np.exp(-sine_x) * np.sin(2 * np.pi * sine_x)
        synth.data = np.convolve(synth.data, damped_sine)
        # Normalize:
        synth.data = synth.data / synth.data.max()
        # maximum_synth = synth[0].data.max()
        # RMS_max = np.sqrt(np.mean(np.square(synth[0].data)))
        traces = [synth.copy() for i in range(10)]

        shifts, ccs = align_traces(traces, shift_len=2, master=False)
        for shift in shifts:
            self.assertEqual(shift, 0)

        # Force shifts and check shifts recovered
        shifts_in = np.arange(10)
        original_length = len(traces[0].data)
        for shift, tr in zip(shifts_in, traces):
            pad = np.zeros(shift)
            tr.data = np.concatenate([pad, tr.data])[0:original_length]

        shifts, ccs = align_traces(traces, shift_len=11, master=False)
        for shift_in, shift_out in zip(shifts_in, shifts):
            self.assertEqual(-1 * shift_in, shift_out)

if __name__ == '__main__':
    """
    Run stacking tests
    """
    unittest.main()
