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
        """Test the utils.stacking.align_traces function."""
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

    def test_known_align(self):
        """Test alignment with a known outcome."""
        from obspy import read
        import os
        import glob
        # from eqcorrscan.utils.stacking import align_traces
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_')
        # testing_path = 'eqcorrscan/tests/test_data/WAV/TEST_/'
        wavefiles = sorted(glob.glob(os.path.join(testing_path, '*')))
        trace_list = []
        for wavfile in wavefiles:
            st = read(wavfile)
            tr = st.select(station='FRAN', channel='SH1')
            if len(tr) == 1:
                tr.detrend('simple').filter('bandpass', freqmin=2, freqmax=20)
                trace_list.append(tr[0])
        shifts, ccs = align_traces(trace_list=trace_list, shift_len=200)
        ccs = [float(str(cc)) for cc in ccs]
        f = open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                              'test_data', 'known_alignment.csv'), 'r')
        known_shifts = [line.rstrip().split(', ') for line in f]
        f.close()
        known_shifts = [(float(a[0]), float(a[1])) for a in known_shifts]
        known_shifts, known_ccs = zip(*known_shifts)
        self.assertEqual(shifts, list(known_shifts))
        ccs = [round(cc, 3) for cc in ccs]
        known_ccs = [round(cc, 3) for cc in known_ccs]
        self.assertEqual(ccs, list(known_ccs))

    def test_known_align_positive(self):
        """Test a known alignment case with forced positive correlation."""
        from obspy import read
        import os
        import glob
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data', 'WAV', 'TEST_')
        wavefiles = sorted(glob.glob(os.path.join(testing_path, '*')))
        trace_list = []
        for wavfile in wavefiles:
            st = read(wavfile)
            tr = st.select(station='FRAN', channel='SH1')
            if len(tr) == 1:
                tr.detrend('simple').filter('bandpass', freqmin=2, freqmax=20)
                trace_list.append(tr[0])
        shifts, ccs = align_traces(trace_list=trace_list, shift_len=200,
                                   positive=True)
        ccs = [float(str(cc)) for cc in ccs]
        f = open(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                              'test_data', 'known_positive_alignment.csv'),
                 'r')
        known_shifts = [line.rstrip().split(', ') for line in f]
        f.close()
        known_shifts = [(float(a[0]), float(a[1])) for a in known_shifts]
        known_shifts, known_ccs = zip(*known_shifts)
        self.assertEqual(shifts, list(known_shifts))
        ccs = [round(cc, 3) for cc in ccs]
        known_ccs = [round(cc, 3) for cc in known_ccs]
        self.assertEqual(ccs, list(known_ccs))

if __name__ == '__main__':
    """
    Run stacking tests
    """
    unittest.main()
