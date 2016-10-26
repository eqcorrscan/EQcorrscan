"""Functions to test the parameter set-up
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class TestParameterSetup(unittest.TestCase):
    """Test setup and read/write parameter files"""
    def test_io(self):
        """Test the IO operations on parameters"""
        from eqcorrscan.utils import parameters
        from obspy import UTCDateTime
        import os
        par = parameters.EQcorrscanParameters(['bob'], 2, 8, 3, 20, 0,
                                              UTCDateTime() - 86400,
                                              UTCDateTime(), '.',
                                              'seishub', 4, False, '.',
                                              '.jpg', False, 8.0, 'MAD', 6.0)
        # Write out
        par.write('test_par')
        par_in = parameters.read_parameters('test_par')
        self.assertEqual(par.template_names, par_in.template_names)
        self.assertEqual(par.lowcut, par_in.lowcut)
        self.assertEqual(par.highcut, par_in.highcut)
        self.assertEqual(par.filt_order, par_in.filt_order)
        self.assertEqual(par.samp_rate, par_in.samp_rate)
        self.assertEqual(par.debug, par_in.debug)
        self.assertEqual(par.startdate, par_in.startdate)
        self.assertEqual(par.enddate, par_in.enddate)
        self.assertEqual(par.archive, par_in.archive)
        self.assertEqual(par.arc_type, par_in.arc_type)
        self.assertEqual(par.cores, par_in.cores)
        self.assertEqual(par.plotvar, par_in.plotvar)
        self.assertEqual(par.plotdir, par_in.plotdir)
        self.assertEqual(par.plot_format, par_in.plot_format)
        self.assertEqual(par.tempdir, par_in.tempdir)
        self.assertEqual(par.threshold, par_in.threshold)
        self.assertEqual(par.threshold_type, par_in.threshold_type)
        self.assertEqual(par.trigger_interval, par_in.trigger_interval)
        os.remove('test_par')

if __name__ == '__main__':
    unittest.main()
