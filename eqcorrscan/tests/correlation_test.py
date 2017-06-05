"""Functions to test how accurate correlations are using different methods.

This should be used to test openCV installs - conda installs give a different
result to installs from source, which is possibly a large issue."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest


class CorrelationTesting(unittest.TestCase):
    def test_cv2(self):
        """Test the cv2 library installed"""
        import numpy as np
        import os
        from obspy import read
        from eqcorrscan.core.match_filter import normxcorr2
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        image = read(os.path.join(testing_path, 'test_image.ms'))[0].\
            data.astype(np.float32)
        template = read(os.path.join(testing_path, 'test_template.ms'))[0].\
            data.astype(np.float32)
        ccc = normxcorr2(image=image, template=template)
        ccc = ccc.flatten()
        expected_ccc = np.load(os.path.join(testing_path, 'test_ccc.npy'))
        matlab_ccc = np.load(os.path.join(testing_path, 'test_ccc_matlab.npy'))
        self.assertTrue(np.allclose(ccc, expected_ccc, atol=0.003))
        self.assertTrue(np.allclose(ccc, matlab_ccc, atol=0.003))


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
