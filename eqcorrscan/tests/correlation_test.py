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
        import cv2
        import os
        from obspy import read
        testing_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                    'test_data')
        image = read(os.path.join(testing_path, 'test_image.ms'))[0].\
            data.astype(np.float32)
        template = read(os.path.join(testing_path, 'test_template.ms'))[0].\
            data.astype(np.float32)
        ccc = cv2.matchTemplate(image=image, templ=template,
                                method=cv2.TM_CCOEFF_NORMED)
        ccc = ccc.flatten()
        expected_ccc = np.load(os.path.join(testing_path, 'test_ccc.npy'))
        matlab_ccc = np.load(os.path.join(testing_path, 'test_ccc_matlab.npy'))
        if (ccc == expected_ccc).all():
            print('Using openCV from source')
            self.assertTrue((ccc == expected_ccc).all())
        elif (ccc == matlab_ccc).all():
            print('Using pre-compiled openCV')
            self.assertTrue((ccc == matlab_ccc).all())
        else:
            print('Well fuck, they are all different!')


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
