"""
Simple scrip to adapt the python path and run tests on travis using py.test.

Needed to ignore the downloaded eqcorrscan directory, should only use the
installed version: the local version will not have the built, compiled code,
which makes the tests fail.
"""

import sys
import pytest
import os

sys.path = [item for item in sys.path
            if item not in ['.', os.path.dirname(os.path.abspath(__file__))]]
print(sys.path)
pytest.main(['-x', 'eqcorrscan'])
