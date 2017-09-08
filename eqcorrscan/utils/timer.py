"""
Useful timer for performance testing.

Not licenced as part of EQcorrscan.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import time


class Timer(object):
    """
    Simple wrapper for timing objects, used for debug.

    .. rubric:: Example

    >>> from time import sleep
    >>> with Timer() as t:
    ...     sleep(0.1)
    >>> print("%.1f" % t.secs) # doctest:+SKIP
    0.1
    >>> with Timer(verbose=True) as t:
    ...     sleep(0.1)  # doctest:+ELLIPSIS
    elapsed time: ...
    """

    def __init__(self, verbose=False):
        """Create timer object."""
        self.verbose = verbose

    def __enter__(self):
        """Enter timer."""
        self.start = time.time()
        return self

    def __exit__(self, *args):
        """Exit timer."""
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print('elapsed time: %.0f ms' % self.msecs)


if __name__ == "__main__":
    """Doc-test."""
    import doctest
    doctest.testmod()
