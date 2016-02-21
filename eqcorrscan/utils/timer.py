#!/usr/bin/python
"""
Useful timer for performance testing.

Not licenced as part of EQcorrscan.
"""

import time


class Timer(object):
    """Simple wrapper for timing objects, used for debug."""

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
            print('elapsed time: %f ms' % self.msecs)


if __name__ == "__main__":
    """Doc-test."""
    import doctest
    doctest.testmod()
