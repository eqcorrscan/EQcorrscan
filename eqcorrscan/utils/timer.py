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


def time_func(func, name, *args, **kwargs):
    """ call a func with args and kwargs, print name of func and how
    long it took. """
    tic = time.time()
    out = func(*args, **kwargs)
    toc = time.time()
    print('%s took %0.2f seconds' % (name, toc - tic))
    return out

if __name__ == "__main__":
    """Doc-test."""
    import doctest
    doctest.testmod()
