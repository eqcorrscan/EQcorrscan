"""
EQcorrscan's simple logging.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import logging

LOG = {0: logging.debug, 1: logging.info, 2: logging.warning, 3: logging.error,
       4: logging.critical}


def debug_print(string, debug_level, print_level, use_logging=True):
    """
    Print the string if the print_level exceeds the debug_level.

    :type string: str
    :param string: String to print
    :type print_level: int
    :param print_level: Print-level for statement
    :type debug_level: int
    :param debug_level: Output level for function
    :type use_logging: bool
    :param use_logging: Uses logging to output information

    .. rubric:: Example
    >>> debug_print("Albert", 2, 0, use_logging=False)
    >>> debug_print("Norman", 0, 2, use_logging=False)
    Norman
    """
    if not use_logging and print_level > debug_level:
        print(string)
    if use_logging:
        try:
            LOG[print_level](string)
        except IndexError:
            LOG[4](string)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
