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


DEBUG_MAP = {0: logging.CRITICAL, 1: logging.ERROR, 2: logging.WARNING,
             3: logging.INFO, 4: logging.DEBUG}


def debug_print(string, debug_level, print_level):
    """
    Print the string if the print_level exceeds the debug_level.

    :type string: str
    :param string: String to print
    :type print_level: int
    :param print_level: Print-level for statement
    :type debug_level: int
    :param debug_level: Output level for function

    .. rubric:: Example
    >>> debug_print("Albert", 2, 0)
    >>> debug_print("Norman", 0, 2)
    Norman
    """
    if print_level > debug_level:
        print(string)


def debug_logger(name, debug_level):
    """
    Create a logger for a function at a certain debug output level.

    :type name: str
    :param name: Function or object name
    :type debug_level: int
    :param debug_level: Debug level - corresponds to level of logger

    :return: logging.logger

    .. rubric:: Example
    >>> logger = debug_logger("Tester", debug_level=3)
    >>> logger.info("There is no output from info at this level")
    >>> logger.warning("But it will log warnings") # doctest:+SKIP
    2...- Tester - WARNING - But it will log warnings
    """
    logger = logging.getLogger(name=name)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(DEBUG_MAP[debug_level])
    logger.addHandler(ch)
    return logger


if __name__ == '__main__':
    import doctest
    doctest.testmod()
