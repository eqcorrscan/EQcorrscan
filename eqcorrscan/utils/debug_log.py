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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
