"""
:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

# Import mapping to let users ignore the refactoring
from eqcorrscan.core.match_filter.party import Party, read_party  # NOQA
from eqcorrscan.core.match_filter.family import Family  # NOQA
from eqcorrscan.core.match_filter.template import Template  # NOQA
from eqcorrscan.core.match_filter.tribe import Tribe, read_tribe  # NOQA
from eqcorrscan.core.match_filter.detection import (  # NOQA
    Detection, read_detections, get_catalog, write_catalog)  # NOQA
from eqcorrscan.core.match_filter.match_filter import MatchFilterError  # NOQA
from eqcorrscan.core.match_filter.helpers import (  # NOQA
    normxcorr2, extract_from_stream, _spike_test)  # NOQA

CAT_EXT_MAP = {"QUAKEML": "xml", "SC3ML": "xml"}  # , "NORDIC": "out"}
# TODO: add in nordic support once bugs fixed upstream - 1.2.0 Obspy PR #2195


__all__ = [
    'detection', 'family', 'helpers', 'match_filter', 'party', 'template',
    'tribe']


if __name__ == '__main__':
    import doctest

    doctest.testmod(exclude_empty=True)
