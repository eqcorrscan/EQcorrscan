"""
Helper functions for seismic data clients.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import logging

from obspy import Stream


Logger = logging.getLogger(__name__)


def get_waveform_client(waveform_client):
    """
    Bind a `get_waveforms_bulk` method to client if it doesn't already have one

    :param waveform_client: Obspy client with a `get_waveforms` method

    :returns: waveform_client with `get_waveforms_bulk`.
    """
    def _get_waveforms_bulk_naive(self, bulk_arg):
        """ naive implementation of get_waveforms_bulk that uses iteration. """
        st = Stream()
        for arg in bulk_arg:
            st += self.get_waveforms(*arg)
        return st

    # add waveform_bulk method dynamically if it doesn't exist already
    if not hasattr(waveform_client, "get_waveforms_bulk"):
        bound_method = _get_waveforms_bulk_naive.__get__(waveform_client)
        setattr(waveform_client, "get_waveforms_bulk", bound_method)

    return waveform_client


if __name__ == "__main__":
    import doctest

    doctest.testmod()
