"""
Helpers to get library names. Edited for our use from obspy.core.util.libnames

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import os
import ctypes
import logging
from distutils import sysconfig
import importlib.machinery


Logger = logging.getLogger(__name__)


def _load_cdll(name):
    """
    Helper function to load a shared library built during installation
    with ctypes.

    :type name: str
    :param name: Name of the library to load (e.g. 'mseed').
    :rtype: :class:`ctypes.CDLL`
    """
    libdir = os.path.join(os.path.dirname(__file__), 'lib')
    # Try and cope with static FFTW libs
    static_fftw = os.path.join(libdir, 'libfftw3-3.dll')
    static_fftwf = os.path.join(libdir, 'libfftw3f-3.dll')
    try:
        fftw_lib = ctypes.CDLL(str(static_fftw))  # noqa: F841
        fftwf_lib = ctypes.CDLL(str(static_fftwf))  # noqa: F841
    except Exception:
        pass

    # Cope with range of possible extensions for different python versions
    errs = []
    for ext in importlib.machinery.EXTENSION_SUFFIXES:
        libpath = os.path.join(libdir, name + ext)
        try:
            cdll = ctypes.CDLL(str(libpath), mode=ctypes.RTLD_LOCAL)
        except Exception as e:
            msg = 'Could not load shared library "%s".\n\n %s' % (name, str(e))
            errs.append(msg)
            Logger.debug(msg)
        else:
            Logger.debug(f"Loaded library from {libpath}")
            return cdll
    raise ImportError(
        "Could not load shared library {0} due to "
        "errors:\n{1}".format(name, "\n".join(errs)))


if __name__ == '__main__':
    cdll = _load_cdll('libutils')
