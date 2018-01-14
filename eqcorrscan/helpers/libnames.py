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
import warnings
from distutils import sysconfig


def _get_lib_name(lib):
    """
    Helper function to get an architecture and Python version specific library
    filename.
    """
    # append any extension suffix defined by Python for current platform
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    # in principle "EXT_SUFFIX" is what we want.
    # "SO" seems to be deprecated on newer python
    # but: older python seems to have empty "EXT_SUFFIX", so we fall back
    if not ext_suffix:
        try:
            ext_suffix = sysconfig.get_config_var("SO")
        except Exception as e:
            msg = ("Empty 'EXT_SUFFIX' encountered while building CDLL "
                   "filename and fallback to 'SO' variable failed "
                   "(%s)." % str(e))
            warnings.warn(msg)
            pass
    if ext_suffix:
        libname = lib + ext_suffix
    return libname


def _load_cdll(name):
    """
    Helper function to load a shared library built during installation
    with ctypes.

    :type name: str
    :param name: Name of the library to load (e.g. 'mseed').
    :rtype: :class:`ctypes.CDLL`
    """
    # our custom defined part of the extension file name
    libname = _get_lib_name(name)
    libdir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'lib')
    libpath = os.path.join(libdir, libname)
    static_fftw = os.path.join(libdir, 'libfftw3-3.dll')
    static_fftwf = os.path.join(libdir, 'libfftw3f-3.dll')
    try:
        fftw_lib = ctypes.CDLL(str(static_fftw))  # noqa: F841
        fftwf_lib = ctypes.CDLL(str(static_fftwf))  # noqa: F841
    except:
        pass
    try:
        cdll = ctypes.CDLL(str(libpath))
    except Exception as e:
        msg = 'Could not load shared library "%s".\n\n %s' % (libname, str(e))
        raise ImportError(msg)
    return cdll


if __name__ == '__main__':
    cdll = _load_cdll('libutils')
