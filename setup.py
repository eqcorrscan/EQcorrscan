
from setuptools import setup, Command, Extension
from setuptools.command.build_ext import build_ext
using_setuptools = True

from distutils.ccompiler import get_default_compiler
from pkg_resources import get_build_platform

import os
import sys
import shutil
import sysconfig
import glob

with open("eqcorrscan/__init__.py", "r") as init_file:
    version_line = [line for line in init_file
                    if '__version__' in line][0]
VERSION = version_line.split()[-1].split("'")[1]

# Check if we are on RTD and don't build extensions if we are.
READ_THE_DOCS = os.environ.get('READTHEDOCS', None) == 'True'

long_description = '''
EQcorrscan: repeating and near-repeating earthquake detection and analysis
in Python.  Open-source routines for: systematic template
creation, matched-filter detection, subspace detection, brightness detection,
clustering of seismic events, magnitude calculation by singular value
decomposition, and more!
'''

# Get a list of all the scripts not to be installed
scriptfiles = glob.glob('eqcorrscan/tutorials/*.py')
scriptfiles += glob.glob('eqcorrscan/scripts/*.py')

if READ_THE_DOCS:
    try:
        environ = os.environb
    except AttributeError:
        environ = os.environ

    environ[b"CC"] = b"x86_64-linux-gnu-gcc"
    environ[b"LD"] = b"x86_64-linux-gnu-ld"
    environ[b"AR"] = b"x86_64-linux-gnu-ar"


def get_package_data():
    """
    Get the library files for appveyor - this shouldn't make any difference
    for install with system libraries.
    """
    from pkg_resources import get_build_platform

    package_data = {}

    if get_build_platform() in ('win32', 'win-amd64'):
        package_data['eqcorrscan.utils.lib'] = [
            'libfftw3-3.dll', 'libfftw3f-3.dll', 'libfftw3l-3.dll']

    return package_data


def get_package_dir():
    from pkg_resources import get_build_platform

    package_dir = {}
    if get_build_platform() in ('win32', 'win-amd64'):
        package_dir['eqcorrscan.utils.lib'] = os.path.join(
            'eqcorrscan', 'utils', 'lib')

    return package_dir


def get_include_dirs():
    import numpy
    from pkg_resources import get_build_platform

    include_dirs = [os.path.join(os.getcwd(), 'include'),
                    os.path.join(os.getcwd(), 'eqcorrscan', 'utils', 'src'),
                    numpy.get_include(),
                    os.path.join(sys.prefix, 'include')]

    if get_build_platform() in ('win32', 'win-amd64'):
        # Add the Library dir
        include_dirs.append(os.path.join(sys.prefix, 'Library', 'include'))

    if get_build_platform().startswith('freebsd'):
        include_dirs.append('/usr/local/include')

    return include_dirs


def get_library_dirs():
    from pkg_resources import get_build_platform

    library_dirs = []
    if get_build_platform() in ('win32', 'win-amd64'):
        library_dirs.append(os.path.join(os.getcwd(), 'eqcorrscan', 'utils',
                                         'lib'))
        library_dirs.append(os.path.join(sys.prefix, 'lib'))
        library_dirs.append(os.path.join(sys.prefix, 'Library', 'lib'))

    library_dirs.append(os.path.join(sys.prefix, 'lib'))
    if get_build_platform().startswith('freebsd'):
        library_dirs.append('/usr/local/lib')

    return library_dirs


def get_mkl():
    from pkg_resources import get_build_platform

    mkl_found = False
    # TODO: not sure about windows so ignoring for now
    if not get_build_platform() in ('win32', 'win-amd64'):
        # look for MKL
        mklroot = os.getenv("MKLROOT")
        if mklroot is not None and os.path.isdir(mklroot):
            if os.path.isdir(os.path.join(mklroot, "mkl")):
                mklroot = os.path.join(mklroot, "mkl")

            # look for MKL lib
            _libdir = os.path.join(mklroot, "lib")
            if os.path.isdir(os.path.join(_libdir, "intel64")):
                _libdir = os.path.join(_libdir, "intel64")
            libs = glob.glob(os.path.join(_libdir, "libmkl_rt.*"))
            if len(libs) == 1:
                _lib = os.path.splitext(os.path.basename(libs[0]))[0][3:]
                mkl_lib = [_lib]
                mkl_libdir = [_libdir]

                # look for include
                if os.path.exists(os.path.join(mklroot, "include", "fftw",
                                               "fftw3.h")):
                    mkl_inc = [os.path.join(mklroot, "include", "fftw")]
                    mkl_found = True

        # check for conda installed MKL too (doesn't set MKLROOT)
        # this assume we are in a conda env (not root)
        if not mkl_found:
            conda = os.getenv("CONDA_PREFIX")
            print(conda)
            if conda is not None:
                # look for MKL lib
                libs = glob.glob(os.path.join(conda, "lib", "libmkl_rt.*"))
                if len(libs) == 1:
                    _lib = os.path.splitext(os.path.basename(libs[0]))[0][3:]
                    mkl_lib = [_lib]
                    mkl_libdir = [os.path.join(conda, "lib")]

                    # look for include too
                    if os.path.exists(os.path.join(conda, "include", "fftw",
                                                   "fftw3.h")):
                        mkl_inc = [os.path.join(conda, "include", "fftw")]
                        mkl_found = True
                    elif os.path.exists(os.path.join(conda, "include",
                                                     "fftw3.h")):
                        mkl_inc = [os.path.join(conda, "include")]
                        mkl_found = True
    if mkl_found:
        print("Found MKL:")
        print("  MKL includes:", mkl_inc)
        print("  MKL lib dirs:", mkl_libdir)
        print("  MKL libs:", mkl_lib)
        return mkl_inc, mkl_libdir, mkl_lib
    else:
        return None


def get_libraries():
    from pkg_resources import get_build_platform

    if get_build_platform() in ('win32', 'win-amd64'):
        # libraries = ['libfftw3-3', 'libfftw3f-3']
        libraries = ['fftw3', 'fftw3f']
    else:
        libraries = ['fftw3', 'fftw3_threads', 'fftw3f', 'fftw3f_threads']

    return libraries


def export_symbols(*path):
    """
    Required for windows systems - functions defined in libutils.def.
    """
    lines = open(os.path.join(*path), 'r').readlines()[2:]
    return [s.strip() for s in lines if s.strip() != '']


def get_extensions(no_mkl=False):
    from distutils.extension import Extension

    if READ_THE_DOCS:
        return []
    # will use static linking if STATIC_FFTW_DIR defined
    static_fftw_path = os.environ.get('STATIC_FFTW_DIR', None)
    link_static_fftw = static_fftw_path is not None

    common_extension_args = {
        'include_dirs': get_include_dirs(),
        'library_dirs': get_library_dirs()}

    sources = [os.path.join('eqcorrscan', 'utils', 'src', 'multi_corr.c'),
               os.path.join('eqcorrscan', 'utils', 'src', 'time_corr.c'),
               os.path.join('eqcorrscan', 'utils', 'src', 'find_peaks.c'),
               os.path.join('eqcorrscan', 'utils', 'src',
                            'distance_cluster.c')]
    exp_symbols = export_symbols("eqcorrscan/utils/src/libutils.def")

    if get_build_platform() not in ('win32', 'win-amd64'):
        if get_build_platform().startswith('freebsd'):
            # Clang uses libomp, not libgomp
            extra_link_args = ['-lm', '-lomp']
        else:
            extra_link_args = ['-lm', '-lgomp']
        extra_compile_args = ['-fopenmp']
        if all(arch not in get_build_platform()
               for arch in ['arm', 'aarch']):
            extra_compile_args.extend(['-msse2', '-ftree-vectorize'])
    else:
        extra_link_args = []
        extra_compile_args = ['/openmp', '/TP']

    libraries = get_libraries()
    if link_static_fftw:
        if get_build_platform() in ('win32', 'win-amd64'):
            lib_pre = ''
            lib_ext = '.lib'
        else:
            lib_pre = 'lib'
            lib_ext = '.a'
        for lib in libraries:
            extra_link_args.append(
                os.path.join(static_fftw_path, lib_pre + lib + lib_ext))

        common_extension_args['extra_link_args'] = extra_link_args
        common_extension_args['libraries'] = []
        common_extension_args['extra_compile_args'] = extra_compile_args
        common_extension_args['export_symbols'] = exp_symbols
    else:
        # otherwise we use dynamic libraries
        common_extension_args['extra_link_args'] = extra_link_args
        common_extension_args['extra_compile_args'] = extra_compile_args
        common_extension_args['export_symbols'] = exp_symbols
        if no_mkl:
            mkl = None
        else:
            mkl = get_mkl()
        if mkl is not None:
            # use MKL if we have it
            common_extension_args['include_dirs'].extend(mkl[0])
            common_extension_args['library_dirs'].extend(mkl[1])
            common_extension_args['libraries'] = mkl[2]
        else:
            common_extension_args['libraries'] = libraries
    ext_modules = [
        Extension('eqcorrscan.utils.lib.libutils', sources=sources,
                  **common_extension_args)]
    return ext_modules


class CustomBuildExt(build_ext):
    def finalize_options(self):

        build_ext.finalize_options(self)

        if self.compiler is None:
            compiler = get_default_compiler()
        else:
            compiler = self.compiler

        cfg_vars = sysconfig.get_config_vars()
        # Hack around OSX setting a -m flag
        if "macosx" in get_build_platform() and "CFLAGS" in cfg_vars:
            print("System C-flags:")
            print(cfg_vars["CFLAGS"])
            cflags = []
            for flag in cfg_vars["CFLAGS"].split():
                if flag in ["-m", "-isysroot"]:
                    continue
                # Remove sdk links
                if flag.endswith(".sdk"):
                    continue
                cflags.append(flag)
            cfg_vars["CFLAGS"] = " ".join(cflags)
            print("Editted C-flags:")
            print(cfg_vars["CFLAGS"])
        # Remove unsupported C-flags
        unsupported_flags = [
            "-fuse-linker-plugin", "-ffat-lto-objects", "-flto-partition=none"]
        for key in ["CFLAGS", "LDFLAGS", "LDSHARED"]:
            if key in cfg_vars:
                print("System {0}:".format(key))
                print(cfg_vars[key])
                flags = []
                for flag in cfg_vars[key].split():
                    if flag in unsupported_flags:
                        continue
                    flags.append(flag)
                cfg_vars[key] = " ".join(flags)
                print("Editted {0}:".format(key))
                print(cfg_vars[key])

        if compiler == 'msvc':
            # Add msvc specific hacks

            # Sort linking issues with init exported symbols
            def _get_export_symbols(self, ext):
                return ext.export_symbols

            build_ext.get_export_symbols = _get_export_symbols

            if (sys.version_info.major, sys.version_info.minor) < (3, 3):
                # The check above is a nasty hack. We're using the python
                # version as a proxy for the MSVC version. 2008 doesn't
                # have stdint.h, so is needed. 2010 does.
                #
                # We need to add the path to msvc includes

                msvc_2008_path = (
                    os.path.join(os.getcwd(), 'include', 'msvc_2008'))

                if self.include_dirs is not None:
                    self.include_dirs.append(msvc_2008_path)

                else:
                    self.include_dirs = [msvc_2008_path]

            elif (sys.version_info.major, sys.version_info.minor) < (3, 5):
                # Actually, it seems that appveyor doesn't have a stdint that
                # works, so even for 2010 we use our own (hacked) version
                # of stdint.
                # This should be pretty safe in whatever case.
                msvc_2010_path = (
                    os.path.join(os.getcwd(), 'include', 'msvc_2010'))

                if self.include_dirs is not None:
                    self.include_dirs.append(msvc_2010_path)

                else:
                    self.include_dirs = [msvc_2010_path]

            # We need to prepend lib to all the library names
            _libraries = []
            for each_lib in self.libraries:
                _libraries.append('lib' + each_lib)

            self.libraries = _libraries


def setup_package():

    # Figure out whether to add ``*_requires = ['numpy']``.
    build_requires = []
    try:
        import numpy
    except ImportError:
        build_requires = ['numpy']

    if not READ_THE_DOCS:
        install_requires = ['matplotlib>=1.3.0', 'scipy>=1.10',
                            'bottleneck', 'obspy>=1.3.0', 'numpy>=1.22.2',
                            'h5py']
    else:
        install_requires = ['matplotlib>=1.3.0', 'obspy>=1.0.3', 'mock']
    install_requires.extend(build_requires)

    setup_args = {
        'name': 'EQcorrscan',
        'version': VERSION,
        'description':
        'EQcorrscan - matched-filter earthquake detection and analysis',
        'long_description': long_description,
        'url': 'https://github.com/eqcorrscan/EQcorrscan',
        'author': 'Calum Chamberlain',
        'author_email': 'calum.chamberlain@vuw.ac.nz',
        'license': 'LGPL-3.0-or-later',
        'license-files': ['LICENCE.txt'],
        'classifiers': [
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            'Programming Language :: Python :: 3.13',
        ],
        'keywords': 'earthquake correlation detection match-filter',
        'scripts': scriptfiles,
        'install_requires': install_requires,
        'setup_requires': ['pytest-runner'],
        'tests_require': ['pytest>=2.0.0', 'pytest-cov', 'pytest-pep8',
                          'pytest-xdist', 'pytest-rerunfailures',
                          'obspy>=1.1.0'],
        'cmdclass': {'build_ext': CustomBuildExt},
        'packages': [
            'eqcorrscan', 'eqcorrscan.utils', 'eqcorrscan.core',
            'eqcorrscan.core.match_filter',
            'eqcorrscan.core.match_filter.helpers', 'eqcorrscan.utils.lib',
            'eqcorrscan.tutorials', 'eqcorrscan.helpers', 'eqcorrscan.tests'],
    }

    if using_setuptools:
        setup_args['setup_requires'] = build_requires
        setup_args['install_requires'] = install_requires

    no_mkl = False
    if "--no-mkl" in sys.argv:
        no_mkl = True
        sys.argv.remove("--no-mkl")
    if len(sys.argv) >= 2 and (
        '--help' in sys.argv[1:] or
        sys.argv[1] in ('--help-commands', 'egg_info', '--version',
                        'clean')):
        # For these actions, NumPy is not required.
        pass
    else:
        setup_args['ext_modules'] = get_extensions(no_mkl=no_mkl)
        setup_args['package_data'] = get_package_data()
        setup_args['package_dir'] = get_package_dir()
    if os.path.isdir("build"):
        shutil.rmtree("build")
    setup(**setup_args)


if __name__ == '__main__':
    setup_package()
