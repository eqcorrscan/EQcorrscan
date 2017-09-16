try:
    # use setuptools if we can
    from setuptools import setup, Command, Extension
    from setuptools.command.build_ext import build_ext
    using_setuptools = True
except ImportError:
    from distutils.core import setup, Command, Extension
    from distutils.command.build_ext import build_ext
    using_setuptools = False

from distutils.ccompiler import get_default_compiler

import os
import sys
import shutil
import glob
import eqcorrscan

VERSION = eqcorrscan.__version__

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
        package_data['eqcorrscan.lib'] = [
            'libfftw3-3.dll', 'libfftw3f-3.dll', 'libfftw3l-3.dll']

    return package_data


def get_package_dir():
    from pkg_resources import get_build_platform

    package_dir = {}
    if get_build_platform() in ('win32', 'win-amd64'):
        package_dir['eqcorrscan.lib'] = os.path.join('eqcorrscan', 'lib')

    return package_dir


def get_include_dirs():
    import numpy
    from pkg_resources import get_build_platform

    include_dirs = [os.path.join(os.getcwd(), 'include'),
                    os.path.join(os.getcwd(), 'eqcorrscan', 'lib'),
                    numpy.get_include(),
                    os.path.join(sys.prefix, 'include')]

    if get_build_platform().startswith('freebsd'):
        include_dirs.append('/usr/local/include')

    return include_dirs


def get_library_dirs():
    from pkg_resources import get_build_platform

    library_dirs = []
    if get_build_platform() in ('win32', 'win-amd64'):
        library_dirs.append(os.path.join(os.getcwd(), 'eqcorrscan', 'lib'))
        library_dirs.append(os.path.join(sys.prefix, 'bin'))

    library_dirs.append(os.path.join(sys.prefix, 'lib'))
    if get_build_platform().startswith('freebsd'):
        library_dirs.append('/usr/local/lib')

    return library_dirs


def get_libraries():
    from pkg_resources import get_build_platform

    if get_build_platform() in ('win32', 'win-amd64'):
        libraries = ['libfftw3-3', 'libfftw3f-3']
    else:
        libraries = ['fftw3', 'fftw3_threads', 'fftw3f', 'fftw3f_threads']

    return libraries


def export_symbols(*path):
    """
    Required for windows systems - functions defined in libutils.def.
    """
    lines = open(os.path.join(*path), 'r').readlines()[2:]
    return [s.strip() for s in lines if s.strip() != '']


def get_extensions():
    from distutils.extension import Extension
    from pkg_resources import get_build_platform

    if READ_THE_DOCS:
        return []
    # will use static linking if STATIC_FFTW_DIR defined
    static_fftw_path = os.environ.get('STATIC_FFTW_DIR', None)
    link_static_fftw = static_fftw_path is not None

    common_extension_args = {
        'include_dirs': get_include_dirs(),
        'library_dirs': get_library_dirs()}

    sources = [os.path.join('eqcorrscan', 'lib', 'multi_corr.c'),
               os.path.join('eqcorrscan', 'lib', 'time_corr.c')]
    exp_symbols = export_symbols("eqcorrscan/lib/libutils.def")

    if get_build_platform() not in ('win32', 'win-amd64'):
        extra_link_args = ['-lm', '-lgomp']
        extra_compile_args = ['-fopenmp', '-ftree-vectorize', '-msse2']
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
        common_extension_args['libraries'] = libraries
        common_extension_args['extra_compile_args'] = extra_compile_args
        common_extension_args['export_symbols'] = exp_symbols
    ext_modules = [
        Extension('eqcorrscan.lib.libutils', sources=sources,
                  **common_extension_args)]
    return ext_modules


class CustomBuildExt(build_ext):
    def finalize_options(self):

        build_ext.finalize_options(self)

        if self.compiler is None:
            compiler = get_default_compiler()
        else:
            compiler = self.compiler

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
        build_requires = ['numpy>=1.6, <2.0']

    if not READ_THE_DOCS:
        install_requires = ['matplotlib>=1.3.0', 'scipy>=0.18', 'LatLon',
                            'bottleneck', 'obspy>=1.0.3', 'numpy>=1.12',
                            'h5py']
    else:
        install_requires = ['matplotlib>=1.3.0', 'LatLon', 'obspy>=1.0.3',
                            'mock']
    install_requires.extend(build_requires)

    setup_args = {
        'name': 'EQcorrscan',
        'version': VERSION,
        'description': 'EQcorrscan - matched-filter earthquake detection and analysis',
        'long_description': long_description,
        'url': 'https://github.com/calum-chamberlain/EQcorrscan',
        'author': 'Calum Chamberlain',
        'author_email': 'calum.chamberlain@vuw.ac.nz',
        'license': 'LGPL',
        'classifiers': [
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering',
            'License :: OSI Approved :: GNU Library or Lesser General Public '
            'License (LGPL)',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],
        'keywords': 'earthquake correlation detection match-filter',
        'scripts': scriptfiles,
        'install_requires': install_requires,
        'setup_requires': ['pytest-runner'],
        'tests_require': ['pytest>=2.0.0', 'pytest-cov', 'pytest-pep8',
                          'pytest-xdist', 'pytest-rerunfailures'],
        'cmdclass': {'build_ext': CustomBuildExt}
    }

    if using_setuptools:
        setup_args['setup_requires'] = build_requires
        setup_args['install_requires'] = install_requires

    if len(sys.argv) >= 2 and (
        '--help' in sys.argv[1:] or
        sys.argv[1] in ('--help-commands', 'egg_info', '--version',
                        'clean')):
        # For these actions, NumPy is not required.
        pass
    else:
        setup_args['packages'] = ['eqcorrscan', 'eqcorrscan.utils',
                                  'eqcorrscan.core', 'eqcorrscan.lib',
                                  'eqcorrscan.tutorials']
        setup_args['ext_modules'] = get_extensions()
        setup_args['package_data'] = get_package_data()
        setup_args['package_dir'] = get_package_dir()
    if os.path.isdir("build"):
        shutil.rmtree("build")
    setup(**setup_args)

if __name__ == '__main__':
    setup_package()
