"""A setuptools based setup module for EQcorrscan package.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import sys
import os
import eqcorrscan
# To use a consistent encoding
from codecs import open
from os import path
import warnings
import glob
from distutils.extension import Extension
import numpy as np
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    msg = ' '.join([])
    msg = ' '.join(["warning: pypandoc module not found,",
                    " could not convert Markdown to RST"])
    print(msg)
    read_md = lambda f: open(f, 'r').read()

READ_THE_DOCS = os.environ.get('READTHEDOCS', None) == 'True'
if not READ_THE_DOCS:
    from Cython.Distutils import build_ext
    ext = [Extension("eqcorrscan.core.subspace_statistic",
                     ["eqcorrscan/core/subspace_statistic.pyx"],
                     include_dirs=[np.get_include()])]
    cmd_class = {'build_ext': build_ext}
else:
    ext = []
    cmd_class = {}

try:
    import cv2  # NOQA
except ImportError:
    print(sys.path)
    msg = '##### No cv2 module, openCV, you need to install this yourself'
    warnings.warn(msg)

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
long_description = "EQcorrscan: matched-filter earthquake detection and " +\
    "analysis in Python.  Open-source routines for: systematic template " +\
    "creation, multi-parallel matched-filter detection, clustering of " +\
    "events, integration with SEISAN, SAC, QuakeML and NonLinLoc, " +\
    "magnitude calculation by singular value decomposition, and more!"

# Get a list of all the scripts not to be installed
scriptfiles = glob.glob('eqcorrscan/tutorials/*.py')
scriptfiles += glob.glob('eqcorrscan/scripts/*.py')

if sys.version_info.major == 2:
    if not READ_THE_DOCS:
        install_requires = ['numpy>=1.8.0', 'obspy>=1.0.0',
                            'matplotlib>=1.3.0', 'joblib>=0.8.4',
                            'scipy>=0.14', 'multiprocessing',
                            'LatLon', 'h5py', 'cython']
    else:
        install_requires = ['numpy>=1.8.0', 'obspy>=1.0.0',
                            'matplotlib>=1.3.0', 'joblib>=0.8.4',
                            'multiprocessing',
                            'LatLon']
else:
    if not READ_THE_DOCS:
        install_requires = ['numpy>=1.8.0', 'obspy>=0.10.2',
                            'matplotlib>=1.3.0', 'joblib>=0.8.4',
                            'scipy>=0.14', 'LatLon', 'h5py', 'cython']
    else:
        install_requires = ['numpy>=1.8.0', 'obspy>=0.10.2',
                            'matplotlib>=1.3.0', 'joblib>=0.8.4',
                            'LatLon']
# install_requires.append('ConfigParser')
setup(
    name='EQcorrscan',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=eqcorrscan.__version__,

    description='EQcorrscan - matched-filter earthquake detection and analysis',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/calum-chamberlain/EQcorrscan',

    # Author details
    author='Calum Chamberlain',
    author_email='goride42@gmail.com',

    # Choose your license
    license='LGPL',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU Library or Lesser General Public ' +
        'License (LGPL)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.4',
    ],

    # What does your project relate to?
    keywords='earthquake correlation detection match-filter',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['docs', 'tests', 'test_data',
                                    'grid', 'detections', 'templates',
                                    'stack_templates', 'par']),

    scripts=scriptfiles,

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=install_requires,

    # Test requirements for using pytest
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-cov', 'pytest-pep8'],

    # Build our extension for subspace detection
    cmdclass=cmd_class,
    ext_modules=ext
)
