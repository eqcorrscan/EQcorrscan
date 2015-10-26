"""A setuptools based setup module for EQcorrscan package.

    EQcorrscan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EQcorrscan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EQcorrscan.  If not, see <http://www.gnu.org/licenses/>.
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import eqcorrscan
# To use a consistent encoding
from codecs import open
from os import path
import warnings
import glob

try:
    import cv2
except:
    warnings.warn('##### No cv2 module, openCV, you need to install this yourself')

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Get a list of all the scripts not to be installed
scriptfiles=glob.glob('eqcorrscan/scripts/*.py')
scriptfiles+=glob.glob('eqcorrscan/*.sl')
scriptfiles+=glob.glob('eqcorrscan/tutorial.py')
scriptfiles+=glob.glob('eqcorrscan/WHATVsearch.py')
scriptfiles+=glob.glob('eqcorrscan/LFE_brightness_search.py')
scriptfiles+=glob.glob('eqcorrscan/synth_test.py')


setup(
    name='EQcorrscan',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=eqcorrscan.__version__,

    description='EQcorrscan - correlation earthquake detection',
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
        'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='earthquake correlation detection match-filter',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['docs', 'tests', 'test_data',\
                                    'grid', 'detections', 'templates',\
                                    'stack_templates', 'par']),

    scripts = scriptfiles,

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['obspy', 'numpy', 'matplotlib', 'joblib',\
     'scipy>=0.14'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #   'tutorial_data': ['test_data'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    data_files=[('tutorial_data', ['eqcorrscan/test_data/tutorial_data.tgz'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    # entry_points={
    #     'console_scripts': [
    #         'sample=sample:main',
    #     ],
    # },
)
