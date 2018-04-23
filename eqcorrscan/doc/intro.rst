Introduction to the EQcorrscan package
======================================

This document is designed to give you an overview of the capabilities and
implementation of the EQcorrscan Python package.

Why EQcorrscan?
---------------
EQcorrscan is designed to compute detections of earthquakes, or any seismic signal
(explosions work *really* well) using more advanced routines than standard
amplitude-ratio methods.

This package was originally based around a matched-filter detection routine
which works by comparing templates with continuous data.
The main benefit of EQcorrscan's matched-filter routine is the level of parallel
processing that can be achieved.  EQcorrscan will run on anything from a 1GB RAM
single-board computer to a multi-hundred-GB RAM, thousand CPU high-performance
computer.  Because the internals of EQcorrscan's matched-filter routine scale
reasonably well, the developers have observed speed-ups of 150x (from 2 months
to 10 hours) by migrating from a small cluster
to a large one (for a 6.5 year long continuous dataset and 800 templates).

The authors of EQcorrscan foresee this project as an open repository for the
development of software for the detection and analysis of repeating and
near-repeating earthquakes.  This repository will continue to grow and develop
and any and all help/criticism will be appreciated.

There are a lot of things that could be added to this project - if you want to
get involved the best place to start, and the most valuable thing for your
understanding, and for the health of this package would be to contribute tests and
documentation.

Installation - Updated for version 0.2.7
----------------------------------------

Recommended install method
~~~~~~~~~~~~~~~~~~~~~~~~~~

In general we recommend users to install EQcorrscan in a virtual environment,
|conda| will simplify your install greatly (install instuction for anaconda 
or miniconda are here: |conda-install|) - we recommend creating a conda
environment with the following:

.. code-block:: bash

    conda config --add channels conda-forge
    conda create -n eqcorrscan colorama numpy scipy matplotlib obspy bottleneck pyproj
    # If using bash run:
    source activate eqcorrscan
    # If not you should get away with:
    activate eqcorrscan
    
To then install EQcorrscan you can simply run:

.. code-block:: bash

    conda install eqcorrscan

Not-recomended but workable methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a myriad of other ways to install EQcorrscan and it's dependancies,
some details for some of those cases can be found :doc:`here</installation>`.
We do note that the other methods are more involved, and can be problematic. If
you do chose not to use conda then you should definitely test your install.

.. |conda| raw:: html

    <a href="https://conda.io/docs/" target="_blank">conda</a>

.. |conda-install| raw:: html

    <a href="https://conda.io/docs/user-guide/install/index.html#installing-conda-on-a-system-that-has-other-python-installations-or-packages" target="_blank">conda installation</a>

.. |fftw-install| raw:: html

    <a href="http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix" target="_blank">fftw installation</a>

.. |fftw-windows| raw:: html

    <a href="http://www.fftw.org/install/windows.html" target="_blank">fftw-windows install</a>

.. |pyasdf| raw:: html

    <a href="http://seismicdata.github.io/pyasdf/index.html" target="_blank">pyASDF</a>

.. |virtualenvwrapper| raw:: html

    <a href="https://virtualenvwrapper.readthedocs.io/en/latest/" target="blank">virtualenvwrapper</a>

.. |pyimagesearch| raw:: html

   <a href="http://www.pyimagesearch.com/" target="_blank">pyimagesearch</a>

.. |cv3_ubuntu| raw:: html

   <a href="http://www.pyimagesearch.com/2015/07/20/install-opencv-3-0-and-python-3-4-on-ubuntu/" target="_blank">install cv3 on ubuntu</a>

Supported environments
----------------------

We support Linux, OSX and Windows environments running Python 2.7, 3.4 and 3.5.
We don't run our tests on other versions of Python so you might have some issues
with other Python 3.x series, if you do, let us know.

We do **not** support Python 2.6.


Functions
---------

This package is divided into sub-directories of :doc:`core </core>` and :doc:`utils </utils>`.  The
:doc:`utils </utils>` directory contains simple functions for integration with |seisan_link|,
these are in the :doc:`sfile_util </submodules/utils.sfile_util>`
module and functions therein which are essentially barebones and do not have the
full functionality that seisan can handle.  :doc:`utils </utils>` also contains a simple
peak-finding algorithm :doc:`findpeaks </submodules/utils.findpeaks>` which looks for peaks within noisy data
above a certain threshold and within windows.

Many other functions have been
added to this module to handle the analysis of repeating and near-repeating
earthquakes, including stacking routines, clustering algorithms, magnitude
calculation both by amplitude picking and by singular value decomposition.  I
recommend you take a look in here to see if any of it is useful.  There are also
some plotting routines that make handling large datasets a little simpler.  Most
recently I have added a simple synthetic seismogram generator, which is currently
my main project focus.

.. |seisan_link| raw:: html

  <a href="http://seisan.info/" target="_blank">Seisan</a>

Since earlier versions the :doc:`core </core>` modules have moved away from using parameter
files, and instead rely on explicit argument calls.  The parameter files are
still included by not documented here (see inside the par files), and remain
useful when generating batch scripts (see the scripts in the github repo).

Within :doc:`core </core>` you will find the core routines to generate templates,
(:doc:`template_gen </submodules/core.template_gen>`) search for likely templates
(:doc:`bright_lights </submodules/core.bright_lights>`) and
compute cross-channel correlations from these templates (:doc:`match_filter </submodules/core.match_filter>`).  The
bright_lights and match_filter submodules have been designed with parallel
computing in mind, to the extent that the more cores and machines you have
running them the better.  These rely on the python multiprocessing module to
handle parallelisation at lower-levels.  You can also do some 'brute-force'
parallelisation on a day level when computing detections over multiple days.
I tend to run one day per node of a cluster computer, with each day running
templates in parallel.

.. _RunningTests:

Running tests
-------------

One of the main goals of EQcorrscan is to improve reliability and reproducibility
of earthquake detection.  To this end, EQcorrscan has a moderate test-base (you
can check how much of our codebase if tested by looked at the badges in the
|github| repository).  You can also run these tests yourself locally to ensure
that everything runs as you would expect in your environment.  Although every
effort has been made to ensure these tests run smoothly on all supported environments
(using the ci bots), if you do find any issues, please let us know on the
|github| page.

.. |github| raw:: html

    <a href="https://github.com/eqcorrscan/EQcorrscan" target="_blank">github</a>

To run the tests you will need to have pytest installed along with a couple of
extras (pytest-pep8 and pytest-cov).  These can be installed by pip:

.. code-block:: bash

    pip install pytest pytest-pep8 pytest-cov

You will also need to have a clone of the github repository:

.. code-block:: bash

    git clone https://github.com/eqcorrscan/EQcorrscan.git

You can then run the tests from within the repository directory:

.. code-block:: bash

    python setup.py test
    
Note that if you needed to prepend CC=gcc (if the default compiler is `clang`)
then you will also need to here for the setup.py command.

If this fails with an error like:

.. code-block:: bash

    error: invalid command 'pytest'

If you already have an install of EQcorrscan then you can instead run the
tests using the following:

.. code-block:: bash

    pip install requirements.txt
    py.test

Otherwise, if you have the pytest error and no previous install of EQcorrscan
you can install from source and run the tests using:

.. code-block:: bash

    python setup.py install
    py.test

Tests will take about half an hour to run (as of v.0.3.0) and will provide
a coverage report at the end and notify you of any failures.
