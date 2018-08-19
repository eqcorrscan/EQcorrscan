Introduction to the EQcorrscan package
======================================

This document is designed to give you an overview of the capabilities and
implementation of the EQcorrscan Python package.

Why EQcorrscan?
---------------
EQcorrscan is designed to compute detections of earthquakes, or any seismic signal
(explosions work *really* well) using more advanced routines than standard
amplitude-ratio methods.

The authors of EQcorrscan foresee this project as an open repository for the
development of software for the detection and analysis of repeating and
near-repeating earthquakes.  This repository will continue to grow and develop
and any and all help/criticism will be appreciated.

There are a lot of things that could be added to this project - if you want to
get involved a good place to start would be to contribute tests and
documentation.

Installation - Updated for versions > 0.2.7
-------------------------------------------

Recommended install method
~~~~~~~~~~~~~~~~~~~~~~~~~~

In general we recommend users to install EQcorrscan in a virtual environment,
|conda| will simplify your install greatly (install instructions for anaconda
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

Not-recommended but workable methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are other ways to install EQcorrscan and it's dependencies,
some details for some of those cases can be found :doc:`here</installation>`.
We do note that the other methods are more involved, and can be problematic. If
you do chose not to use conda then you should definitely test your install (see below).

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

We support Linux, OSX and Windows environments running Python 2.7 and 3.x.

We do **not** support Python 2.6.

We will stop support for Python 2.7 in a forthcoming release, for more information
see |#242|.

.. |#242| raw:: html

   <a href="https://github.com/eqcorrscan/EQcorrscan/issues/242" target="_blank">issue #242</a>


Functionality
-------------

Within :doc:`core </core>` you will find the core routines to generate templates,
(:doc:`template_gen </submodules/core.template_gen>`) search for likely templates
(:doc:`bright_lights </submodules/core.bright_lights>`),
compute cross-channel correlations from these templates
(:doc:`match_filter </submodules/core.match_filter>`), generate cross-correlation
corrected pick-times (:doc:`lag_calc </submodules/core.lag_calc>`),
and run subspace detection (:doc:`subspace </submodules/core.subspace>`).

The bright_lights and match_filter submodules have been designed with parallel
computing in mind, to the extent that the more cores and machines you have
running them (generally) the better.  These rely on the python multiprocessing
module, and some C extensions using openMP to handle parallelisation at
lower-levels.

.. _RunningTests:

Running tests
-------------

One of the main goals of EQcorrscan is to improve reliability and reproducibility
of earthquake detection.  To this end, EQcorrscan has a set of tests (you
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

From version 0.3.2 onwards EQcorrscan includes a test script that will be installed
onto your path when you install EQcorrscan.  This test-script will download the
test data and run the tests (you no longer have to clone the git repository). Just run
(from anywhere):

.. code-block:: bash

    test_eqcorrscan.py

Tests will take about half an hour to run (as of v.0.3.2) and will provide
a coverage report at the end and notify you of any failures.
