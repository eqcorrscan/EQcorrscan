Introduction to the EQcorrscan package
======================================

This document is designed to give you an overview of the capabilities and
implementation of the EQcorrscan Python module.

Why EQcorrscan?
---------------
EQcorrscan is designed to compute detections of earthquakes, or any seismic signal
(explosions work *really* well) by comparing templates with continuous data.
The main benefit of EQcorrscan's matched-filter routine is the level of parallel
processing that can be achieved.  By exploiting the fact that each template
does not rely on any other template, detections from a single template through
a day of seismic data can be computed in parallel.  By computing these in parallel
rather than a single template through multiple days we reduce IO load.  At a low
level, each time-step is computed in parallel by using the openCV matchTemplate
function.  The net result is that these functions are *very* scalable, we have
obtained a speed-up from 2 months to 10 hours by migrating from a small cluster
to a large one (for a 6.5 year long continuous dataset and 800 templates).

The authors of EQcorrscan foresee this project as an open repository for the
development of software for the detection and analysis of repeating and
near-repeating earthquakes.  This repository will continue to grow and develop
and any and all help/criticism will be appreciated.

There are a lot of things that could be added to this project - if you want to
get involved the best place to start, and the most valuable thing for your
understanding, and for the health of this package would be to contribute tests and
documentation.

Installation
------------

In general we recommend users to install EQcorrscan in a virtual environment,
for this the |virtualenvwrapper| package is handy.

Within a virtual environment, a fresh install should be as simple as:

**pip install eqcorrscan**

Most codes should work without any effort on your part.  However you will need to
install the openCV-python package yourself.  We recommend installing openCV version
3, and we recommend installing it from source - it is available via anaconda, but
it will run faster if you compile it yourself, and it will give more consistent
results.  See |pyimagesearch| for details for install on all operating systems
(including raspberry pi, which EQcorrscan runs on too :) ).

On Linux with Python 2.7:

**apt-get install python-opencv**

On OSX with Python 2.7:

**port install py27-numpy**
**port install opencv +python27**
or
**brew install opencv**

You can also install from source; for Python 3 this is a must as you will have
to install openCV 3.  |pyimagesearch| has lots of lovely tutorials like this
|cv3_ubuntu|.

.. |virtualenvwrapper| raw:: html

    <a href="https://virtualenvwrapper.readthedocs.io/en/latest/" target="blank">virtualenvwrapper</a>

.. |pyimagesearch| raw:: html

   <a href="http://www.pyimagesearch.com/" target="_blank">pyimagesearch</a>

.. |cv3_ubuntu| raw:: html

   <a href="http://www.pyimagesearch.com/2015/07/20/install-opencv-3-0-and-python-3-4-on-ubuntu/" target="_blank">install cv3 on ubuntu</a>

On Windows you can follow nice instructions |windows_opencv|.

.. |windows_opencv| raw:: html

   <a href="http://docs.opencv.org/3.1.0/d5/de5/tutorial_py_setup_in_windows.html#gsc.tab=0" target="_blank">here</a>

Note you may have issues with these installs if you don't have numpy installed: but if
you don't have numpy installed then you have bigger issues...

If you plan to run the *bright_lights* or generating a synthetic grid of
templates you will need to have grid csv files, which the authors have
previously used NonLinLoc to generate.  This is not provided here and should
be sourced from |NLLoc_link|. This will provide
the Grid2Time routine which is required to set-up a lag-time grid for your
velocity model.  You should read the NonLinLoc documentation for more
information regarding how this process works and the input files you are
required to give.

.. |NLLoc_link| raw:: html

  <a href="http://alomax.free.fr/nlloc/" target="_blank">NonLinLoc</a>

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

Tests will take about half an hour to run (as of v.0.1.4) and will provide
a coverage report at the end and notify you of any failures.