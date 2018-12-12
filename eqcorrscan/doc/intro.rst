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

To test your installed version of EQcorrscan we provide a |test-script|.  For
version<=0.3.2 you should download the script and run it. In later versions this
script is included in the package.

.. |test-script| raw:: html

    <a href="https://gist.github.com/calum-chamberlain/0887455551862a363a43887f0195ec06" target="_blank">test-script</a>

This test-script will download the test data and run the tests (you no longer
have to clone the git repository). Just run (from anywhere):

.. code-block:: bash

    test_eqcorrscan.py

Tests will take about half an hour to run (as of version 0.3.2) and will provide
a coverage report at the end and notify you of any failures.
