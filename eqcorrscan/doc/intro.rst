Introduction to the EQcorrscan package
======================================

Supported environments
----------------------

We support Linux, OSX and Windows environments running Python 3.x.

Functionality
-------------

Within :doc:`core </core>` you will find the core routines to generate templates
(:doc:`template_gen </submodules/core.template_gen>`),
compute cross-channel correlations from these templates
(:doc:`match_filter </submodules/core.match_filter>`), generate cross-correlation
corrected pick-times (:doc:`lag_calc </submodules/core.lag_calc>`),
and run subspace detection (:doc:`subspace </submodules/core.subspace>`).

.. _RunningTests:

Running tests
-------------

You can run tests yourself locally to ensure that everything runs as you would expect
in your environment.  Although we try to ensure that these tests run smoothly on all
supported environments (using the ci bots), if you do find any issues, please let us
know on the |github| page.

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
