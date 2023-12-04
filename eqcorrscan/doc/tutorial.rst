EQcorrscan tutorials
====================
Welcome to EQcorrscan - this package is designed to compute earthquake detections
using a paralleled matched-filter network cross-correlation routine, and analyse the
results.

EQcorrscan is divided into two main sub-modules, the
**core** and **utils** sub-modules.
The core sub-module contains the main, high-level functions:

:template_gen:
        A series of routines to generate templates for match-filter detection
        from continuous or cut data, with pick-times either defined manually,
        or defined in event files;
:match_filter:
        The main matched-filter routines, this is split into several
        smaller functions to allow python-based parallel-processing;
:subspace:
        Subspace detection routine based on |Harris2006|.
:lag_calc:
        Routines for calculating optimal lag-times for events detected
        by the match-filter routine, these lags can then be used to define new picks
        for high accuracy re-locations.

Some other high-level functions are included in the :doc:`utils </api>` sub-module
and are documented here with tutorials:

:mag_calc:
        Simple local magnitude calculation and high-precision relative moment
        calculation using singular-value decomposition.
:clustering:
        Routines for clustering earthquakes based on a range of metircs using
        agglomorative clustering methods.

The **utils** sub-module contains useful, but small functions.
These functions are rarely cpu intensive, but perform vital operations, such
as finding peaks in noisy data (:doc:`findpeaks </submodules/utils.findpeaks>`),
converting a database to hypoDD formatted files and computing cross-correlations between
detections for |hypoDD| (a double difference relocation software)
(:doc:`catalog_to_dd </submodules/utils.catalog_to_dd>`), calculating
magnitudes (:doc:`mag_calc </submodules/utils.mag_calc>`),
clustering detections (:doc:`clustering </submodules/utils.clustering>`),
stacking detections (:doc:`stacking </submodules/utils.stacking>`),
making plots (:doc:`plotting </submodules/utils.plotting>`),
and processing seismic data in the same way repeatedly using |Obspy|'s
functionality (:doc:`pre_processing </submodules/utils.pre_processing>`).

As of EQcorrscan v.0.4.0, output is controlled by the |logging| module.
This allows the user to decide now much output is needed (via the `level`
option), and where output should be sent (either to a file, stdout, or both).
By default, if no logging parameters are set before calling EQcorrscan functions
only messages of WARNING or above are printed to stdout. To make things a little
prettier you can begin your work using something like:

.. code-block:: python

  import logging

  logging.basicConfig(
      level=logging.INFO,
      format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

This will output messages at level INFO and above (quite verbose), in messages
like:

.. code-block:: python

  2018-06-25 22:48:17,113	eqcorrscan.utils.pre_processing	INFO	Working on: NZ.CPWZ.10.EHZ

For more options, see the |logginghowto| guide.

The following is an expanding set of tutorials that should take you
through some of the key functionality of the EQcorrscan package.  The tutorials
are a (slow) work in progress to convert to jupyter notebook.  If you have anything
you want to see in the tutorials please let us know on github.

.. toctree::
  :numbered:
  :titlesonly:

  tutorials/quick_start.ipynb
  tutorials/local_files.ipynb
  tutorials/matched-filter.rst
  tutorials/template-creation.rst
  tutorials/subspace.rst
  tutorials/lag-calc.rst
  tutorials/mag-calc.rst
  tutorials/clustering.rst

.. |Harris2006| raw:: html

    <a href="https://e-reports-ext.llnl.gov/pdf/335299.pdf" target="_blank">Harris (2006)</a>

.. |HypoDD| raw:: html

    <a href="http://www.ldeo.columbia.edu/~felixw/hypoDD.html" target="_blank">HypoDD</a>

.. |Obspy| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">Obspy</a>

.. |logging| raw:: html

  <a href="https://docs.python.org/3/library/logging.html" target="_blank">Logging</a>

.. |logginghowto| raw:: html

  <a href="https://docs.python.org/3/howto/logging.html#logging-basic-tutorial" target="_blank">Logging HOWTO</a>
