EQcorrscan tutorials
====================
Welcome to EQcorrscan - this package is designed to compute earthquake detections
using a paralleled matched-filter network cross-correlation routine, and analyse the
results.

Before continuing with this tutorial please check that you have installed all
the pre-requisite modules, as not all will be installed by the setup.py file.
The list of these is in the :doc:`Introduction <intro>` section of this documentation.

As you will see, this package is divided into two main sub-modules, the
:doc:`core </core>` and :doc:`utils </utils>` sub-modules.
The core sub-module contains the main, high-level functions:

:bright_lights:
        A brightness based template detection routine;
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

Some other high-level functions are included in the :doc:`utils </utils>` sub-module
and are documented here with tutorials:

:mag_calc:
        Simple local magnitude calculation and high-precision relative moment
        calculation using singular-value decomposition.
:clustering:
        Routines for clustering earthquakes based on a range of metircs using
        agglomorative clustering methods.

The :doc:`utils </utils>` sub-module contains useful, but small functions.
These functions are rarely cpu intensive, but perform vital operations, such
as reading |Seisan| s-files (:doc:`sfile_util </submodules/utils.sfile_util>`),
finding peaks in noisy data (:doc:`findpeaks </submodules/utils.findpeaks>`),
converting a seisan database to hypoDD formatted files and computing cross-correlations between
detections for |hypoDD| (a double difference relocation software)
(:doc:`catalog_to_dd </submodules/utils.catalog_to_dd>`), calculating
magnitudes (:doc:`mag_calc </submodules/utils.mag_calc>`),
clustering detections (:doc:`clustering </submodules/utils.clustering>`),
stacking detections (:doc:`stacking </submodules/utils.stacking>`),
making pretty plots (:doc:`plotting </submodules/utils.plotting>`),
and processing seismic data in the same way repeatedly using |Obspy|'s
functionality (:doc:`pre_processing </submodules/utils.pre_processing>`).

What follows is an expanding set of tutorials that should take you
through some of the key functionality of the EQcorrscan package.

.. toctree::
  :numbered:
  :titlesonly:

  tutorials/template-creation.rst
  tutorials/matched-filter.rst
  tutorials/subspace.rst
  tutorials/lag-calc.rst
  tutorials/mag-calc.rst
  tutorials/clustering.rst

.. |Harris2006| raw:: html

    <a href="https://e-reports-ext.llnl.gov/pdf/335299.pdf" target="_blank">Harris (2006)</a>

.. |HypoDD| raw:: html

    <a href="http://www.ldeo.columbia.edu/~felixw/hypoDD.html" target="_blank">HypoDD</a>

.. |Seisan| raw:: html

  <a href="http://seisan.info/" target="_blank">Seisan</a>

.. |Obspy| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">Obspy</a>