EQcorrscan API
==============

EQcorrscan contains two main modules, core and utils. Core contains most of the
most useful functions (including matched-filtering), and utils contains a range
of possibly helpful functions that may be useful to you and/or are used by the
core routines.

Core
----

Core routines of the EQcorrscan project.

.. toctree::
   :maxdepth: 1

   submodules/core.lag_calc
   submodules/core.match_filter
   submodules/core.subspace
   submodules/core.template_gen

Utils
-----

Various utility functions to help the core routines, and/or for use to analyse the results
of the core routines.

.. toctree::
   :maxdepth: 1

   submodules/utils.archive_read
   submodules/utils.catalog_to_dd
   submodules/utils.catalog_utils
   submodules/utils.clustering
   submodules/utils.correlate
   submodules/utils.despike
   submodules/utils.findpeaks
   submodules/utils.mag_calc
   submodules/utils.picker
   submodules/utils.plotting
   submodules/utils.pre_processing
   submodules/utils.sac_util
   submodules/utils.stacking
   submodules/utils.synth_seis
   submodules/utils.trigger