Core
====

Core programs for the EQcorrscan project.  To be added: Lag calculation by
cross-correlation to find optimum phase arrival picks for near-repeating
earthquakes.

bright_lights contains a series of functions to detect events using the
brightness-based beamforming method of Frank et. al (2014).  This has been
tested significantly, but has failed to detect events unambiguously in the
central Southern Alps.  As such development of these functions has ceased.

template_gen contains routines for cutting waveforms around picks for use as
templates in match_filter.  Included in this are wrappers to directly read in
Seisan formattaed pick files and waveforms associated with the picks, and
generate templates from these.  There are also wrappers for quakeML events
and catalogs, and seishub databases.

match_filter contains the core routines for earthquake detection by
cross-correlation.  This is optimized for large-scale, multi-paralleled
detection, with large numbers of templates.  Because we are unsure of your
architecture we have not written functions for the top level of possible
parallel computing, which would be to compute detections for multiple days
in parallel in a High-Performance Computing, cluster environment.  If you
want to know more about doing this please contact the authors.  We use
a cluster running SLURM for job scheduling and handle multiple days using
the batch job submission capability which distributes daily detections across
multiple nodes.  This allows us to detect earthquakes through > 6 years of
multi-channel data using > 600 templates in less than 36 hours.

.. toctree::
   :maxdepth: 4

   submodules/core.bright_lights
   submodules/core.template_gen
   submodules/core.match_filter
