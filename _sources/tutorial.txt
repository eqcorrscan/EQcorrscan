EQcorrscan tutorial
===================
Welcome to EQcorrscan - this package is designed to compute earthquake detections
using a paralleled match-filter network cross-correlation routine.  The inner
loop of this package is the cross-correaltion of templates of seismic data
with daylong seismic data.  This inner function is the openCV.match_template
function - this appears to be a well optimized cross-correlation function written
in c++.  Cross-correlations are computed in the frequency domain for large
datasets, for which a day of seismic data usually qualifies.

Before continuing with this tutorial please check that you have installed all
the pre-requisite modules, as this won't be done for you.  The list of these is
in the Introduction section of this documentation.

As you will see, this package is divided into three main sub-modules, the
Core, Utils and Par sub-modules.  The Core sub-module contains the main, high-level
functions:

:bright_lights: 
        A brightness based template detection routine;
:template_gen:
        A series of routines to generate templates for match-filter detection
        from continuous or cut data, with pick-times defined either manually, or from a
        *Seian* s-file;
:match_filter:
        The main match-filter routines, this is split into several
        smaller functions to allow python based parallelisation;
:lag_calc:
        Routines for calculating optimal lag-times for events detected
        by the match-filter routine, these lags can then be used to define new picks
        for high accuracy reloactions.

The Par sub-module contains parameter files which are provided to allow for simple
bulk processing of large datasets.  These *MUST* be edited by the user for their
dataset.

The Utils sub-module contains useful, but small functions.  These functions are
rarely cpu intensive, but perform vital operations, such as reading *Seisan* s-files,
finding peaks in noisy data, converting a seisan database to hypoDD formatted
files and computing cross-correlations between detections for hypoDD (a double
difference reloaction software), calculating magnitudes, clustering detections,
stacking detections, making pretty plots, and processing seismic data in the
same way repeatedly using *Obspy*'s functionality.


Match-filter detection
----------------------

In this section we will discuss generating a template from a *Seisan* s-file, and
using this template to scan for similar earthquakes within a day of data.  This single
template and single day of data does not truely exploit the parallel operations within
this package however, so you would be encouraged to think about where parallel operations
occur (*hint, the code can run one template per cpu*), and why there are --instance and
--splits flags (*hint, if you have heaps of memory and cpus you can do some brute
force day parallelisation!*).

The following script is included in the top-level directory alongside the full-scripts
used by the author to generate a 6.5 year long catalogue of low-frequency earthquakes
for the central Southern Alps of New Zealand.

This tutorial script highlights the ability of the match-filter method in detecting
earthquakes of near-repeating nature.  The dataset is a day of data taken from the
New Zealand national database, and the Southern Alp Microearthquake Borehole Array
(SAMBA) network (Boese et al. 2012).  This day was found to contain a swarm of
earthquakes, as published by Boese et al. (2014), the s-file provided is one of
these events.

The main processing flow is outlined in the figure below, note the main speedups
in this process are achioeved by running multiple templates at once, however this
increaces memory usage.  If memory is problem there are flags (mem_issue) in the
match_filter.py source that can be turned on - the codes will then write temporary
files, which is slower, but can allow for more data crunching at once, your trade-off,
your call.

.. image:: processing_flow.png
     :width: 600px
     :align: center
     :alt: processing_flow.png

References
----------
* CM Boese, J Townend, E Smith, T Stern (2012). `Microseismicity and stress in the vicinity of the Alpine Fault, central Southern Alps, New Zealand <http://onlinelibrary.wiley.com/doi/10.1029/2011JB008460/full>`_, *JGR*, doi:10.1029/2011JB008460
* CM Boese, KM Jacobs, EGC Smith, TA Stern, J Townend (2014). `Background and delayed-triggered swarms in the central Southern Alps, South Island, New Zealand <http://onlinelibrary.wiley.com/doi/10.1002/2013GC005171/full>`_, *G-cubed*, doi:10.1002/2013GC005171

.. literalinclude:: ../tutorial.py
