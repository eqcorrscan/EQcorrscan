Introduction to the EQcorrscan package
======================================

This document is designed to give you an overview of the capabilities and
implementation of the EQcorrscan python module.

Why EQcorrscan?
---------------
EQcorrscan is designed to compute matched-filter detections of earthquakes,
or any seismic signal (explosions work *really* well) by comparing templates
with continuous data.  The main benefit of EQcorrscan is the level of
parallel processing that can be achieved.  By exploiting the fact that each template
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

We have a long way to go with this project - if you want to get involved the
best place to start, and the most valuable thing for your understanding, and
for the health of this repository would be to contribute tests and
documentation.  Ideally we would like to have one test for every function!

Installation
------------
A fresh install should be as simple as:

**pip install eqcorrscan**

Most codes should work without any effort on your part.  However you may need to
install the openCV-python package yourself.

This install has only been tested on Linux and OSX machines.  You
should be prepared for small differences in the results of your correlations
relating to foating-point truncation differences between 32 and 64-Bit
machines.

If you plan to run the bright_lights or generating a synthetic grid of
templates you will need to have grid csv files, which the authors have
previously used NonLinLoc to generate.  This is not provided here and should
be sourced from |NLLoc_link|. This will provide
the Grid2Time routine which is required to set-up a lag-time grid for your
velocity model.  You should read the NonLinLoc documentation for more
information regarding how this process works and the input files you are
required to give.

.. |NLLoc_link| raw:: html

  <a href="http://alomax.free.fr/nlloc/" target="_blank">NonLinLoc</a>

Functions
---------

This package is divided into sub-directories of *core* and *utils*.  The
*utils* directory contains simple functions for integration with |seisan_link|,
these are in the *Sfile_util.py*
module and functions therein which are essentially barebones and do not have the
full functionality that seisan can handle.  *utils* also contains a simple
peak-finding algorithm *find_peaks.py* which looks for peaks within noisy data
above a certain threshold and within windows.  Many other functions have been
added to this module to handle the analysis of repeating and near-repeating
earthquakes, including stacking routines, clustering algorithms, magnitude
calculation both by amplitude picking and by singular value decomposition.  I
recommend you take a look in here to see if any of it is useful.  There are also
some plotting routines that make handling large datasets a little simpler.  Most
recently I have added a simple synthetic seismogram generator, which is currently
my main project focus.

.. |seisan_link| raw:: html

  <a href="http://seisan.info/" target="_blank">Seisan</a>

Since earlier versions the *core* modules have moved away from using parameter
files, and instead rely on explicit argument calls.  The parameter files are
still included by not documented here (see inside the par files), and remain
useful when generating batch scripts (see the scripts in the github repo).

Within *core* you will find the core routines to generate templates,
*(template_gen)* search for likely templates *(bright_lights)* and
compute cross-channel correlations from these templates *(match_filter)*.  The
bright_lights and match_filter submodules have been designed with parallel
computing in mind, to the extent that the more cores and machines you have
running them the better.  These rely on the python multiprocesisng module to
handle parallelisation at lower-levels.  You can also do some 'brute-force'
parallelisation on a day level when computing detections over multiple days.
I tend to run one day per node of a cluster computer, with each day running
templates in parallel.
