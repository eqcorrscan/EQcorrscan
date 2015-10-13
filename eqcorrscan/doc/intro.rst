Introduction to the EQcorrscan package
======================================

This document is designed to give you an overview of the capabilities and
implimentation of the EQcorrscan python module.

Installation
------------

Most codes should work without any effort on your part.  However you must
install the packages this package relies on yourself, this includes the follwing
packages:
        * matplotlib
        * numpy
        * scipy
        * obspy
        * joblib
        * openCV (2)

This install has only been tested on Linux machines and even then has some
issues when installing on 32-Bit versus 64-Bit machines.  In this instance you
should be prepared for small differences in the results of your correlations
relating to foating-point truncation differences between 32 and 64-Bit
machines.

If you plan to run the bright_lights.py routines you will need to have
NonLinLoc installed on your machine.  This is not provided here and should
be sourced from `NonLinLoc <http://alomax.free.fr/nlloc/>`_ This will provide
the Grid2Time routine which is required to set-up a lag-time grid for your
velocity model.  You should read the NonLinLoc documentation for more
information regarding how this process works and the input files you are
required to give.

Functions
---------

This package is divided into sub-directories of *core*, *par* and *utils*.  The
*utils* directory contains simple functions for integration with
`seisan <http://seisan.info/>`_, these is the *Sfile_util.py* 
module and functions therein which are essentially barebones and do not have the
full functionality that seisan can handle.  *utils* also contains a simple
peak-finding algorithm *find_peaks.py* which looks for peaks within noisy data
above a certain threshold and within windows.  

Within *par* you will find parameter files which you will need to edit for
each of the *core* scripts.  *core* scripts often call on multiple *par* files
so be sure to set them all up.  The *template_gen_par.py* script is used by all
*core* modules and must be set-up.  Within this you will define all your
template parameters.  Currently the templates must all be of the same length,
but this may change in a future release.

Within *core* you will find the core routines to generate templates,
*(template_gen.py)* search for likely templates *(bright_lights.py)* and
compute cross-channel correlations from these templates *(match_filter.py)*.
