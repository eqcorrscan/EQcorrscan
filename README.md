# EQcorrscan
## A python package for the detection and analysis of repeating and near-repeating earthquakes.

## Citation:
We have a manuscript on the development of EQcorrscan, if you make use of EQcorrscan please cite the following paper:

Chamberlain, C. J., Hopp, C. J., Boese, C. M., Warren-Smith, E., Chambers, D., Chu, S. X., Michailos, K., Townend, J., [EQcorrscan: Repeating and near-repeating earthquake detection and analysis in Python.](https://pubs.geoscienceworld.org/ssa/srl/article/89/1/173/524875/eqcorrscan-repeating-and-near-repeating-earthquake) Seismological Research Letters *2017*

If you want to you should also cite the version number:
[![DOI](https://zenodo.org/badge/35918157.svg)](https://zenodo.org/badge/latestdoi/35918157)

# Installation

The easiest way to install EQcorrscan is through anaconda:
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/eqcorrscan/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)

Instructions for installing EQcorrscan and the required dependency, fftw
are linked from the
[docs](http://eqcorrscan.readthedocs.io/en/latest/intro.html#installation)


# Updates

If you want to be kept informed about releases, bug-tracking and enhancements
without having to keep looking on github, subscribe to our [google group](https://groups.google.com/forum/#!forum/eqcorrscan-users).

# Documentation

The full documentation for this package can be found here:
[Docs](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest).
Any errors including typos and just missing bits can either be fixed by you,
or flagged in the issues tab here.  We host our docs on readthedocs, which
uses sphinx to scrape the docstrings in the codes, so it is simple to
match the docs to the codes and change the docstrings.

# Contributing

Please fork this project and work on it there then create a pull request to
merge back to this main repository.  Please create a branch from *develop*.

When you make changes please run the tests in the test directory to ensure
everything merges with minimum effort.  If there is not yet a test to cope
with your changes then please write one.

Please document your functions following the other documentation within the
functions, these doc-scripts will then be built into the main documentation
using Sphinx.

# Functionality

This package contains routines to enable the user to conduct matched-filter earthquake
detections using [obspy](https://github.com/obspy/obspy/wiki) bindings when reading
and writing seismic data, as well as subspace detection, brightness source-scanning,
relative moment calculation using singular-value decomposition,
and correlation pick-adjustment for similar events.

Also within this package are:
* Clustering routines for seismic data;
* Peak finding algorithm (basic, but appropriate for noisy data);
* Automatic amplitude picker for local magnitude scale;
* Obspy.core.event integration, which opens up lots of other functions (Seishub, hypoDDpy etc.);
* Stacking routines including phase-weighted stacking based on Thurber at al. (2014);
* Brightness based template creation based on the work of Frank et al. (2014);
* Singular Value Decomposition derived magnitude calculations based on Rubinstein & Ellsworth (2010).

The code-base has grown to be quite large - it is probably worth
having a look at the docs to check what functions we have.
We are writing a series of tutorials included on the EQcorrscan API
to highlight key functions.

*A note on correlation precision*
*EQcorrscan* computes normalised cross-correlations in the frequency-domain using the
[fftw](www.fftw.org) (Fastest Fourier Transform in the West).  Internally
the C routines enforce double-precision (64-Bit floating point numbers)
for all aspects of the cross-correlations (despite requiring 32-Bit float
input and output). Results in testing are accurate to within ~0.0001 of
time-domain cross-correlation results.

# Test status
Note that tests for travis and appveyor are run daily on master as cron jobs, and may reflect time-out issues.

| Service tests | Badge |
|---------------|-------|
| CI checks | ![test](https://github.com/eqcorrscan/EQcorrscan/workflows/test/badge.svg)
| Code coverage | [![codecov](https://codecov.io/gh/eqcorrscan/EQcorrscan/branch/master/graph/badge.svg)](https://codecov.io/gh/eqcorrscan/EQcorrscan) 

# Licence

This package is written  and maintained by the EQcorrscan developers,
and is distributed under the LGPL GNU License, 
Copyright EQcorrscan developers 2018.
