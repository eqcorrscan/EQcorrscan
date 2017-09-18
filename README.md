# EQcorrscan
## A python package for the detection and analysis of repeating and near-repeating earthquakes.

# Citation:
We have a manuscript in review, if you make use of EQcorrscan please cite the folloing paper:

Chamberlain, C. J., Hopp, C. J., Boese, C. M., Warren-Smith, E., Chambers, D., Chu, S. X., Michailos, K., Townend, J., EQcorrscan: Repeating and near-repeating earthquake detection and analysis in Python. Seismological Research Letters *in review*

If you want to you should also cite the version number:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.893621.svg)](https://doi.org/10.5281/zenodo.893621)


# Test status
Note that tests for travis and appveyor are run daily on master as cron jobs, and may reflect time-out issues.

| Service tests | Badge |
|---------------|-------|
| OSX & Linux | [![TravisCIStatus](https://travis-ci.org/eqcorrscan/EQcorrscan.svg?branch=master)](https://travis-ci.org/eqcorrscan/EQcorrscan)
| Windows | [![Build status](https://ci.appveyor.com/api/projects/status/b0924mp0uwwyap3d/branch/master?svg=true)](https://ci.appveyor.com/project/calum-chamberlain/eqcorrscan-jsycv/branch/master)
| Code coverage | [![codecov](https://codecov.io/gh/eqcorrscan/EQcorrscan/branch/master/graph/badge.svg)](https://codecov.io/gh/eqcorrscan/EQcorrscan)
| Documentation | [![DocumentationStatus](http://readthedocs.org/projects/eqcorrscan/badge/?version=latest)](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest)
| Dependency status | [![Dependency Status](https://dependencyci.com/github/eqcorrscan/EQcorrscan/badge)](https://dependencyci.com/github/eqcorrscan/EQcorrscan)
| Network tests | [![CircleCI](https://circleci.com/gh/eqcorrscan/EQcorrscan/tree/master.svg?style=svg)](https://circleci.com/gh/eqcorrscan/EQcorrscan/tree/master)
| Issues ready | [![Stories in Ready](https://badge.waffle.io/eqcorrscan/EQcorrscan.png?label=ready&title=Ready)](http://waffle.io/eqcorrscan/EQcorrscan)
| Chat on gitter | [![Join the chat at https://gitter.im/eqcorrscan/EQcorrscan](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/eqcorrscan/EQcorrscan?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# Installation

The easiest way to install EQcorrscan is through anaconda:
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/eqcorrscan/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)

Installation has been tested on both OSX and Linux (Ubuntu), and
Windows systems.  We support Python versions 2.7, 3.4, 3.5 and 3.6.
Note that, although we support Windows, EQcorrscan is optimized for
linux style distributions, and the developers are not extensive Windows
users.

*OSX with gcc-4.9 from homebrew doesn't appear to compile properly, all other gcc versions
 seem to work*

Instructions for installing EQcorrscan and the required dependency, fftw
are linked from the
[docs](http://eqcorrscan.readthedocs.io/en/latest/intro.html#installation)

*A note on correlation precision*
*EQcorrscan* computes normalised cross-correlations in the frequency-domain using the
[fftw](www.fftw.org) (Fastest Fourier Transform in the West).  Internally
the C routines enforce double-precision (64-Bit floating point numbers)
for all aspects of the cross-correlations (despite requiring 32-Bit float
input and output). Results in testing are accurate to within ~0.0001 of
time-domain cross-correlation results.

## Updates

If you want to be kept informed about releases, bug-tracking and enhancements
without having to keep looking on github, subscribe to our [google group](https://groups.google.com/forum/#!forum/eqcorrscan-users).

# Documentation

The full documentation for this package can be found here:
[Docs](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest).
Any errors including typos and just missing bits can either be fixed by you,
or flagged in the issues tab here.  We host our docs on readthedocs, which
uses sphinx to scrape the docstrings in the codes, so it is simple to
match the docs to the codes and change the docstrings.

We also have a github-pages site [EQcorrscan](http://calum-chamberlain.github.io/EQcorrscan/),
which uses jekyll to build the site.  Changes or additions to this site can be made on
the gh-pages branch.

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

# Licence

This package was initially written by Calum Chamberlain and Chet Hopp
of Victoria University of Wellington, and is distributed under the
LGPL GNU License, Copyright EQcorrscan developers 2015, 2016, 2017.


# Contributing

Please fork this project and work on it there then create a pull request to
merge back to this main repository.  Please create a branch from *develop*.

When you make changes please run the tests in the test directory to ensure
everything merges with minimum effort.  If there is not yet a test to cope
with your changes then please write one.

Please document your functions following the other documentation within the
functions, these doc-scripts will then be built into the main documentation
using Sphinx.

# References
* CJ Chamberlain, DR Shelly, J Townend, TA Stern (2014) [Low‐frequency earthquakes reveal punctuated slow slip on the deep extent of the Alpine Fault, New Zealand](http://onlinelibrary.wiley.com/doi/10.1002/2014GC005436/full), __G-cubed__,doi:10.1002/2014GC005436
* Thurber, C. H., Zeng, X., Thomas, A. M., & Audet, P. (2014). [Phase‐Weighted Stacking Applied to Low‐Frequency Earthquakes](http://www.bssaonline.org/content/early/2014/08/12/0120140077.abstract), __BSSA__, doi:10.1785/0120140077.
* Frank, W. B., & Shapiro, N. M. (2014). [Automatic detection of low-frequency earthquakes (LFEs) based on a beamformed network response](http://gji.oxfordjournals.org/content/197/2/1215.short), __Geophysical Journal International__, 197(2), 1215-1223, doi:10.1093/gji/ggu058.
* Rubinstein, J. L., & Ellsworth, W. L. (2010). [Precise estimation of repeating earthquake moment: Example from Parkfield, California](http://www.bssaonline.org/content/100/5A/1952.short), __BSSA__, doi:10.1785/0120100007
fa66e06971b2dd1e4bf766
