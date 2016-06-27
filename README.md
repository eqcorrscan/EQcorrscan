# EQcorrscan
## A python package to conduct matched-filter earthquake detections.

[![Join the chat at https://gitter.im/calum-chamberlain/EQcorrscan](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/calum-chamberlain/EQcorrscan?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![TravisCIStatus](https://travis-ci.org/calum-chamberlain/EQcorrscan.svg?branch=master)](https://travis-ci.org/calum-chamberlain/EQcorrscan)
[![Build status](https://ci.appveyor.com/api/projects/status/69bpa53loaq473w7?svg=true)](https://ci.appveyor.com/project/calum-chamberlain/eqcorrscan)
[![Coverage Status](https://coveralls.io/repos/github/calum-chamberlain/EQcorrscan/badge.svg?branch=develop)](https://coveralls.io/github/calum-chamberlain/EQcorrscan?branch=develop)
[![DOI](https://zenodo.org/badge/18852/calum-chamberlain/EQcorrscan.svg)](https://zenodo.org/badge/latestdoi/18852/calum-chamberlain/EQcorrscan)
[![DocumentationStatus](http://readthedocs.org/projects/eqcorrscan/badge/?version=latest)](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest)

# Installation
Installation has been tested on both OSX and Linux (Ubuntu), and now
Windows systems.  We support Python versions 2.7 and 3.5.  The codes likely
work on Py 3.4 too, but we currently don't test this and recommend users to
work in Py 3.5.

Instructions for installing EQcorrscan and the required dependency, openCV
are linked from the [docs](http://eqcorrscan.readthedocs.io/en/latest/intro.html#installation)

*A note on correlation precision*
OpenCV computes cross-correlations in the frequency-domain for normal seismic
datasets (if the dataset is very small then the cross-correlation will be
computed in the time-domain, but this is rare for seismic data).  In testing we
have found that different methods of installing openCV provide different results
for cross-correlations at the very low-end of cross-correlations.  We think this
comes down to how the ffts are computed.  However, for moderate to high cross-correlations
(above 0.05 normalised cross-correlation), all methods provide the same result.

The outcome of this is that for very low thresholds, you may see changes in
your results, however for standard operations this is not an issue.  We have found
that differences are, on average, 0.0024 - which shifts the mean of a single
channel cross-correlation from very close to zero, to 0.0024, and alters the
median.  However we have found that this results in no change in the median
absolute deviation of the data, so thresholds based on this will be the same,
although the cross-correlations themselves will be shifted.  You would have to be
running a very low threshold to see the result of this (0.5 * MAD, rather than
commonly used values around 8 * MAD).

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

This package contains routines to enable the user to conduct match-filter earthquake
detections using [obspy](https://github.com/obspy/obspy/wiki) bindings when reading
and writing seismic data, and the correlation routine in [openCV](http://opencv.org/).
The OpendCV package is not installed by this software, due to a need to build from
source.  The user should follow the instructions above for OpenCV install.
This package was written to implement the Matlab routines
used by Chamberlain et al. (2014) for the detection of low-frequency earthquakes.

Also within this package are:
* Clustering routines for seismic data;
* Peak finding algorithm (basic, but appropriate for noisy data);
* Automatic amplitude picker for local magnitude scale;
* [Seisan](http://seisan.info/) S-file integration for database management and routine earthquake location;
* Obspy.core.event integration, which opens up lots of other functions (Seishub, hypoDDpy etc.);
* Stacking routines including phase-weighted stacking based on Thurber at al. (2014);
* Brightness based template creation based on the work of Frank et al. (2014);
* Singular Value Decomposition derived magnitude calculations based on Rubinstein & Ellsworth (2010).

We are currently hovering around 9,000 lines of code (including doc-strings) - it is probably worth
having a look at the docs to check what functions we have.  We plan to write a series of tutorials to be
included on the EQcorrscan API to highlight key functions, currently our tutorials only show
how to do the core matched-filter detection.

# Licence

This package is written by Calum Chamberlain and Chet Hopp of Victoria University of Wellington, and
is distributed under the LGPL GNU License, Copyright Calum Chamberlain and Chet Hopp 2015, 2016.


# Contributing

Please fork this project and work on it there then create a pull request to
merge back to this main repository.  If you are working on a bug-fix then
use the *develop* branch, otherwise, create a feature branch and work
on your addition there.

When you make changes please run the tests in the test directory to ensure
everything merges with minimum effort.  If there is not yet a test to cope
with your changes then please write one.

Please document your functions following the other documentation within the
functions, these doc-scripts will then be built into the main documentation
using Sphinx.

We are trying to implement a better branching model, following that found [here](http://nvie.com/posts/a-successful-git-branching-model/).
To this end, please fork the development branch if you want to develop
things, and flag issues in the master for us to bugfix.
If you have a feature you want to develop please create a new branch
from the development branch and work on it there, we can then merge
it back in to the development branch when it is stable enough.

This branching model (git-flow) is pretty well established, and I would recommend
you to install [git-flow](https://github.com/nvie/gitflow/wiki/Installation) and
read their [documentation](https://github.com/nvie/gitflow). It seems pretty intuitive and
will keep us all branching in the same way.

# References
* CJ Chamberlain, DR Shelly, J Townend, TA Stern (2014) [Low‐frequency earthquakes reveal punctuated slow slip on the deep extent of the Alpine Fault, New Zealand](http://onlinelibrary.wiley.com/doi/10.1002/2014GC005436/full), __G-cubed__,doi:10.1002/2014GC005436
* Thurber, C. H., Zeng, X., Thomas, A. M., & Audet, P. (2014). [Phase‐Weighted Stacking Applied to Low‐Frequency Earthquakes](http://www.bssaonline.org/content/early/2014/08/12/0120140077.abstract), __BSSA__, doi:10.1785/0120140077.
* Frank, W. B., & Shapiro, N. M. (2014). [Automatic detection of low-frequency earthquakes (LFEs) based on a beamformed network response](http://gji.oxfordjournals.org/content/197/2/1215.short), __Geophysical Journal International__, 197(2), 1215-1223, doi:10.1093/gji/ggu058.
* Rubinstein, J. L., & Ellsworth, W. L. (2010). [Precise estimation of repeating earthquake moment: Example from Parkfield, California](http://www.bssaonline.org/content/100/5A/1952.short), __BSSA__, doi:10.1785/0120100007
fa66e06971b2dd1e4bf766
