# EQcorrscan
## A python package to conduct match-filter earthquake detections.

[![Join the chat at https://gitter.im/calum-chamberlain/EQcorrscan](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/calum-chamberlain/EQcorrscan?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![TravisCIStatus](https://travis-ci.org/calum-chamberlain/EQcorrscan.svg?branch=master)](https://travis-ci.org/calum-chamberlain/EQcorrscan)
[![DOI](https://zenodo.org/badge/18852/calum-chamberlain/EQcorrscan.svg)](https://zenodo.org/badge/latestdoi/18852/calum-chamberlain/EQcorrscan)
[![DocumentationStatus](http://readthedocs.org/projects/eqcorrscan/badge/?version=latest)](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest)

# Installation
Installation has been tested on both OSX and Linux (Ubuntu), we currently do not support
Windows systems, but plan to in the future.  Installation for Linux and OS X should be as simple as:

```pip install EQcorrscan```

If upgrading from a previous version, rather than running install --upgrade, I recommend the following:

```pip install -U --no-deps EQcorrscan```

This will not try to upgrade your dependencies, which is not needed for 0.0.8.

*You will likely need sudo permissions to run this command.*

If you have any issues installing please let me know.  You will need to install openCV
seperately using (on Linux):

```apt-get install python-opencv```

Or, for Mac users, this is available on Macports or other similar package managers.

## Updates

If you want to be kept informed about releases, bug-tracking and enhancements
without having to keep looking on github, subscribe to our [google group](https://groups.google.com/forum/#!forum/eqcorrscan-users).

# Documenation

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
Neither of these packages are installed by this software, due to a range of
licenses being implemented.  However, both are open-source and should be installed
before using this package.  This package was written to implement the Matlab routines
used by Chamberlain et al. (2014) for the detection of low-frequency earthquakes.

Also within this package are:
* Clustering routines for seismic data;
* Peak finding algorithm (basic);
* Automatic amplitude picker for local magnitude scale;
* [Seisan](http://seisan.info/) S-file integration for database management and routine earthquake location;
* Obspy.core.event integration, which opens up lots of other functions (Seishub, hypoDDpy etc.);
* Stacking routines including phase-weighted stacking based on Thurber at al. (2014);
* Brightness based template creation based on the work of Frank et al. (2014);
* Singular Value Decomposition derived magnitude calculations based on Rubinstein & Ellsworth (2010).

We are currently hovering around 9,000 lines of code (including doc-strings) - it is probably worth
having a look at the docs to check what functions we have.  We plan to write a series of tutorials to be
included on the EQcorrscan API to highlight key functions, currently our tutorial only shows
how to do the core matched-filter detection.

# Licence

This package is written by Calum Chamberlain of Victoria University of Wellington, and
is distributed under the LGPL GNU License, Copyright Calum Chamberlain 2015.


# Contributing

Please fork this project and work on it there then create a pull request to
merge back into develop.

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
