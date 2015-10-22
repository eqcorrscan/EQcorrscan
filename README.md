# EQcorrscan
A python package to conduct match-filter earthquake detections.

[![TravisCI Status](https://travis-ci.org/calum-chamberlain/EQcorrscan.svg?branch=develop)](https://travis-ci.org/calum-chamberlain/EQcorrscan)
[![DOI](https://zenodo.org/badge/18852/calum-chamberlain/EQcorrscan.svg)](https://zenodo.org/badge/latestdoi/18852/calum-chamberlain/EQcorrscan)
[![Documentation Status](http://readthedocs.org/projects/eqcorrscan/badge/?version=latest)](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest)

If you are looking for a (relatively) stable release, please download the latest
[release](https://github.com/calum-chamberlain/EQcorrscan/releases/tag/v0.0-alpha.1).
This is currently the first alpha release and things will change rapidly, however
changes will be implemented first to the master branch, then pushed to the release
and the release number updated.  Please raise issues on this repository, preferably
with helpful comments, like line numbers, the command run and the error messages.
Alternatively, fork the master and help me develop this!

The full documentation for this package can be found here:
[EQcorrscan](http://calum-chamberlain.github.io/EQcorrscan/) - you are also
welcome to suggest documentation updates or update the doc in the master branch, please
do not work on the documentation in the gh-pages branch!

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
* Stacking routines including phase-weighted stacking based on Thurber at al. (2014);
* Brightness based template creation based on the work of Frank et al. (2014)

This package is written by Calum Chamberlain of Victoria University of Wellington, and
is distributed under the LGPL GNU License, Copyright Calum Chamberlain 2015.

# Installation
If you understand that this package is in **heavy** development, but don't plan
on contributing to your own fork, or if you just want to give this package a spin,
then this package can be installed by pip:

*pip install git+https://github.com/calum-chamberlain/EQcorrscan.git*

You will likely need sudo permissions to run this command.  This installation
method is quite new to the package (as of v0.0.4), so there are some things that
are changing in the install process at the moment.  This should be smoothed out
by v0.1.0 (maybe mid 2016).

# Parameter files
To use this package you will need to set up default parameters in the parameter
file. Currently these are located in the source directory, this will change in
the future, hopefully. It is recommended that you copy these default parameter
files before adding your own to allow you to easily transfer back to other
parameter set ups.

# Contributing
Please fork this project and work on it there then create a pull request to
merge back into the master.

When you make changes please run the tests in the test directory to ensure
everything merges with minimum effort.

Please document your functions following the other documentation within the
functions, these doc-scripts will then be built into the main documentation
using Sphinx.

We are trying to implement a better branching model, following that found here:
http://nvie.com/posts/a-successful-git-branching-model/
To this end, please fork the development branch if you want to develop
things, and flag issues in the master for us to bugfix.
If you have a feature you want to develop please create a new branch
of the development branch for this and work in there, we can then merge
it back in to the development branch when it is stable enough.

This branching model (git-flow) is pretty well established, and I would recommend
you to install [git-flow](https://github.com/nvie/gitflow/wiki/Installation) and
read their [docs](https://github.com/nvie/gitflow). It seems pretty intuitive and
will keep us all branching in the same way.

# References
* CJ Chamberlain, DR Shelly, J Townend, TA Stern (2014) [Low‐frequency earthquakes reveal punctuated slow slip on the deep extent of the Alpine Fault, New Zealand](http://onlinelibrary.wiley.com/doi/10.1002/2014GC005436/full), __G-cubed__,doi:10.1002/2014GC005436
* Thurber, C. H., Zeng, X., Thomas, A. M., & Audet, P. (2014). [Phase‐Weighted Stacking Applied to Low‐Frequency Earthquakes](http://www.bssaonline.org/content/early/2014/08/12/0120140077.abstract), __BSSA__, doi:10.1785/0120140077.
* Frank, W. B., & Shapiro, N. M. (2014). [Automatic detection of low-frequency earthquakes (LFEs) based on a beamformed network response](http://gji.oxfordjournals.org/content/197/2/1215.short), __Geophysical Journal International__, 197(2), 1215-1223, doi:10.1093/gji/ggu058.
