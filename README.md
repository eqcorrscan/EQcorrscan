# EQcorrscan
A python package to conduct match-filter earthquake detections.

This package contains routines to enable the user to conduct match-filter earthquake
detections using [obspy](https://github.com/obspy/obspy/wiki) bindings when reading
and writing seismic data, and the correlation routine in [openCV](http://opencv.org/).
Neither of these packages are installed by this software, due to a range of
licences being implimented.  However, both are open-source and should be installed
before using this package.  This package was written to impliment the matlab routines
used by Chamberlain et al. (2014) for the detection of low-frequency earthquakes.

Also within this package are:
* Clustering routines for seismic data;
* Peak finding algorithm (basic);
* Automatic amplitude picker for local magnitude scale;
* [Seisan](http://seisan.info/) S-file integration for database management and routine earthquake location;
* Stacking routines including phase-weighted stacking based on Thurber at al. (2014);
* Brightness based template creation based on the work of Frank et al. (2014)

This package is written by Calum Chamberlain of Victoria University of Wellington, and
is distributed under the LGPL GNU Licence, Copyright Calum Chamberlain 2015.

# References
* CJ Chamberlain, DR Shelly, J Townend, TA Stern (2014) [Low‐frequency earthquakes reveal punctuated slow slip on the deep extent of the Alpine Fault, New Zealand](http://onlinelibrary.wiley.com/doi/10.1002/2014GC005436/full), __G-cubed__,doi:10.1002/2014GC005436
* Thurber, C. H., Zeng, X., Thomas, A. M., & Audet, P. (2014). [Phase‐Weighted Stacking Applied to Low‐Frequency Earthquakes](http://www.bssaonline.org/content/early/2014/08/12/0120140077.abstract), __BSSA__, doi:10.1785/0120140077.
* Frank, W. B., & Shapiro, N. M. (2014). [Automatic detection of low-frequency earthquakes (LFEs) based on a beamformed network response](http://gji.oxfordjournals.org/content/197/2/1215.short), __Geophysical Journal International__, 197(2), 1215-1223, doi:10.1093/gji/ggu058.
