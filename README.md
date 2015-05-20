# EQcorrscan
Python package to conduct match-filter earthquake correlation detections

This package conducts a match-filter routine using seismic waveform data as a template.
Within this package are fucntion to generate templates from raw waveform data using
picks from Seisan S-files, functoins to process daylong data using the routines in the
ObsPy project, and a correlation routine using the openCV match_template python
bindings.  

This project makes use of the joblib parallelisation package to run multiple templates
through a single day of siemsic data.  The overarching routine in this code can run
multiple days at once, allowing exploitation of many cores.

This package is written by Calum Chamberlain of Victoria University of Wellington, and
is based on the method outlined in Chamberlain et al. 2014 (G-cubed).

THIS BRANCH IS FOR THE DEVELOPMENT OF THE BRIGHTNESS FUNCTION USED IN FRANK ET AL. 2013
