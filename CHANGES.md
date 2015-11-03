## 0.0-a.1 First alpha release
* First release to allow students to use the core functions
* Not fully functional, bigs being found daily
* match_filter functions work well - could be improved by using openCV matchTemplate ability to work with multiple templates, which could be optimized on GPUs.

## 0.0-a.2
* Bug fixes to Sfile_util when writing S-files - converting from milliseconds to decimal seconds
* Fix catalogue2DD to weight by correlation value
* Fix lagging in bright_lights
* Now tested on supercomputer (cluster) computers

## 0.0.3
* Many small bug-fixes
* Added greater parallel ability
* Change directory structure to be more like a true python package - now pip installable

## 0.0.4
* Travis.CI integration implemented
* Tests run (needs more tests)
* Now has synthetic template generation (major step, very much in alpha)

## 0.0.5
* Update all paths in functions when calling EQcorrscan functions to use the installed version
* Remove parameter file usage in core functions and replace with variables, parameter files remain for scripts

## 0.0.6 & 0.0.7
Note, double release due to failed 0.0.6 release
* Properly installable via pip

## Development
* Added SVD magnitude inversion to utils.mag_calc, tested for multi-channel
* Bug-fix in s-file printing when printing AIN, convert to int now before print
* Add master option to clustering.align_traces, alignment can be forced to this
* Add plot directory output option to match_filter
* Change plot downsampling to 10 Hz
* Clustering.cluster now output groups properly when computing clustering by cross-correlation
* Add plot_synth_real function to EQcorrscan_plotting - use for plotting a synthetic template on top of a real template
* Add space-time clustering fucntion to clustering, use to group repeating events
* Add re_thresh_csv to clustering, can be used to increase the detection threshold after detection run
