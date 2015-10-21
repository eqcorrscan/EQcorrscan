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
