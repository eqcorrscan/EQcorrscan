## Development:
* Added SVD magnitude inversion to [utils.mag_calc](EQcorrscan/eqcorrscan/utils/mag_calc.py#L530),
tested for multi-channel;
* Bug-fix in s-file printing when printing AIN,
[convert to int](EQcorrscan/eqcorrscan/utils/Sfile_util.py#L75) now before print;
* Add master option to [stacking.align_traces](EQcorrscan/eqcorrscan/utils/stacking.py#L93),
alignment can be forced to this;
* Add plot directory output option to [match_filter](EQcorrscan/eqcorrscan/core/match_filter.py#L341);
* Change plot downsampling to 10 Hz in [match_filter](EQcorrscan/eqcorrscan/core/match_filter.py#L489);
* [Clustering.cluster](EQcorrscan/eqcorrscan/utils/clustering.py#L81)
now output groups properly when
computing clustering by cross-correlation;
* Add [plot_synth_real](EQcorrscan/eqcorrscan/utils/EQcorrscan_plotting.py#L765)
function to EQcorrscan_plotting -
use for plotting a synthetic template on top of a real template;
* Add [space-time](EQcorrscan/eqcorrscan/utils/clustering.py#L513)
clustering fucntion to clustering,
use to group repeating events;
* Add [re_thresh_csv](EQcorrscan/eqcorrscan/utils/clustering.py#L551)
to clustering, can be used to increase
the detection threshold after detection run;
* Add mock modules to conf.py for ReadTheDocs, also removed
requirements.txt file as it was unused;
* Added parallel options to distance_matrix computation in clustering.
* Change matplotlib import location in match_filter to allow
other functions to be called without turning off interactive plots.
* Add **both** option to utils.synth_seis to allow individual creation
of both P and S phases.
* Add sanity check for nyquist frequency to utils.pre_processing, now
breaks if highcut >= samp_rate
* Add plot_format option to core.match_filter
* Add multiple phase capability to utils.synth_seis

:volcano:

## 0.0.6 & 0.0.7;
*Note, double release due to failed 0.0.6 release*
* Properly installable via pip.

## 0.0.5:
* Update all paths in functions when calling EQcorrscan
functions to use the installed version;
* Remove parameter file usage in core functions and
replace with variables, parameter files remain for scripts.

## 0.0.4:
* Travis.CI integration implemented;
* Tests run (needs more tests);
* Now has synthetic template generation (major step, very much in alpha).

## 0.0.3:
* Many small bug-fixes;
* Added greater parallel ability;
* Change directory structure to be more like a true python
package - now pip installable.

## 0.0-a.2:
* Bug fixes to Sfile_util when writing S-files - converting
from milliseconds to decimal seconds;
* Fix catalogue2DD to weight by correlation value;
* Fix lagging in bright_lights;
* Now tested on supercomputer (cluster) computers.

## 0.0-a.1 First alpha release:
* First release to allow students to use the core functions;
* Not fully functional, bugs being found daily;
* match_filter functions work well - could be improved by using
openCV matchTemplate ability to work with multiple templates,
which could be optimised on GPUs **GPU capability not applicable as far as I know**
