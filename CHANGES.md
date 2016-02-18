## 0.1.0
* Implimented tests for synthetic generation and match-filter functions
* Developed new tutorials that download from GeoNet and are much clearer
* Changed from PICK and EVENTINFO classes to obspy.core.event classes, note
this will break some previous scripts, however wrappers are included for this
* Added synthetic seismogram generation, very basic seismograms, but they work
as templates for the detection of *some* seismicity.

## 0.0.9
* Working towards following pep8
* Change how match_filter handles missing data - now remove template traces (which have been copied)
rather than adding null traces to the data stream - more memory efficient, and faster
* Change float handling for large amplitudes in Sfile_utils
* Change distance decimal handling in Sfile_utils
* Add magnitude-frequency plotting to EQcorrscan_plotting
* Update tutorial and docs for tutorial
* **Major bug-fix** Change in match_filter - debug and index were the wrong way
round from version 0.0.5 onwards, hence for multiple templates, cross-correlation
vectors were not matched for templates. Will push a release because of this.

## 0.0.8:
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
* **BUG-FIX** Change match-filter data stream handling to copy the stream and
keep it safe before adding null-traces or removing excess traces.  Match_filter
will now remove excess, un-needed traces from the copied stream.  This is specifically
necessary when splitting a large template set, with a range of channels, into smaller
groups to be run serially.

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
