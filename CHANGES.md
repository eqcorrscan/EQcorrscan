## 0.1.4
* Bug-fix in plot_repicked removed where data were not normalized
properly;
* Bug-fix in lag_calc where data were missing in the continuous data 
fixed (this led to incorrect picks, **major bug!**);
* Output cross-channel correlation sum in lag-calc output;
* Add id to DETECTION objects, which is consistent with the events
within DETECTION objects and catalog output, and used in lag_calc to
allow linking of detections to catalog events;
* Add lots of logging and error messages to lag-calc to ensure user 
understands limits;
* Add error to day-proc to ensure user is aware of risks of padding;
* Change utils.pre_processing.process to accept different length of 
data enforcement, not just full day (allow for overlap in processing,
which might be useful for reducing day start and end effects);
* Bug-fix in mag_calc.amp_pick_event, broke loop if data were missing;
* Lots of docs adjustment to sort order of doc-strings and hyper-links;
* Allow multiple uses of the same channel in templates (e.g. you can
now use a template with two windows from the same channel, such as a P
and an S);
* Add evaluation mode filter to utils.catalog_utils.filter_picks;
* Update subspace plot to work when detector is not partitioned;
* Make tests run a little faster;
* Add pep8 testing for all code.

## 0.1.3
* Now testing on OSX (python 2.7 and 3.5) - also added linux python 3.4;
* Add lag-calculation and tests for it;
* Change how lag-calc does the trace splitting to reduce memory usage;
* Added pick-filtering utility to clean up tutorials;
* Change template generation function names for clarity (wrappers for
depreciated names);
* Add more useful error messages when picks are not associated with
waveforms;
* Add example plots for more plotting functions;
* Add subspace detector including docs and tutorial.
* Add *delayed* option to all template_gen functions, set to True by
default which retains old behaviour.

## 0.1.2
* Add handling for empty location information in sfiles.
* Added project setup script which creates a useful directory structure and copies
a default match-filter script to the directory.
* Add archive reader helper for default script, and parameter classes and
definitions for default script.
* Re-write history to make repository smaller, removed trash files that had
been added carelessly.
* Now tested on appveyor, so able to be run on Windows.
* Added ability to read hypoDD/tomoDD phase files to obspy events.
* Added simple despiking algorithm - not ideal for correlation as spikes are
interpolated around when found: eqcorrscan.utils.despike.
* Option to output catalog object from match_filter - this will become the
default once we introduce meta-data to templates - currently the picks for
events are the template trace start-times, which will be before the phase-pick
by the lag defined in the template creation - also added event into detection
class, so you can access the event info from the detections, or create a
catalog from a list of detections.
* Add option to extract detections at run-time in match_filter.match_filter.
* Edited multi_event_singlechan to take a catalog with multiple picks, but
requires you to specify the station and channel to plot.
* Add normalize option to stacking routines.
* Add tests for stacking - PWS test needs more checks.
* Add many examples to doc-strings, not complete though.
* Change docs to have one page per function.
* Python 3.5 testing underway, all tests pass, but only testing about 65% of
codebase.
* Add io functions to match_filter to simplify detection handling including
writing detections to catalog and to text file.
* Stricter match_filter testing to enforce exactly the same result with a
variety of systems.
* Add hack to template_gen tutorial to fix differences in sorting between python 3.x
and python 2.
* Added advanced network triggering routine from Konstantinos, allows
different parameters for individual stations - note only uses recursive
sta-lta triggering at the moment.  Useful for template generations alongside
pickers.
* Added magnitude of completeness and b-value calculators to utils.mag_calc

## 0.1.1
* Cope with events not always having time_errors in them in eventtoSfile;
* Convert Quakeml depths from m to km;
* Multiple little fixes to make Sfile conversion play well with GeoNet QuakeML files;
* Add function to convert from obspy.core.inventory.station.Station to string format
for Seisan STATION0.HYP file;
* Merged feature branch - hypoDD into develop, this provides mappings for the
hypoDD location program, including generation of dt.cc files;
* Added tests for functions in catalog_to_dd;
* Implemented unittest tests;
* Changed name of EQcorrscan_plotting to plotting;
* Added depreciation warnings;
* Changed internal structure of pre-processing to aid long-term upkeep;
* Added warnings in docs for template_gen relating to template generation from
set length files;
* Updated template_creation tutorial to use day-long data;
* Renamed Sfile_util to sfile_util, and functions there-in: will warn about name changes;
* Updated template plotting to include pick labels;
* Updated template_creation tutorial to download S-picks as well as P-picks;
* Update sfile_util to cope with many possible unfilled objects;
* Added sac_util to convert from sac headers to useful event information - note,
does not convert all things, just origin and pick times;
* Added from_sac function to template_gen.

## 0.1.0
* Implimented tests for synthetic generation and match-filter functions
* Developed new tutorials that download from GeoNet and are much clearer
* Changed from PICK and EVENTINFO classes to obspy.core.event classes, note
this will break some previous scripts, however wrappers are included for this,
this ended up being the biggy, and is the reason for ceasing the 0.0.x line.
* Added synthetic seismogram generation, very basic seismograms, but they work
as templates for the detection of *some* seismicity.
* Added compatibility with Obspy v.1.0.0.
* All files now follow pep8.
* Removed parameter files completely, and redundant scripts.

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
