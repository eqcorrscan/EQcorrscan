What's new
==========

Version 0.3.0
-------------
* Compiled peak-finding routine written to speed-up peak-finding.
* Change default match-filter plotting to not decimate unless it has to.
* BUG-FIX: changed minimum variance for fftw correlation backend.
* Do not try to process when no processing needs to be done in 
  core.match_filter._group_process.
* Length checking in core.match_filter._group_process done in samples rather
  than time.
* BUG-FIX: Fix bug where data lengths were not correct in 
  match_filter.Tribe.detect when sampling time-stamps were inconsistent between
  channels, which previously resulted in error.
* BUG-FIX: Fix memory-leak in tribe.construct
* Add plotting options for plotting rate to Party.plot
* Add filtering detections by date as Party.filter
* BUG-FIX: Change method for Party.rethreshold: list.remove was not reliable.
* Add option `full_peaks` to detect methods to map to find_peaks.
* pre-processing (and match-filter object methods) are now gap-aware and will
  accept gappy traces and can return gappy traces. By default gaps are filled to
  maintain backwards compatibility. Note that the fftw correlation backend
  requires gaps to be padded with zeros.
* **Removed sfile_utils** This support for Nordic IO has been upgraded and moved
  to obspy for obspy version 1.1.0.  All functions are there and many bugs have
  been fixed. This also means the removal of nordic-specific functions in
  EQcorrscan - the following functions have been removed:
  * template_gen.from_sfile
  * template_gen.from_contbase
  * mag_calc.amp_pick_sfile
  * mag_calc.pick_db
  All removed functions will error and tell you to use obspy.io.nordic.core.
  This now means that you can use obspy's `read_events` to read in sfiles.
* Added `P_all` and `S_all` options to template generation functions
  to allow creation of multi-channel templates starting at the P and S
  times respectively.
* Refactored `template_gen`, all options are available via 
  `template_gen(method=...)`, and depreciation warnings are in place.
* Added some docs for converting older templates and detections into Template
  and Party objects.

Version 0.2.7
-------------
* Patch multi_corr.c to work with more versions of MSVC;
* Revert to using single-precision floats for correlations (as in previous,
  < 0.2.x versions) for memory efficiency.

Version 0.2.6
-------------
* Added the ability to change the correlation functions used in detection
  methods through the parameter xcorr_func of match_filter, Template.detect
  and Tribe.detect, or using the set_xcorr context manager in
  the utils.correlate module. Supported options are:
    * numpy
    * fftw
    * time-domain
    * or passing a function that implements the xcorr interface.
* Added the ability to change the concurrency strategy of xcorr functions
  using the paramter concurrency of match_filter, Template.detect
  and Tribe.detect. Supported options are:
    * None - for single-threaded execution in a single process
    * multithread - for multi-threaded execution
    * multiprocess- for multiprocess execution
    * concurrent - allows functions to describe their own preferred currency
    methods, defaults to multithread
* Change debug printing output, it should be a little quieter;
* Speed-up time-domain using a threaded C-routine - separate from frequency
  domain C-routines;
* Expose useful parallel options for all correlation routines;
* Expose cores argument for match-filter objects to allow limits to be placed
  on how much of your machine is used;
* Limit number of workers created during pre-processing to never be more than
  the number of traces in the stream being processed;
* Implement openMP parallelisation of cross-correlation sum routines - memory
  consumption reduced by using shared memory, and by computing the
  cross-correlation sums rather than individual channel cross-correlations.
  This also leads to a speed-up.  This routine is the default concurrent
  correlation routine;
* Test examples in rst doc files to ensure they are up-to-date;
* Tests that were prone to timeout issues have been migrated to run on circleci
  to allow quick re-starting of fails not due to code errors


Version 0.2.5
-------------
* Fix bug with \_group_process that resulted in stalled processes.
* Force NumPy version
* Support indexing of Tribe and Party objects by template-name.
* Add tests for lag-calc issue with preparing data
* Change internals of *eqcorrscan.core.lag_calc._prepare_data* to use a
  dictionary for delays, and to work correctly! Issues arose from not checking
  for masked data properly and not checking length properly.
* Fix bug in match_filter.match_filter when checking for equal length traces,
  length count was one sample too short.

Version 0.2.4
-------------
* Increase test coverage (edge-cases) in template_gen;
* Fix bug in template_gen.extract_from_stack for duplicate channels in
template;
* Increase coverage somewhat in bright_lights, remove non-parallel
option (previously only used for debugging in development);
* Increase test coverage in lag_calc;
* Speed-up tests for brightness;
* Increase test coverage for match_filter including testing io of
detections;
* Increase subspace test coverage for edge cases;
* Speed-up catalog_to_dd_tests;
* Lag-calc will pick S-picks on channels ending E, N, 1 and 2, change
from only picking on E and N before; warning added to docs;
* Add full tests for pre-processing;
* Run tests in parallel on ci, speed-up tests dramatically;
* Rename singular-value decomposition functions (with depreciation
warnings);
* Rename SVD_moments to lower-case and add depreciation warning;
* Increase test coverage in utils.mag_calc;
* Add Template, Tribe, Family, Party objects and rename DETECTION to
Detection;
    * Template objects maintain meta-data associated with their creation
    to stream-line processing of data (e.g. reduce chance of using the
    wrong filters).
    * Template events have a detect method which takes unprocessed data
    and does the correct processing using the Template meta-data, and
    computes the matched-filter detections.
    * Tribe objects are containers for multiple Templates.
    * Tribe objects have a detect method which groups Templates with
    similar meta-data (processing information) and runs these templates
    in parallel through the matched-filter routine. Tribe.detect outputs
    a Party of Family objects.
    * The Party object is a container for many Family objects.
    * Family objects are containers for detections from the same
    Template.
    * Family and Party objects have a lag_calc method which computes
    the cross-correlation pick-refinements.
    * The upshot of this is that it is possible to, in one line,
    generate a Tribe of templates, compute their matched-filter
    detections, and generate cross-correlation pick refinements, which
    output Event objects, which can be written to a catalog:
        Tribe.construct(method, **kwargs).detect(st, **kwargs).lag_calc(**kwargs).write()
    * Added 25 tests for these methods.
    * Add parameters *threshold_type* and *threshold_input* to Detection
    class.  Add support for legacy Detection objects via NaN and unset
    values.
* Removed support for obspy < 1.0.0
* Update / correct doc-strings in template-gen functions when describing
processing parameters.
* Add warning message when removing channels from continuous data in
match_filter;
* Add min_snr option for template generation routines, if the
signal-to-noise ratio is below a user-defined threshold, the channel
will not be used.
* Stop enforcing two-channel template channel names.
* Fix bug in detection_multiplot which didn't allow streams with
fewer traces than template;
* Update internals to custom C fftw-based correlation rather than openCV (Major change);
    * OpenCV has been removed as a dependancy;
    * eqcorrscan.core.match_filter.normxcorr2 now calls a compiled C routine;
    * Parallel workflows handled by openMP rather than Python Multiprocessing
      for matched-filter operations to allow better memory handling.
        * It is worth noting that we tried re-writing using SciPy internals
        which led to a significant speed-up, but with high memory costs,
        we ended up going with this option, which was the more difficult
        option, because it allows effective use on SLURM managed systems
        where python multiprocessing results in un-real memory spikes
        (issue #88).
        
Version 0.2.0-0.2.3
-------------------
* See 0.2.4: these versions were not fully released while trying to get
  anaconda packages to build properly.

Version 0.1.6
-------------
* Fix bug introduced in version 0.1.5 for match_filter where looping
through multiple templates did not correctly match image and template
data: 0.1.5 fix did not work;
* Bug-fix in catalog_to_dd for events without magnitudes;
* Amend match-filter to not edit the list of template names in place.
Previously, if a template was not used (due to no matching continuous
data) then the name of the template was removed: this now copies the
list of template_names internally and does not change the external list.

Version 0.1.5
-------------
* Migrate coverage to codecov;
* Fix bug introduced in version 0.1.5 for match_filter where looping
through multiple templates did not correctly match image and template
data.

Version 0.1.4
-------------
* Bug-fix in plot_repicked removed where data were not normalized properly;
* Bug-fix in lag_calc where data were missing in the continuous data fixed (this led to incorrect picks, **major bug!**);
* Output cross-channel correlation sum in lag-calc output;
* Add id to DETECTION objects, which is consistent with the events within DETECTION objects and catalog output, and used in lag_calc to allow linking of detections to catalog events;
* Add lots of logging and error messages to lag-calc to ensure user understands limits;
* Add error to day-proc to ensure user is aware of risks of padding;
* Change utils.pre_processing.process to accept different length of data enforcement, not just full day (allow for overlap in processing, which might be useful for reducing day start and end effects);
* Bug-fix in mag_calc.amp_pick_event, broke loop if data were missing;
* Lots of docs adjustment to sort order of doc-strings and hyper-links;
* Allow multiple uses of the same channel in templates (e.g. you can now use a template with two windows from the same channel, such as a P and an S);
* Add evaluation mode filter to utils.catalog_utils.filter_picks;
* Update subspace plot to work when detector is not partitioned;
* Make tests run a little faster;
* Add pep8 testing for all code.


Version 0.1.3
-------------
* Now testing on OSX (python 2.7 and 3.5) - also added linux python 3.4;
* Add lag-calculation and tests for it;
* Change how lag-calc does the trace splitting to reduce memory usage;
* Added pick-filtering utility to clean up tutorials;
* Change template generation function names for clarity (wrappers for depreciated names);
* Add more useful error messages when picks are not associated with waveforms;
* Add example plots for more plotting functions;
* Add subspace detector including docs and tutorial.
* Add *delayed* option to all template_gen functions, set to True by default which retains old behaviour.


Version 0.1.2
-------------

* Add handling for empty location information in sfiles;
* Added project setup script which creates a useful directory structure and copies a default match-filter script to the directory;
* Add archive reader helper for default script, and parameter classes and definitions for default script;
* Re-write history to make repository smaller, removed trash files that had been added carelessly;
* Now tested on appveyor, so able to be run on Windows;
* Added ability to read hypoDD/tomoDD phase files to obspy events;
* Added simple despiking algorithm - not ideal for correlation as spikes are interpolated around when found: eqcorrscan.utils.despike;
* Option to output catalog object from match_filter - this will become the default once we introduce meta-data to templates - currently the picks for events are the template trace start-times, which will be before the phase-pick by the lag defined in the template creation - also added event into detection class, so you can access the event info from the detections, or create a catalog from a list of detections;
* Add option to extract detections at run-time in match_filter.match_filter;
* Edited multi_event_singlechan to take a catalog with multiple picks, but requires you to specify the station and channel to plot;
* Add normalize option to stacking routines;
* Add tests for stacking - PWS test needs more checks;
* Add many examples to doc-strings, not complete though;
* Change docs to have one page per function.
* Python 3.5 testing underway, all tests pass, but only testing about 65% of codebase.
* Add io functions to match_filter to simplify detection handling including writing detections to catalog and to text file.
* Stricter match_filter testing to enforce exactly the same result with a variety of systems.
* Add hack to template_gen tutorial to fix differences in sorting between python 3.x and python 2.
* Added advanced network triggering routine from Konstantinos, allows different parameters for individual stations - note only uses recursive sta-lta triggering at the moment.  Useful for template generations alongside pickers.
* Added magnitude of completeness and b-value calculators to utils.mag_calc

Version 0.1.1
-------------

* Cope with events not always having time_errors in them in eventtoSfile;
* Convert Quakeml depths from m to km;
* Multiple little fixes to make Sfile conversion play well with GeoNet QuakeML files;
* Add function to convert from obspy.core.inventory.station.Station to string format for Seisan STATION0.HYP file;
* Merged feature branch - hypoDD into develop, this provides mappings for the hypoDD location program, including generation of dt.cc files;
* Added tests for functions in catalog_to_dd;
* Implemented unittest tests;
* Changed name of EQcorrscan_plotting to plotting;
* Added depreciation warnings;
* Changed internal structure of pre-processing to aid long-term upkeep;
* Added warnings in docs for template_gen relating to template generation from set length files;
* Updated template_creation tutorial to use day-long data;
* Renamed Sfile_util to sfile_util, and functions there-in: will warn about name changes;
* Updated template plotting to include pick labels;
* Updated template_creation tutorial to download S-picks as well as P-picks;
* Update sfile_util to cope with many possible unfilled objects;
* Added sac_util to convert from sac headers to useful event information - note, does not convert all things, just origin and pick times;
* Added from_sac function to template_gen.
