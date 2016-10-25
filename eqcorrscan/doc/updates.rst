What's new
==========

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
