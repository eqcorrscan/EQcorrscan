## Current
* core.match_filter.template
  - new quick_group_templates function for 50x quicker template grouping.
* utils.catalog_to_dd._prepare_stream
  - Now more consistently slices templates to length = extract_len * samp_rate
    so that user receives less warnings about insufficient data.
* utils.cluster.decluster_distance_time
  - Bug-fix: fix segmentation fault when declustering more than 46340 detections
    with hypocentral_separation.

## 0.4.4
* core.match_filter
  - Bug-fix: peak-cores could be defined twice in _group_detect through kwargs.
    Fix: only update peak_cores if it isn't there already.
* core.match_filter.tribe
 - Detect now allows passing of pre-processed data
* core.match_filter.template
 - Remove duplicate detections from overlapping windows using `._uniq()`
* core.lag_calc._xcorr_interp
 - CC-interpolation replaced with resampling (more robust), old method
   deprecated. Use new method with use_new_resamp_method=True as **kwarg.
* core.lag_calc:
 - Fixed bug where minimum CC defined via min_cc_from_mean_cc_factor was not
   set correctly for negative correlation sums.
* utils.correlate
 - Fast Matched Filter now supported natively for version >= 1.4.0
 - Only full correlation stacks are returned now (e.g. where fewer than than
   the full number of channels are in the stack at the end of the stack, zeros
   are returned).
* utils.mag_calc.relative_magnitude
 - fixed bug where S-picks / traces were used for relative-magnitude calculation
   against user's choice.
 - implemented full magnitude bias-correction for CC and SNR
* utils.mag_calc.relative_amplitude:
 - returns dicts for SNR measurements
* utils.catalog_to_dd.write_correlations
 - Fixed bug on execution of parallel execution.
 - Added parallel-options for catalog-dt measurements and for stream-preparation
   before cross correlation-dt measurements.
 - Default parallelization of dt-computation is now across events (loads CPUs
   more efficiently), and there is a new option ``max_trace_workers` to use
   the old parallelization strategy across traces.
 - Now includes `all_horiz`-option that will correlate all matching horizontal
   channels no matter to which of these the S-pick is linking.
* utils.clustering
 - Allow to handle indirect comparison of event-waveforms when (i.e., events
   without matching traces which can be compared indirectly via a third event)
 - Allows to set clustering method, metric, and sort_order from
   scipy.cluster.hierarchy.linkage.
* tribe, template, template_gen, archive_read, clustering: remove option to read
  from seishub (deprecated in obspy).

## 0.4.3
* core.match_filter
 - match_filter:
   - Provide option of exporting the cross-correlation sums for additional later
     analysis.
* core.match_filter.party.write
  - BUG-FIX: When `format='tar'` is selected, added a check for .tgz-file
    suffix before checking the filename against an existing file. Previously,
    when a filename without '.tgz'-suffix was supplied, then the file was
    overwritten against the function's intention.
  - Add option `overwrite=True` to allow overwriting of existing files.
* core.match_filter.party.read
  - BUG-FIX: Ensure wildcard reading works as expected: #453
* core.match_filter.party.rethreshold:
  - added option to rethreshold based on absolute values to keep relevant
    detections with large negative detect_val.
* core.lag_calc:
  - Added option to set minimum CC threshold individually for detections based
    on: min(detect_val / n_chans * min_cc_from_mean_cc_factor, min_cc).
  - Added the ability of saving correlation data of the lag_calc.
* core.template_gen:
  - Added support for generating templates from any object with a
    get_waveforms method. See #459.
* utils.mag_calc.calc_b_value:
  - Added useful information to doc-string regarding method and meaning of
    residuals
  - Changed the number of magnitudes used to an int (from a string!?)
* utils.mag_calc.relative_magnitude:
  - Refactor so that `min_cc` is used regardless of whether 
    `weight_by_correlation` is set. See issue #455.
* utils.archive_read
  - Add support for wildcard-comparisons in the list of requested stations and
    channels.
  - New option `arctype='SDS'` to read from a SeisComp Data Structure (SDS).
    This option is also available in `utils.clustering.extract_detections` and
    in `utils.archive_read._check_available_data`.
* utils.catalog_to_dd
  - Bug-fixes in #424:
    - only P and S phases are used now (previously spurious amplitude picks 
      were included in correlations);
    - Checks for length are done prior to correlations and more helpful error
      outputs are provided.
    - Progress is not reported within dt.cc computation
  - `write_station` now supports writing elevations: #424.
* utils.clustering
  - For `cluster`, `distance_matrix` and `cross_chan_correlation`, implemented
    full support for `shift_len != 0`. The latter two functions now return, in
    addition to the distance-matrix, a shift-matrix (both functions) and a
    shift-dictionary (for `distance_matrix`). New option for shifting streams
    as a whole or letting traces shift individually
    (`allow_individual_trace_shifts=True`).
* utils.plotting
  - Function added (twoD_seismplot) for plotting seismicity (#365).
  
## 0.4.2
* Add seed-ids to the _spike_test's message.
* utils.correlation
  - Cross-correlation normalisation errors no-longer raise an error
  - When "out-of-range" correlations occur a warning is given by the C-function
    with details of what channel, what template and where in the data vector
    the issue occurred for the user to check their data.
  - Out-of-range correlations are set to 0.0
  - After extensive testing these errors have always been related to data issues
    within regions where correlations should not be computed (spikes, step
    artifacts due to incorrectly padding data gaps).
  - USERS SHOULD BE CAREFUL TO CHECK THEIR DATA IF THEY SEE THESE WARNINGS
* utils.mag_calc.amp_pick_event
  - Added option to output IASPEI standard amplitudes, with static amplification
    of 1 (rather than 2080 as per Wood Anderson specs).
  - Added `filter_id` and `method_id` to amplitudes to make these methods more
    traceable.
* core.match_filter
  - Bug-fix - cope with data that are too short with `ignore_bad_data=True`.
    This flag is generally not advised, but when used, may attempt to trim all
    data to zero length.  The expected behaviour is to remove bad data and run
    with the remaining data.
  - Party:
    - decluster now accepts a hypocentral_separation argument. This allows
      the inclusion of detections that occur close in time, but not in space.
      This is underwritten by a new findpeaks.decluster_dist_time function
      based on a new C-function.
  - Tribe:
    - Add monkey-patching for clients that do not have a `get_waveforms_bulk`
      method for use in `.client_detect`. See issue #394.
* utils.pre_processing
  - Only templates that need to be reshaped are reshaped now - this can be a lot
    faster.
  
## 0.4.1
* core.match_filter
  - BUG-FIX: Empty families are no longer run through lag-calc when using 
    Party.lag_calc().  Previously this resulted in a "No matching data" error,
    see #341.
* core.template_gen
  - BUG-FIX: Fix bug where events were incorrectly associated with templates
    in `Tribe().construct()` if the given catalog contained events outside
    of the time-range of the stream. See issue #381 and PR #382.
* utils.catalog_to_dd
  - Added ability to turn off parallel processing (this is turned off by 
    default now) for `write_correlations` - parallel processing for moderate
    to large datasets was copying far too much data and using lots of memory.
    This is a short-term fix - ideally we will move filtering and resampling to
    C functions with shared-memory parallelism and GIL releasing.
    See PR #374.
  - Moved parallelism for `_compute_dt_correlations` to the C functions to
    reduce memory overhead. Using a generator to construct sub-catalogs rather
    than making a list of lists in memory. See issue #361.
* utils.mag_calc:
  - `amp_pick_event` now works on a copy of the data by default
  - `amp_pick_event` uses the appropriate digital filter gain to correct the
    applied filter. See issue #376.
  - `amp_pick_event` rewritten for simplicity.
  - `amp_pick_event` now has simple synthetic tests for accuracy.
  - `_sim_wa` uses the full response information to correct to velocity
    this includes FIR filters (previously not used), and ensures that the
    wood-anderson poles (with a single zero) are correctly applied to velocity
    waveforms.
  - `calc_max_curv` is now computed using the non-cumulative distribution.
* Some problem solved in _match_filter_plot. Now it shows all new detections.
* Add plotdir to eqcorrscan.core.lag_calc.lag_calc function to save the images.

## 0.4.0
* Change resampling to use pyFFTW backend for FFT's.  This is an attempt to
  alleviate issue related to large-prime length transforms.  This requires an
  additional dependency, but EQcorrscan already depends on FFTW itself (#316).
* Refactor of catalog_to_dd functions (#322):
  - Speed-ups, using new correlation functions and better resource management
  - Removed enforcement of seisan, arguments are now standard obspy objects.
* Add plotdir to lag-calc, template construction and matched-filter detection
  methods and functions (#330, #325).
* Wholesale re-write of lag-calc function and methods. External interface is
  similar, but some arguments have been depreciated as they were unnecessary (#321).
  - This was done to make use of the new internal correlation functions which
    are faster and more memory efficient.
  - Party.lag_calc and Family.lag_calc now work in-place on the events in 
    the grouping.
  - Added relative_mags method to Party and Family; this can be called from
    lag-calc to avoid reprocessing data.
  - Added lag_calc.xcorr_pick_family as a public facing API to implement
    correlation re-picking of a group of events.
* Renamed utils.clustering.cross_chan_coherence to 
  utils.clustering.cross_chan_correlation to better reflect what it actually 
  does.
* Add --no-mkl flag for setup.py to force the FFTW correlation routines not
  to compile against intels mkl.  On NeSI systems mkl is currently causing
  issues.
* BUG-FIX: `eqcorrscan.utils.mag_calc.dist_calc` calculated the long-way round
  the Earth when changing hemispheres. We now use the Haversine formula, which
  should give better results at short distances, and does not use a flat-Earth
  approximation, so is better suited to larger distances as well.
* Add C-openmp parallel distance-clustering (speed-ups of ~100 times).
* Allow option to not stack correlations in correlation functions.
* Use compiled correlation functions for correlation clustering (speed-up).
* Add time-clustering for catalogs and change how space-time cluster works
  so that it uses the time-clustering, rather than just throwing out events
  outside the time-range.
* Changed all prints to calls to logging, as a result, debug is no longer
  an argument for function calls.
* `find-peaks` replaced by compiled peak finding routine - more efficient
  both in memory and time #249 - approx 50x faster
  * Note that the results of the C-func and the Python functions are slightly
    different.  The C function (now the default) is more stable when peaks
    are small and close together (e.g. in noisy data).
* multi-find peaks makes use of openMP parallelism for more efficient
  memory usage #249
* enforce normalization of continuous data before correlation to avoid float32
  overflow errors that result in correlation errors (see pr #292).
* Add SEC-C style chunked cross-correlations.  This is both faster and more
  memory efficient.  This is now used by default with an fft length of
  2 ** 13.  This was found to be consistently the fastest length in testing.
  This can be changed by the user by passing the `fft_len` keyword argument.
  See PR #285.
* Outer-loop parallelism has been disabled for all systems now. This was not
  useful in most situations and is hard to maintain.
* Improved support for compilation on RedHat systems
* Refactored match-filter into smaller files. Namespace remains the same.
  This was done to ease maintenance - the match_filter.py file had become
  massive and was slow to load and process in IDEs.
* Refactored `_prep_data_for_correlation` to reduce looping for speed, 
  now approximately six times faster than previously (minor speed-up)
  * Now explicitly doesn't allow templates with different length traces - 
    previously this was ignored and templates with different length 
    channels to other templates had their channels padded with zeros or 
    trimmed.
* Add `skip_short_channels` option to template generation.  This allows users 
  to provide data of unknown length and short channels will not be used, rather
  than generating an error. This is useful for downloading data from 
  datacentres via the `from_client` method.
* Remove pytest_namespace in conftest.py to support pytest 4.x
* Add `ignore_bad_data` kwarg for all processing functions, if set to True
  (defaults to False for continuity) then any errors related to bad data at 
  process-time will be supressed and empty traces returned.  This is useful 
  for downloading data from  datacentres via the `from_client` method when
  data quality is not known.
* Added relative amplitude measurements as
  `utils.mag_calc.relative_amplitude` (#306).
* Added relative magnitude calculation using relative amplitudes weighted by
  correlations to `utils.mag_calc.relative_magnitude`.
* Added `relative_magnitudes` argument to 
  `eqcorrscan.core.match_filter.party.Party.lag_calc` to provide an in-flow
  way to compute relative magnitudes for detected events.
* Events constructed from detections now include estimated origins alongside
  the picks. These origins are time-shifted versions of the template origin and
  should be used with caution. They are corrected for prepick (#308).
* Picks in detection.event are now corrected for prepick *if* the template is
  given. This is now standard in all Tribe, Party and Family methods. Picks will
  not be corrected for prepick in match_filter (#308).
* Fix #298 where the header was repeated in detection csv files. Also added
  a `write_detections` function to `eqcorrscan.core.match_filter.detection`
  to streamline writing detections.
* Remove support for Python 2.7.
* Add warning about unused data when using `Tribe.detect` methods with data that
  do not fit into chunks. Fixes #291.
* Fix #179 when decimating for cccsum_hist in `_match_filter_plot`
* `utils.pre_processing` now uses the `.interpolate` method rather than
  `.resample` to change the sampling rate of data. This is generally more
  stable and faster than resampling in the frequency domain, but will likely
  change the quality of correlations.
* Removed depreciated `template_gen` functions and `bright_lights` and
  `seismo_logs`. See #315
* BUG-FIX: `eqcorrscan.core.template_gen.py` fix conflict with special character on windows
  output-filename. See issue #344

## 0.3.3
* Make test-script more stable.
* Fix bug where `set_xcorr` as context manager did not correctly reset
  stream_xcorr methods.
* Correct test-script (`test_eqcorrscan.py`) to find paths properly.
* BUG-FIX in `Party.decluster` when detections made at exactly the same time
  the first, rather than the highest of these was taken.
* Catch one-sample difference in day properly in pre-processing.dayproc
* Shortproc now clips and pads to the correct length asserted by starttime and
  endtime.
* Bug-fix: Match-filter collection objects (Tribe, Party, Family) implemented
  addition (`__add__`) to alter the main object. Now the main object is left
  unchanged.
* `Family.catalog` is now an immutable property.

## 0.3.2
* Implement reading Party objects from multiple files, including wildcard
  expansion. This will only read template information if it was not 
  previously read in (which is a little more efficient).
* Allow reading of Party objects without reading the catalog files.
* Check quality of downloaded data in `Tribe.client_detect()` and remove it if it
  would otherwise result in errors.
* Add `process_cores` argument to `Tribe.client_detect()` and `Tribe.detect()`
  to provide a separate number of cores for processing and peak-finding - both
  functions are less memory efficient that fftw correlation and can result in
  memory errors if using lots of cores.
* Allow passing of `cores_outer` kwarg through to fftw correlate functions to
  control inner/outer thread numbers. If given, `cores` will define the number
  of inner-cores (used for parallel fft calculation) and `cores_outer` sets
  the number of channels to process in parallel (which results in increased
  memory usage).
* Allow Tribe and Party IO to use QUAKEML or SC3ML format for catalogs (NORDIC
  to come once obspy updates).
* Allow Party IO to not write detection catalogs if so desired, because 
  writing and reading large catalogs can be slow.
* If detection-catalogs are not read in, then the detection events will be
  generated on the fly using `Detection._calculate_event`.
* BUG-FIX: When one template in a set of templates had a channel repeated,
  all detections had an extra, spurious pick in their event object. This
  should no-longer happen.
* Add `select` method to `Party` and `Tribe` to allow selection of a 
  specific family/template.
* Use a compiled C peak-finding function instead of scipy ndimage - speed-up
  of about 2x in testing.
* BUG-FIX: When `full_peaks=True` for `find_peaks2_short` values that were not
  above their neighbours were returned. Now only values greater than their two
  neighbours are returned.
* Add ability to "retry" downloading in `Tribe.client_detect`.
* Change behaviour of template_gen for data that are daylong, but do not start
  within 1 minute of a day-break - previous versions enforced padding to
  start and end at day-breaks, which led to zeros in the data and undesirable 
  behaviour.
* BUG-FIX: Normalisation errors not properly passed back from internal fftw
  correlation functions, gaps not always properly handled during long-period
  trends - variance threshold is now raised, and Python checks for low-variance
  and applies gain to stabilise correlations if needed.
* Plotting functions are now tested and have a more consistent interface:
  * All plotting functions accept the keyword arguments `save`, `savefile`,
    `show`, `return_figure` and `title`.
  * All plotting functions return a figure.
  * `SVD_plot` renamed to `svd_plot`
* Enforce pre-processing even when no filters or resampling is to be done
  to ensure gaps are properly processed (when called from `Tribe.detect`,
  `Template.detect` or `Tribe.client_detect`)
* BUG-FIX in `Tribe.client_detect` where data were processed from data 
  one sample too long resulting in minor differences in data processing
  (due to difference in FFT length) and therefore minor differences 
  in resulting correlations (~0.07 per channel).
  * Includes extra stability check in fftw_normxcorr which affects the
    last sample before a gap when that sample is near-zero.
* BUG-FIX: fftw correlation dot product was not thread-safe on some systems.
  The dot-product did not have the inner index protected as a private variable.
  This did not appear to cause issues for Linux with Python 3.x or Windows, but
  did cause issues for on Linux for Python 2.7 and Mac OS builds.
* KeyboardInterrupt (e.g. ctrl-c) should now be caught during python parallel
  processes.
* Stopped allowing outer-threading on OSX, clang openMP is not thread-safe
  for how we have this set-up. Inner threading is faster and more memory
  efficient anyway.
* Added testing script (`test_eqcorrscan.py`, which will be installed to your
  path on installation of EQcorrscan) that will download all the relevant 
  data and run the tests on the installed package - no need to clone 
  EQcorrscan to run tests!


## 0.3.1
* Cleaned imports in utils modules
* Removed parallel checking loop in archive_read.
* Add better checks for timing in lag-calc functions (#207)
* Removed gap-threshold of twice the template length in `Tribe.client_detect`, see
  issue #224.
* Bug-fix: give multi_find_peaks a cores kwarg to limit thread
  usage.
* Check for the same value in a row in continuous data when computing
  correlations and zero resulting correlations where the whole window
  is the same value repeated (#224, #230).
* BUG-FIX: template generation `from_client` methods for swin=P_all or S_all
  now download all channels and return them (as they should). See #235 and #206
* Change from raising an error if data from a station are not long enough, to
  logging a critical warning and not using the station.
* Add ability to give multiple `swin` options as a list. Remains backwards
  compatible with single `swin` arguments.
* Add option to `save_progress` for long running `Tribe` methods. Files
  are written to temporary files local to the caller.
* Fix bug where if gaps overlapped the endtime set in pre_processing an error
  was raised - happened when downloading data with a deliberate pad at either
  end.

## 0.3.0
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

## 0.2.7
* Patch multi_corr.c to work with more versions of MSVC;
* Revert to using single-precision floats for correlations (as in previous,
  < 0.2.x versions) for memory efficiency.

## 0.2.6
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


## 0.2.5
* Fix bug with \_group_process that resulted in stalled processes.
* Force NumPy version
* Support indexing of Tribe and Party objects by template-name.
* Add tests for lag-calc issue with preparing data
* Change internals of *eqcorrscan.core.lag_calc._prepare_data* to use a
  dictionary for delays, and to work correctly! Issues arose from not checking
  for masked data properly and not checking length properly.
* Fix bug in match_filter.match_filter when checking for equal length traces,
  length count was one sample too short.

## 0.2.4
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
        
## 0.2.0-0.2.3
* See 0.2.4: these versions were not fully released while trying to get
  anaconda packages to build properly.

## 0.1.6
* Fix bug introduced in version 0.1.5 for match_filter where looping
through multiple templates did not correctly match image and template
data: 0.1.5 fix did not work;
* Bug-fix in catalog_to_dd for events without magnitudes;
* Amend match-filter to not edit the list of template names in place.
Previously, if a template was not used (due to no matching continuous
data) then the name of the template was removed: this now copies the
list of template_names internally and does not change the external list.

## 0.1.5
* Migrate coverage to codecov;
* Fix bug introduced in version 0.1.5 for match_filter where looping
through multiple templates did not correctly match image and template
data.

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
* Implemented tests for synthetic generation and match-filter functions
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
