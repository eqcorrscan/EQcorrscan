What's new
==========

**Note: We are updating the namespace at the moment, you might find that your scripts break, but you should get a helpful error message telling you how to fix it.**

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
