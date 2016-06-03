Utils
=====

Utility functions for integration with other software,
and for the analysis of waveforms detected by cross-correlation.

Sfile_util allows for integration with Seisan via event and pick information
writing to nordic format s-files.  Not all information is currently converted,
however pick and origin information are.  We now internally store pick and
event information in the obspy.core.event classes, which support writing and
reading to and from multiple formats, including quakeML and NonLinLoc.
The authors recommend using this integration with NonLinLoc for location
of events.  Because NonLinLoc files can be read in by obspy, and then written
out to nordic format files, which can be easily converted to hypoDD files,
there are now multiple options for location of events.

.. toctree::
   :maxdepth: 1

   submodules/utils.archive_read
   submodules/utils.clustering
   submodules/utils.plotting
   submodules/utils.findpeaks
   submodules/utils.despike
   submodules/utils.picker
   submodules/utils.mag_calc
   submodules/utils.pre_processing
   submodules/utils.seismo_logs
   submodules/utils.sfile_util
   submodules/utils.stacking
   submodules/utils.synth_seis
   submodules/utils.trigger
   submodules/utils.catalog_to_dd
   submodules/utils.sac_util
