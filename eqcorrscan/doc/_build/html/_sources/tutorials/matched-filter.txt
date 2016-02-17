Matched-filter detection
========================

In this section we will outline using the templates generated in the first tutorial
to scan for similar earthquakes within a day of data.  This small example does not truly exploit the parallel
operations within this package however, so you would be encouraged to think
about where parallel operations occur (*hint, the code can run one template
per cpu*), and why there are --instance and--splits flags in the other
scripts in the github repository (*hint, if you have heaps of memory
and cpus you can do some brute force day parallelisation!*).

The main processing flow is outlined in the figure below, note the main speedups
in this process are achieved by running multiple templates at once, however this
increases memory usage.  If memory is a problem there are flags (mem_issue) in the
match_filter.py source that can be turned on - the codes will then write temporary
files, which is slower, but can allow for more data crunching at once, your trade-off,
your call.

.. image:: processing_flow.png
     :width: 600px
     :align: center
     :alt: processing_flow.png

.. literalinclude:: ../../tutorials/match_filter.py

References
----------
* CM Boese, J Townend, E Smith, T Stern (2012). `Microseismicity and stress in the vicinity of the Alpine Fault, central Southern Alps, New Zealand <http://onlinelibrary.wiley.com/doi/10.1029/2011JB008460/full>`_, *JGR*, doi:10.1029/2011JB008460
* CM Boese, KM Jacobs, EGC Smith, TA Stern, J Townend (2014). `Background and delayed-triggered swarms in the central Southern Alps, South Island, New Zealand <http://onlinelibrary.wiley.com/doi/10.1002/2013GC005171/full>`_, *G-cubed*, doi:10.1002/2013GC005171
