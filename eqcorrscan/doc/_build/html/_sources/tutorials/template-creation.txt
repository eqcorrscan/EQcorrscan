Template creation
=================

In this example we will download some data and picks from the New Zealand GeoNet
database and make use of the functions in EQcorrscan to quickly and simply
generate templates for use in matched-filter detection.  In the example we
are looking at an earthquake sequence on the east coast of New Zealand's North
Island that occurred on the 4th of January 2016.  We will take a set of four
template events from the sequence that have been picked by GeoNet, of a range
of magnitudes.  You can decide if these were *good* templates or not.  You could
easily extend this by choosing more template events (the mainshock in the
sequence is a M 5 and you can get the information by clicking |GeoNet_link|).

.. |GeoNet_link| raw:: html

  <a href="http://www.geonet.org.nz/quakes/region/newzealand/2016p008122" target="_blank">here</a>

You do not need to use data from an online server, many pick formats can be
parsed, either by obspy, or (for seisan pick files) through the Sfile_utils
module in EQcorrscan.

This template script is written to be general, so you should be able to give
command line arguments to the script to generate templates from other
FDSN databases.  Note that some data-centers do not support full FDSN quakeml
files, and working out which do is quite painful.

Try this example for another, Northern California Data Center earthquake:

``python template_creation.py NCEDC 72572665``

This will (unfortunately for you) generate a warning about un-labelled picks,
as many data-centers do not provide phase-hints with their picks.  We care about
which phase is which when we have to run our cross-correlation re-picker
(yet to be completed), which means that we, by convention, assign P-picks
to the vertical channel and S-picks to both horizontal channels.

**Note:** To run this script and for all map plotting you will need to install
matplotlib-toolkits basemap package.  Install instructions and a link to the
download are |basemap_link|.

.. |basemap_link| raw:: html

  <a href="http://matplotlib.org/basemap/users/installing.html" target="_blank">here</a>

Important considerations
------------------------
In this tutorial we enforce downloading of day-long data for the template
generation.  This is to ensure that the data we make the template from, and
the data we use for detection are processed in exactly the same way.  If we
were to only download a short segment of data around the event and process this
we would find that the resampling process would result in minor differences
between the templates and the continuous data.  This has the effect that, for
self-detections, the cross-correlation values are less than 1.

This is an important effect and something that you should consider when generating
your own templates.  You **MUST** process your templates in the exact same way
(using the same routines, same filters, same resampling, and same data length)
as your continuous data.  It can have a very significant impact to your results.

The functions provided in eqcorrscan.core.template_gen are there to aid you,
but if you look at the source code, all they are doing is:

* Detrending;
* Resampling;
* Filtering;
* and cutting.

If you want to do these things another way you are more then welcome to!


.. literalinclude:: ../../tutorials/template_creation.py
