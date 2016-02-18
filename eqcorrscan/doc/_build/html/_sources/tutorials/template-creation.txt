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
sequence is a M 5 and you can get the information by clicking
`here <http://www.geonet.org.nz/quakes/region/newzealand/2016p008122>`_).

You do not need to use data from an online server, many pick formats can be
parsed, either by obspy, or (for seisan pick files) through the Sfile_utils
module in EQcorrscan.

.. literalinclude:: ../../tutorials/template_creation.py
