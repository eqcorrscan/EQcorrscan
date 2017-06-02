Template creation
=================

Simple example
--------------

Within :doc:`template_gen <../submodules/core.template_gen>` there are lots
of methods for generating templates depending on the type of pick data and
waveform data you are using.  Have a look at those to check their simple
examples.

Example of how to generate a useful template from data available via
FDSN (see |obspy_fdsn| for a list of possible clients):

.. |obspy_fdsn| raw:: html

   <a href="http://docs.obspy.org/packages/obspy.clients.fdsn.html#module-obspy.clients.fdsn" target="_blank">obspy.clients.fdsn</a>

.. code-block:: python

    from obspy.clients.fdsn import Client
    from obspy.core.event import Catalog
    from eqcorrscan.core.template_gen import from_client
    client = Client('NCEDC')
    catalog = client.get_events(eventid='72572665', includearrivals=True)
    templates = from_client(catalog=catalog, client_id='NCEDC',
                            lowcut=2.0, highcut=9.0, samp_rate=20.0,
                            filt_order=4, length=3.0, prepick=0.15,
                            swin='all', process_len=200)

This will download data for a single event (given by eventid) from the NCEDC
database, then use that information to download relevant waveform data.  These
data will then be filtered and cut according the parameters set.

It is often wise to set these parameters once as variables at the start of your
scripts to ensure that the same parameters are used for both template creation and
matched-filtering.

The method *from_client* can cope with multiple events (hence the use of a
Catalog object), so you can download a lot of events from your chosen client and
generate templates for them all.

Useful considerations
---------------------

With these data-mining techniques, memory consumption is often an issue, as well
as speed.  To reduce memory consumption and increase efficiency it is often
worth using a subset of picks.  The advanced example below shows one way to do
this.  In practice, five picks (and therefore traces in a template) is often
sufficient for matched-filter detections.  However, you should test this on your
own data.

Some other things that you might want to consider when generating templates
include:

* Template-length, you probably only want to include the real earthquake signal in your template,
  so really long templates are probably not the best idea.
* On the same note, don't include much (if any) data before the P-phase, unless you have
  good reason to - assuming your noise is random, including noise will reduce the
  correlations.
* Consider your frequency band - look for peak power in the chosen waveform
  **relative to the noise**.
* Coda waves often describe scatterers - scattered waves are very interesting,
  but may reduce the generality of your templates.  If this is what you want, include
  coda, if you want a more general template, I would suggest not including coda.
  For examples of this you could try generating a lot of templates from a sequence
  and computing the SVD of the templates to see where the most coherent energy is
  (in the first basis vector), or just computing the stack of the waveforms.

Storing templates
-----------------

Templates are returned as obspy Stream objects.  You will likely want to store
these templates on disk.  This is usually best done by saving them as miniseed
files.  Miniseed uses an efficient compression algorithm and allows multiplexing,
which is useful for template storage.  However we do not constrain you to this.

.. code-block:: python

    template.write('template.ms', format="MSEED")


Advanced example
----------------

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

This will also show template waveforms for both the automatic and the manual
picks - you can change this by removing either automatic or manual picks
y setting the *evaluation_mode* flag in :func:`eqcorrscan.utils.catalog_utils.filter_picks`.

**Note:** To run this script and for all map plotting you will need to install
matplotlib-toolkits basemap package.  Install instructions and a link to the
download are |basemap_link|. We recommend that you install the package using
your system package manager if it is available.

.. |basemap_link| raw:: html

  <a href="http://matplotlib.org/basemap/users/installing.html" target="_blank">here</a>

**Important considerations**

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
