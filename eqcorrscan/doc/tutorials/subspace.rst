Subspace Detection
==================

EQcorrscan's subspace detection methods are closely modelled on the method
described by |Harris2006|, Subspace Detectors: Theory.  We offer options to
multiplex data or leave as single-channels (multiplexing is significantly
faster).

Subspace detection is implemented in an object-oriented style, whereby individual
detectors are constructed from data, then used to detect within continuous data.
At the core of the subspace routine is a Cython-based, static-typed routine to
calculate the detection statistics.  We do this to make use of numpy's vectorized
calculations, while taking advantage of the speed-ups afforded by compiling
the sliding window loop.

WARNING
-------

Subspace in EQcorrscan is in **beta**, you must check your results match what
you expect - if you find errors please report them.  Although our test-cases run
correctly, changes in data quality may affect the routines in ways we have not
accounted for.

Important
---------

How you generate you detector is likely to be the most important thing, careful
selection and alignment is key, and because of this we haven't provided a total
*cookie-cutter* system for doing this.  You have freedom to chose your parameters,
how to process, how you align, what traces to keep, whether you multiplex or not,
etc.  This also means you have a lot of freedom to **get it wrong**. You will
have to do significant testing with your own dataset to work out what works and
what doesn't.  Anything that you find that doesn't work well in EQcorrscans
system, it would be great to hear about so that we can make it better.

The following examples demonstrate some of the options, but not all of them.
The advanced example is the example used to test and develop subspace and took
a fair amount of effort over a number of weeks.

Simple example
--------------

To begin with you will need to create a **Detector**:

.. code-block:: python

    from eqcorrscan.core import subspace
    detector = subspace.Detector()

This will create an empty *detector* object.  These objects have various attributes,
including the data to be used as a detector (*detector.data*), alongside the full
input and output basis vector matrices (*detector.u* and *detector.v* respectively)
and the vector of singular-values (*detector.sigma*).  Meta-data are also included,
including whether the detector is multiplexed or not (*detector.multiplex*), the
filters applied (*detector.lowcut*, *detector.highcut*, *detection.filt_order*,
*detector.sampling_rate*), the dimension of the subspace (*detector.dimension*),
and the name of the detector, which you can use for book-keeping
(*detector.name*).

To populate the empty detector you need a design set of streams that have been
aligned (see clustering submodule for alignment methods).

.. code-block:: python

    detector.construct(streams=streams, lowcut=2, highcut=9, filt_order=4,
                       sampling_rate=20, multiplex=True, name='Test_1',
                       align=True, shift_len=0.5)

This will populate all the attributes of your *detector* object, and fill the
*detector.data* with the full input basis vector matrix.

You will want to reduce the dimensions of your subspace detector, such that
you are just describing the signal, preferably with a lot of generality.  Details
for selecting dimensionality should be found in |Harris2006|.  To do this in
EQcorrscan simply use the *partition* method:

.. code-block:: python

    detector.partition(4)

This will populate *detector.data* with the first four, left-most input basis
vectors.  You can test to see how much of your original design set is
described by this detector by using the *energy_capture* method:

.. code-block:: python

    percent_capture = detector.energy_capture()

This will return a percentage capture, you can run this for multiple dimensions
to test what dimension best suits your application.  Again, details for this
selection can be found in |Harris2006|.

Finally, to use your detector to detect within continuous data you should use
the *detect* method.  This requires a stream with the same stations and channels
used in the detector, and a threshold from 0-1, where 0 is no signal, and 1 is
totally described by your detector.  You can extract streams for the detections
at the same time as the detections by setting the extract_detections flag to
True.

.. code-block:: python

    detections = detector.detect(st=stream, threshold=0.5, trig_int=3)


Advanced Example
----------------

This example computes detections for a short data-period during an earthquake
sequence in the Wairarapa region of New Zealand's North Island.  This example only
shows one subspace detector, but could be extended, using the various :doc:`clustering <../submodules/utils.clustering>`
routines in EQcorrscan, to create many subspace detectors.  These could be run
using the :doc:`subspace_detect <../submodules/autogen/eqcorrscan.core.subspace.subspace_detect>`
function, which runs similar
detectors in parallel through the given data.

.. literalinclude:: ../../tutorials/subspace.py

.. |Harris2006| raw:: html

    <a href="https://e-reports-ext.llnl.gov/pdf/335299.pdf" target="_blank">Harris (2006)</a>