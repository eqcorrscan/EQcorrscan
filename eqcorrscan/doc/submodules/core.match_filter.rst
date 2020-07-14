match_filter
============

.. currentmodule:: eqcorrscan.core.match_filter

Classes and functions for matched-filtering.

This is designed for large-scale, multi-paralleled
detection, with moderate to large numbers (hundreds to a few thousands)
of templates.

Object-oriented API
-------------------

EQcorrscan's matched-filter object oriented API is split into two halves, input
(:doc:`Template <core.match_filter.template>` and :doc:`Tribe <core.match_filter.tribe>` objects)
and output (:doc:`Detection <core.match_filter.detection>`,
:doc:`Family <core.match_filter.family>` and :doc:`Party <core.match_filter.party>` objects).

- A :doc:`Tribe <core.match_filter.tribe>` is a collection of :doc:`Templates <core.match_filter.template>`.
- A :doc:`Family <core.match_filter.family>` is a collection of :doc:`Detections <core.match_filter.detection>` for a single :doc:`Template <core.match_filter.template>`.
- A :doc:`Party <core.match_filter.party>` is a collection of :doc:`Families <core.match_filter.family>` for a collection of Templates (a :doc:`Tribe <core.match_filter.tribe>`).

These objects retain useful meta-data including the obspy Event associated with the :doc:`Template <core.match_filter.template>`, and how the
:doc:`Template <core.match_filter.template>` was processed. A core aim of the object-oriented API is reproducibility, and we encourage first-time
users to adopt this from the off.

The methods of :doc:`Template <core.match_filter.template>` and :doc:`Tribe <core.match_filter.tribe>` objects mirror
each other and in general you are best-off working with :doc:`Tribes <core.match_filter.tribe>` .
Both the :doc:`Template <core.match_filter.template>` and :doc:`Tribe <core.match_filter.tribe>` objects
have :meth:`detect <core.match_filter.tribe.Tribe.detect>` methods, which allow you to conduct matched-filtering.

The :doc:`Detection <core.match_filter.detection>`, :doc:`Family <core.match_filter.family>` and :doc:`Party <core.match_filter.party>`
objects contain all the metadata needed to re-create your detections. Furthermore,
the :doc:`Family <core.match_filter.family>` and :doc:`Party <core.match_filter.party>` objects have a
:meth:`lag_calc <core.match_filter.family.Family.lag_calc>` method for conducting cross-correlation phase-picking based on
correlation with the :doc:`Template <core.match_filter.template>`.

=============================================== ===================================================
Object                                          Purpose
=============================================== ===================================================
:doc:`Template <core.match_filter.template>`    To contain the template waveform and meta-data used
                                                to create the template.
:doc:`Tribe <core.match_filter.tribe>`          A collection of multiple Templates. Use the `detect`
                                                method to run matched-filter detections!
:doc:`Detection <core.match_filter.detection>`  Root of the detection object tree - contains
                                                information relevant to a single detection from a
                                                single template.
:doc:`Family <core.match_filter.family>`        Collection of detections for a single Template.
:doc:`Party <core.match_filter.party>`          Collection of Family objects.
=============================================== ===================================================

Function-based API
------------------

core.match_filter.matched_filter
....

.. currentmodule:: eqcorrscan.core.match_filter.matched_filter
.. automodule:: eqcorrscan.core.match_filter.matched_filter

    .. comment to end block

    Functions
    ---------
    .. autosummary::
       :toctree: autogen
       :nosignatures:

       match_filter

    .. comment to end block

core.match_filter.helpers
....

.. currentmodule:: eqcorrscan.core.match_filter.helpers
.. automodule:: eqcorrscan.core.match_filter.helpers

    .. comment to end block

    Functions
    ---------
    .. autosummary::
       :toctree: autogen
       :nosignatures:

       extract_from_stream
       normxcorr2
       temporary_directory


