subspace
--------

.. currentmodule:: eqcorrscan.core.subspace
.. automodule:: eqcorrscan.core.subspace

Subspace detection for either single-channel cases, or
network cases.  This is modelled on the methods described by Harris_.  This method
allows for slightly more variation in detected waveforms than the traditional
matched-filter method.  In this method, templates are constructed either by
using the empirical subspace method, or by computing the basis vectors by
singular-value decomposition.  Both methods are provided as part of EQcorrscan
in the clustering_ module.

.. _Harris: https://e-reports-ext.llnl.gov/pdf/335299.pdf
.. _clustering: utils.clustering.html

    .. comment to end block

    Classes
    -------
    .. toctree::
        :maxdepth: 1

        core.subspace.Detector

    Functions
    ---------
    .. autosummary::
       :toctree: autogen
       :nosignatures:

       read_detector
       multi
       subspace_detect

    .. comment to end block
