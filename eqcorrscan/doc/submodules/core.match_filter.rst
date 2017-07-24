match_filter
------------

.. currentmodule:: eqcorrscan.core.match_filter
.. automodule:: eqcorrscan.core.match_filter

    .. comment to end block

    Classes
    -------
    .. toctree::
        :maxdepth: 1

        core.match_filter.Detection
        core.match_filter.Family
        core.match_filter.Party
        core.match_filter.Template
        core.match_filter.Tribe

    Functions
    ---------
    .. autosummary::
       :toctree: autogen
       :nosignatures:

       extract_from_stream
       get_catalog
       match_filter
       normxcorr2
       read_detections
       read_tribe
       read_party
       read_template
       write_catalog

    .. comment to end block

    Private Functions
    -----------------
    .. autosummary::
       :toctree: autogen
       :nosignatures:

       _group_process
       _group_detect
       _write_family
       _read_family
       _total_microsec
       _test_event_similarity
