Magnitude calculation
=====================
EQcorrscan contains both an automatic amplitude picker and a
singular-value decomposition derived magnitude calculation, which
is very accurate but requires high levels of event similarity.

Relative moment by singular-value decomposition
-----------------------------------------------

This method closely follows the method outlined by |SVD_mag_link|.

.. |SVD_mag_link| raw:: html

    <a href="https://pubs.geoscienceworld.org/ssa/bssa/article/100/5A/1952/325099" target="_blank">Rubinstein & Ellsworth 2010</a>

This example requires data downloaded from the eqcorrscan github repository.

.. code-block:: python

    >>> from eqcorrscan.utils.mag_calc import svd_moments
    >>> from obspy import read
    >>> import glob
    >>> from eqcorrscan.utils.clustering import svd
    >>> import numpy as np
    >>> from eqcorrscan import tests
    >>> import os
    >>> # Get the path for the test-data so we can test this
    >>> testing_path = os.path.dirname(
    ...     tests.__file__) + '/test_data/similar_events_processed'
    >>> stream_files = glob.glob(os.path.join(testing_path, '*'))
    >>> stream_list = [read(stream_file) for stream_file in stream_files]
    >>> event_list = []
    >>> for i, stream in enumerate(stream_list): # doctest:+ELLIPSIS
    ...     st_list = []
    ...     for tr in stream:
    ...         # Only use the vertical channels of sites with known high similarity.
    ...         # You do not need to use this step for your data.
    ...         if (tr.stats.station, tr.stats.channel) not in\
    ...                 [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]:
    ...             stream.remove(tr)
    ...             continue
    ...         st_list.append(i)
    ...     event_list.append(st_list)
    <obspy.core.stream.Stream object at ...>
    >>> event_list = np.asarray(event_list).T.tolist()
    >>> SVectors, SValues, Uvectors, stachans = svd(stream_list=stream_list)
    >>> M, events_out = svd_moments(u=Uvectors, s=SValues, v=SVectors,
    ...                             stachans=stachans, event_list=event_list) # doctest:+ELLIPSIS
