Clustering and stacking
=======================

Prior to template generation, it may be beneficial to cluster earthquake
waveforms.  Clusters of earthquakes with similar properties can then be
stacked to create higher signal-to-noise templates that describe the dataset
well (and require fewer templates to reduce computational cost).  Clusters
could also be reduced to their singular-vectors and employed in a subspace
detection routine (coming soon to eqcorrscan).

The following outlines a few examples of clustering and stacking.

Cluster in space
----------------

Download a catalog of global earthquakes and cluster in space, set the distance
threshold to 1,000km

.. code-block:: python

    from eqcorrscan.utils.clustering import space_cluster
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    client = Client("IRIS")
    starttime = UTCDateTime("2002-01-01")
    endtime = UTCDateTime("2002-02-01")
    cat = client.get_events(starttime=starttime, endtime=endtime,
                            minmagnitude=6, catalog="ISC")
    groups = space_cluster(catalog=cat, d_thresh=1000, show=False)

Download a local catalog of earthquakes and cluster much finer (distance
threshold of 2km).

.. code-block:: python

    from eqcorrscan.utils.clustering import space_cluster
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    client = Client("NCEDC")
    starttime = UTCDateTime("2002-01-01")
    endtime = UTCDateTime("2002-02-01")
    cat = client.get_events(starttime=starttime, endtime=endtime,
                            minmagnitude=2)
    groups = space_cluster(catalog=cat, d_thresh=2, show=False)


Setting show to true will plot the dendrogram for grouping with individual
groups plotted in different colours.  The clustering is performed using scipy's
|clust_link|.  Specifically clustering is performed using the linkage method,
which is an agglomorative clustering routine. EQcorrscan uses the average method
with the euclidean distance metric.

.. |clust_link| raw:: html

    <a href="http://docs.scipy.org/doc/scipy-0.16.0/reference/cluster.hierarchy.html" target="_blank">hierachical clustering routines</a>
Cluster in time and space
-------------------------

EQcorrscan's space-time clustering routine first computes groups in space, using
the space_cluster method, then splits the returned groups based on their
inter-event time.

The following example extends the example of the global catalog with a 1,000km
distance threshold and a one-day temporal limit.

.. code-block:: python

    from eqcorrscan.utils.clustering import space_time_cluster
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    client = Client("IRIS")
    starttime = UTCDateTime("2002-01-01")
    endtime = UTCDateTime("2002-02-01")
    cat = client.get_events(starttime=starttime, endtime=endtime,
                            minmagnitude=6, catalog="ISC")
    groups = space_time_cluster(catalog=cat, t_thresh=86400, d_thresh=1000)


Cluster according to cross-correlation values
---------------------------------------------

Waveforms from events are often best grouped based on their similarity.
EQcorrscan has a method to compute clustering based on average cross-correlations.
This again uses scipy's |clust_link|, however in this case clusters are computed
using the *single* method.  Distances are computed from the average of the
multi-chanel cross-correlation values.

The following example uses data stored in the EQcorrscan github repository,
in the tests directory.

.. code-block:: python

    from obspy import read
    import glob
    import os
    from eqcorrscan.utils.clustering import cluster
    # You will need to edit this line to the location of your eqcorrscan repo.
    testing_path = 'eqcorrscan/tests/test_data/similar_events'
    stream_files = glob.glob(os.path.join(testing_path, '*'))
    stream_list = [(read(stream_file), i)
                   for i, stream_file in enumerate(stream_files)]
    for stream in stream_list:
        for tr in stream[0]:
            if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                stream[0].remove(tr)
                continue
            tr.detrend('simple')
            tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
            tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
    groups = cluster(template_list=stream_list, show=False,
                     corr_thresh=0.3)

Stack waveforms (linear)
------------------------

Following from clustering, similar waveforms can be stacked.  EQcorrscan includes
two stacking algorithms, a simple linear stacking method, and a phase-weighted
stacking method.

The following examples use the test data in the eqcorrscan github repository.

.. code-block:: python

    from obspy import read
    import glob
    import os
    from eqcorrscan.utils.clustering import cluster
    from eqcorrscan.utils.stacking import linstack
    # You will need to edit this line to the location of your eqcorrscan repo.
    testing_path = 'eqcorrscan/tests/test_data/similar_events'
    stream_files = glob.glob(os.path.join(testing_path, '*'))
    stream_list = [(read(stream_file), i)
                   for i, stream_file in enumerate(stream_files)]
    for stream in stream_list:
        for tr in stream[0]:
            if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                stream[0].remove(tr)
                continue
            tr.detrend('simple')
            tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
            tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
    groups = cluster(template_list=stream_list, show=False,
                     corr_thresh=0.3)
    # groups[0] should contain 3 streams, which we can now stack
    # Groups are returned as lists of tuples, of the stream and event index
    group_streams = [st_tuple[0] for st_tuple in groups[0]]
    stack = linstack(streams=group_streams)



Stack waveforms (phase-weighted)
--------------------------------

The phase-weighted stack method closely follows the method outlined by
|Thurber_PWS_link|. In this method the linear stack is weighted by the stack
of the instantaneous phase.  In this manor coherent signals are amplified.

.. |Thurber_PWS_link| raw:: html

    <a href="http://www.bssaonline.org/content/early/2014/08/12/0120140077.abstract" target="_blank">Thurber et al. 2014</a>

.. code-block:: python

    from obspy import read
    import glob
    import os
    from eqcorrscan.utils.clustering import cluster
    from eqcorrscan.utils.stacking import PWS_stack
    # You will need to edit this line to the location of your eqcorrscan repo.
    testing_path = 'eqcorrscan/tests/test_data/similar_events'
    stream_files = glob.glob(os.path.join(testing_path, '*'))
    stream_list = [(read(stream_file), i)
                   for i, stream_file in enumerate(stream_files)]
    for stream in stream_list:
        for tr in stream[0]:
            if tr.stats.station not in ['WHAT2', 'WV04', 'GCSZ']:
                stream[0].remove(tr)
                continue
            tr.detrend('simple')
            tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
            tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
    groups = cluster(template_list=stream_list, show=False,
                     corr_thresh=0.3)
    # groups[0] should contain 3 streams, which we can now stack
    # Groups are returned as lists of tuples, of the stream and event index
    group_streams = [st_tuple[0] for st_tuple in groups[0]]
    stack = PWS_stack(streams=group_streams)