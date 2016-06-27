Magnitude calculation
=====================
EQcorrscan contains both an automatic amplitude picker and a
singular-value decomposition derived magnitude calculation, which
is very accurate but requires high levels of event similarity.

Amplitude picker for local magnitudes
-------------------------------------

Currently this is only implemented for seisan s-files, however we have no plans
to extend this to other formats (although it would be simple as seisan files
are read in as obspy events not).

This example requires data downloaded from the eqcorrscan github repository.

.. code-block:: python

    from eqcorrscan.utils.mag_calc import amp_pick_sfile
    from obspy.core.event import Event
    import os
    testing_path = 'eqcorrscan/tests/test_data'
    sfile = os.path.join(testing_path, 'REA', 'TEST_',
                         '01-0411-15L.S201309')
    datapath = os.path.join(testing_path, 'WAV', 'TEST_')
    respdir = testing_path
    event = amp_pick_sfile(sfile=sfile, datapath=datapath,
                           respdir=respdir, chans=['Z'], var_wintype=True,
                           winlen=0.9, pre_pick=0.2, pre_filt=True,
                           lowcut=1.0, highcut=20.0, corners=4)

Relative moment by singular-value decomposition
-----------------------------------------------

This method closely follows the method outlined by |SVD_mag_link|.

.. |SVD_mag_link| raw:: html

    <a href="http://www.bssaonline.org/content/100/5A/1952.short" target="_blank">Rubinstein & Ellsworth 2010</a>

This example requires data downloaded from the eqcorrscan github repository.

.. code-block:: python

    from eqcorrscan.utils.mag_calc import SVD_moments
    from obspy import read
    import glob
    import os
    from eqcorrscan.utils.clustering import SVD
    import numpy as np
    # Do the set-up
    testing_path = 'eqcorrscan/tests/test_data/similar_events'
    stream_files = glob.glob(os.path.join(testing_path, '*'))
    stream_list = [read(stream_file) for stream_file in stream_files]
    event_list = []
    for i, stream in enumerate(stream_list):
        st_list = []
        for tr in stream:
            # Only use the vertical channels of sites with known high similarity.
            # You do not need to use this step for your data.
            if (tr.stats.station, tr.stats.channel) not in\
                    [('WHAT2', 'SH1'), ('WV04', 'SHZ'), ('GCSZ', 'EHZ')]:
                stream.remove(tr)
                continue
            tr.detrend('simple')
            tr.filter('bandpass', freqmin=5.0, freqmax=15.0)
            tr.trim(tr.stats.starttime + 40, tr.stats.endtime - 45)
            st_list.append(i)
        event_list.append(st_list)
    event_list = np.asarray(event_list).T.tolist()
    SVectors, SValues, Uvectors, stachans = SVD(stream_list=stream_list)
    M, events_out = SVD_moments(U=Uvectors, s=SValues, V=SVectors,
                                stachans=stachans, event_list=event_list)