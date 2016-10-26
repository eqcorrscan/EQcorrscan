"""
Advanced subspace tutorial to show some of the capabilities of the method.

This example uses waveforms from a known earthquake sequence (in the Wairarapa
region north of Wellington, New Zealand). The catalogue locations etc can
be downloaded from this link:

http://quakesearch.geonet.org.nz/services/1.0.0/csv?bbox=175.37956,-40.97912,175.53097,-40.84628&startdate=2015-7-18T2:00:00&enddate=2016-7-18T3:00:00

"""


def run_tutorial(plot=False, multiplex=True, return_streams=False):
    """
    Run the tutorial.

    :return: detections
    """
    # We are going to use data from the GeoNet (New Zealand) catalogue. GeoNet
    # do not implement the full FDSN system yet, so we have a hack to get
    # around this.  It is not strictly part of EQcorrscan, so we haven't
    # included it here, but you can find it in the tutorials directory of the
    # github repository
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client
    from eqcorrscan.tutorials.get_geonet_events import get_geonet_events
    from obspy import UTCDateTime
    from eqcorrscan.utils.catalog_utils import filter_picks
    from eqcorrscan.utils.clustering import space_cluster
    from eqcorrscan.core import subspace

    cat = get_geonet_events(minlat=-40.98, maxlat=-40.85, minlon=175.4,
                            maxlon=175.5, startdate=UTCDateTime(2016, 5, 1),
                            enddate=UTCDateTime(2016, 5, 20))
    # This gives us a catalog of events - it takes a while to download all
    # the information, so give it a bit!
    # We will generate a five station, multi-channel detector.
    cat = filter_picks(catalog=cat, top_n_picks=5)
    stachans = list(set([(pick.waveform_id.station_code,
                          pick.waveform_id.channel_code) for event in cat
                         for pick in event.picks]))
    # In this tutorial we will only work on one cluster, defined spatially.
    # You can work on multiple clusters, or try to whole set.
    clusters = space_cluster(catalog=cat, d_thresh=2, show=False)
    # We will work on the largest cluster
    cluster = sorted(clusters, key=lambda c: len(c))[-1]
    # This cluster contains 32 events, we will now download a trim the
    # waveforms.  Note that each chanel must start at the same time and be the
    # same length for multiplexing.  If not multiplexing EQcorrscan will
    # maintain the individual differences in time between channels and delay
    # the detection statistics by that amount before stacking and detection.
    client = Client('GEONET')
    design_set = []
    bulk_info = []
    for event in cluster:
        t1 = event.origins[0].time
        t2 = t1 + 25
        for station, channel in stachans:
            bulk_info.append(('NZ', station, '*', channel[0:2] + '?', t1, t2))
    st = client.get_waveforms_bulk(bulk=bulk_info)
    for event in cluster:
        t1 = event.origins[0].time
        t2 = t1 + 25
        design_set.append(st.copy().trim(t1, t2))
    # Construction of the detector will process the traces, then align them,
    # before multiplexing.
    detector = subspace.Detector()
    detector.construct(streams=design_set, lowcut=2.0, highcut=9.0,
                       filt_order=4, sampling_rate=20, multiplex=multiplex,
                       name='Wairarapa1', align=True, reject=0.2,
                       shift_len=6, plot=plot).partition(9)
    if plot:
        detector.plot()
    # We also want the continuous stream to detect in.
    t1 = UTCDateTime(2016, 5, 11, 19)
    t2 = UTCDateTime(2016, 5, 11, 20)
    # We are going to look in a single hour just to minimize cost, but you can \
    # run for much longer.
    bulk_info = [('NZ', stachan[0], '*',
                  stachan[1][0] + '?' + stachan[1][-1],
                  t1, t2) for stachan in detector.stachans]
    st = client.get_waveforms_bulk(bulk_info)
    st.merge().detrend('simple').trim(starttime=t1, endtime=t2)
    # We set a very low threshold because the detector is not that great, we
    # haven't aligned it particularly well - however, at this threshold we make
    # two real detections.
    detections, det_streams = detector.detect(st=st, threshold=0.005,
                                              trig_int=2,
                                              extract_detections=True)
    if return_streams:
        return detections, det_streams
    else:
        return detections


if __name__ == '__main__':
    run_tutorial()
