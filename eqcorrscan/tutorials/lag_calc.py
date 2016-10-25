"""Tutorial to illustrate the lag_calc usage."""


def run_tutorial(min_magnitude=2, shift_len=0.2, num_cores=4, min_cc=0.5):
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client
    from obspy.core.event import Catalog
    from obspy import UTCDateTime
    from eqcorrscan.core import template_gen, match_filter, lag_calc
    from eqcorrscan.utils import pre_processing, catalog_utils

    client = Client('NCEDC')
    t1 = UTCDateTime(2004, 9, 28)
    t2 = t1 + 86400
    print('Downloading catalog')
    catalog = client.get_events(starttime=t1, endtime=t2,
                                minmagnitude=min_magnitude,
                                minlatitude=35.7, maxlatitude=36.1,
                                minlongitude=-120.6, maxlongitude=-120.2,
                                includearrivals=True)
    # We don't need all the picks, lets take the information from the
    # five most used stations - note that this is done to reduce computational
    # costs.
    catalog = catalog_utils.filter_picks(catalog, channels=['EHZ'],
                                         top_n_picks=5)
    # There is a duplicate pick in event 3 in the catalog - this has the effect
    # of reducing our detections - check it yourself.
    for pick in catalog[3].picks:
        if pick.waveform_id.station_code == 'PHOB' and pick.onset == 'emergent':
            catalog[3].picks.remove(pick)
    print('Generating templates')
    templates = template_gen.from_client(catalog=catalog, client_id='NCEDC',
                                         lowcut=2.0, highcut=9.0,
                                         samp_rate=50.0, filt_order=4,
                                         length=3.0, prepick=0.15,
                                         swin='all', process_len=3600)
    # In this section we generate a series of chunks of data.
    start_time = UTCDateTime(2004, 9, 28, 17)
    end_time = UTCDateTime(2004, 9, 28, 20)
    process_len = 3600
    chunks = []
    chunk_start = start_time
    while chunk_start < end_time:
        chunk_end = chunk_start + process_len
        if chunk_end > end_time:
            chunk_end = end_time
        chunks.append((chunk_start, chunk_end))
        chunk_start += process_len

    all_detections = []
    picked_catalog = Catalog()
    template_names = [str(template[0].stats.starttime)
                      for template in templates]
    for t1, t2 in chunks:
        print('Downloading and processing for start-time: %s' % t1)
        # Download and process the data
        bulk_info = [(tr.stats.network, tr.stats.station, '*',
                      tr.stats.channel[0] + 'H' + tr.stats.channel[1],
                      t1, t2) for tr in templates[0]]
        # Just downloading a chunk of data
        st = client.get_waveforms_bulk(bulk_info)
        st.merge(fill_value='interpolate')
        st = pre_processing.shortproc(st, lowcut=2.0, highcut=9.0,
                                      filt_order=4, samp_rate=50.0,
                                      debug=0, num_cores=num_cores)
        detections = match_filter.match_filter(template_names=template_names,
                                               template_list=templates,
                                               st=st, threshold=8.0,
                                               threshold_type='MAD',
                                               trig_int=6.0, plotvar=False,
                                               plotdir='.', cores=num_cores)
        # Extract unique detections from set.
        unique_detections = []
        for master in detections:
            keep = True
            for slave in detections:
                if not master == slave and\
                   abs(master.detect_time - slave.detect_time) <= 1.0:
                    # If the events are within 1s of each other then test which
                    # was the 'best' match, strongest detection
                    if not master.detect_val > slave.detect_val:
                        keep = False
                        break
            if keep:
                unique_detections.append(master)
        all_detections += unique_detections

        picked_catalog += lag_calc.lag_calc(detections=unique_detections,
                                            detect_data=st,
                                            template_names=template_names,
                                            templates=templates,
                                            shift_len=shift_len, min_cc=min_cc,
                                            interpolate=False, plot=False,
                                            parallel=True, debug=3)
    # Return all of this so that we can use this function for testing.
    return all_detections, picked_catalog, templates, template_names

if __name__ == '__main__':
    from multiprocessing import cpu_count
    run_tutorial(min_magnitude=4, num_cores=cpu_count())
