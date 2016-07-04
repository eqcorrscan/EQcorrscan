"""Tutorial to illustrate the lag_calc usage."""

def run_tutorial(shift_len=0.2):
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    from eqcorrscan.core import template_gen, match_filter, lag_calc
    from eqcorrscan.utils import pre_processing, catalog_utils

    client = Client('NCEDC')
    t1 = UTCDateTime(2004, 9, 28)
    t2 = t1 + 86400
    catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=4,
                                minlatitude=35.7, maxlatitude=36.1,
                                minlongitude=-120.6, maxlongitude=-120.2,
                                includearrivals=True)
    # We don't need all the picks, lets take the information from the
    # five most used stations - note that this is done to reduce computational
    # costs.
    catalog = catalog_utils.filter_picks(catalog, channels=['EHZ'],
                                         top_n_picks=5)
    templates = template_gen.from_client(catalog=catalog, client_id='NCEDC',
                                         lowcut=2.0, highcut=9.0,
                                         samp_rate=50.0, filt_order=4,
                                         length=3.0, prepick=0.15,
                                         swin='all')
    # Download and process the day-long data
    bulk_info = [(tr.stats.network, tr.stats.station, '*',
                  tr.stats.channel[0] + 'H' + tr.stats.channel[1],
                  t1, t2) for tr in templates[0]]
    st = client.get_waveforms_bulk(bulk_info)
    st.merge(fill_value='interpolate')
    st = pre_processing.dayproc(st, lowcut=2.0, highcut=9.0,
                                filt_order=4, samp_rate=50.0,
                                debug=0, starttime=t1, num_cores=4)
    template_names = [str(template[0].stats.starttime)
                      for template in templates]
    detections = match_filter.match_filter(template_names=template_names,
                                           template_list=templates,
                                           st=st, threshold=8.0,
                                           threshold_type='MAD',
                                           trig_int=6.0, plotvar=False,
                                           plotdir='.', cores=4)
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

    picked_catalog = lag_calc.lag_calc(detections=unique_detections,
                                       detect_data=st,
                                       template_names=template_names,
                                       templates=templates,
                                       shift_len=shift_len, min_cc=0.5,
                                       interpolate=True, plot=False)
    # Return all of this so that we can use this function for testing.
    return unique_detections, picked_catalog, templates, template_names

if __name__ == '__main__':
    run_tutorial()