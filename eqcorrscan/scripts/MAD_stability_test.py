"""
Simple functions to assert the stability of MAD thresholding at short-time
intervals.  MAD relies on the assumption that large outliers have little
affect on otherwise normally distributed data.  Previous versions of EQcorrscan
have forced the use of day-long data, which has resulted in high memory
consumption, however there is no theoretical need for this.  It would be
desirable to use shorter chunks of data to be more memory efficient,
this would allow for even more parallel processing, and, potentially, allow
for data to be processing in near real-time.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

def test_stability():
    """Test various threshold window lengths."""
    from eqcorrscan.core.match_filter import _channel_loop
    from eqcorrscan.utils import pre_processing, catalog_utils, plotting
    from eqcorrscan.core import template_gen
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime, Trace
    import numpy as np
    import matplotlib.pyplot as plt

    # Do some set-up
    client = Client('NCEDC')
    t1 = UTCDateTime(2004, 9, 28)
    t2 = t1 + 86400
    catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=2,
                                minlatitude=35.7, maxlatitude=36.1,
                                minlongitude=-120.6, maxlongitude=-120.2,
                                includearrivals=True)
    catalog = catalog_utils.filter_picks(catalog, channels=['EHZ'],
                                         top_n_picks=5)
    templates = template_gen.from_client(catalog=catalog, client_id='NCEDC',
                                         lowcut=2.0, highcut=9.0,
                                         samp_rate=20.0, filt_order=4,
                                         length=3.0, prepick=0.15,
                                         swin='all')
    bulk_info = [(tr.stats.network, tr.stats.station, '*',
                  tr.stats.channel[0] + 'H' + tr.stats.channel[1],
                  t2 - 3600, t2) for tr in templates[0]]
    st = client.get_waveforms_bulk(bulk_info)
    st.merge(fill_value='interpolate')
    st = pre_processing.shortproc(st, lowcut=2.0, highcut=9.0,
                                  filt_order=4, samp_rate=20.0,
                                  debug=0, num_cores=4)
    i = 0
    cccsums, no_chans, chans = _channel_loop(templates, st)

    cccsum = cccsums[0]
    MAD_thresh = 8
    MAD_daylong = MAD_thresh * np.median(np.abs(cccsum))
    MAD_hours = []
    MAD_five_mins = []
    MAD_minutes = []
    for hour in range(24):
        ccc_hour_slice = cccsum[hour * 3600 * st[0].stats.sampling_rate:
                                (hour + 1) * 3600 *
                                st[0].stats.sampling_rate]
        MAD_hour_slice = MAD_thresh * np.median(np.abs(ccc_hour_slice))
        MAD_hours.append(MAD_hour_slice)
        for five_min in range(12):
            ccc_five_slice = ccc_hour_slice[five_min * 300 *
                                            st[0].stats.sampling_rate:
                                            (five_min + 1) * 300 *
                                            st[0].stats.sampling_rate]
            MAD_five_slice = MAD_thresh * np.median(np.abs(ccc_five_slice))
            MAD_five_mins.append(MAD_five_slice)
        for minute in range(60):
            ccc_min_slice = ccc_hour_slice[minute * 60 *
                                            st[0].stats.sampling_rate:
                                            (minute + 1) * 60 *
                                            st[0].stats.sampling_rate]
            MAD_min_slice = MAD_thresh * np.median(np.abs(ccc_min_slice))
            MAD_minutes.append(MAD_min_slice)
    plotting_cccsum = Trace(cccsum)
    plotting_cccsum.stats.sampling_rate = st[0].stats.sampling_rate
    plotting_cccsum = plotting.chunk_data(plotting_cccsum, 1, 'Maxabs')
    x = np.arange(0, 24, 1.0 / (3600 * plotting_cccsum.stats.sampling_rate))
    x = x[0:len(plotting_cccsum.data)]
    plt.plot(x, plotting_cccsum.data, linewidth=0.7, color='k')
    plt.plot(np.arange(0, 24, 1.0/60), MAD_minutes, label='1 minute MAD',
             color='y')
    plt.plot(np.arange(0, 24, 1.0/60), [-1 * m for m in MAD_minutes],
             color='y')
    plt.plot(np.arange(0, 24, 1.0/12), MAD_five_mins, label='5 minute MAD',
             color='b')
    plt.plot(np.arange(0, 24, 1.0/12), [-1 * m for m in MAD_five_mins],
             color='b')
    plt.plot(np.arange(24), MAD_hours, label='Hourly MAD',
             color='r', linewidth=1.4)
    plt.plot(np.arange(24), [-1 * m for m in MAD_hours],
             color='r', linewidth=1.4)
    plt.plot([0, 24], [MAD_daylong, MAD_daylong], label='Day-long MAD',
             linewidth=1.5, color='g')
    plt.plot([0, 24], [-1 * MAD_daylong, -1 * MAD_daylong],
             linewidth=1.5, color='g')
    plt.legend()
    plt.show()
