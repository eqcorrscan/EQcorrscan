"""
Functions for common operations on whole obspy.Catalog classes

Copyright 2015 Calum Chamberlain

This file is part of EQcorrscan.

    EQcorrscan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EQcorrscan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EQcorrscan.  If not, see <http://www.gnu.org/licenses/>.

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import warnings
import seaborn as sns
import matplotlib.pyplot as plt


def sta_res_dict(catalog):
    r"""
    Function to extract all travel time residuals and weights and sort them\
    into a dictionary keyed by station name.
    """
    sta_dict = {}
    for event in catalog:
        origin = event.preferred_origin() or event.origins[0]
        for arrival in origin.arrivals:
            assoc_pick = arrival.pick_id.get_referred_object()
            pick_sta = assoc_pick.waveform_id.station_code
            if pick_sta not in sta_dict:
                sta_dict[pick_sta] = {'res': [arrival.time_residual],
                                      'wt': [arrival.time_weight]}
            else:
                sta_dict[pick_sta]['res'].append(arrival.time_residual)
                sta_dict[pick_sta]['wt'].append(arrival.time_weight)
    return sta_dict


def extract_all_residuals(catalog):
    r"""
    Function to create a list of all travel time residuals in a catalog.
    """
    for event in catalog:
        origin = event.preferred_origin() or event.origins[0]
        arr_time_resds = [abs(x.time_residual) for x in origin.arrivals]
    return arr_time_resds


def jeffreys_weighting(catalog, plotvar=True):
    r"""
    Function to perform Jeffreys Weighting on a set of picks per the\
    methodology detailed here:

    Jeffreys, H (1973), On travel times in seismology, in Collected Papers\
        of Sir Harold Jeffreys on Geophysics and Other Sciences, edited by\
        H Jeffreys, pp 36 120, Gordon and Breach, London

    Steps are as follows:
        1. Remove all picks with residuals greater than 2 sigma
        2. For the resulting catalog, again remove all picks with residuals\
            greater than 2 sigma
        3. Apply weighting to each pick as per eqn 3 pp 39 of the above\
            citation

    :type catalog: :class: obspy.Catalog
    :param catalog: Catalog for which to assign pick weights
    :type plotvar: bool
    :param plotvar: Whether or not to plot the pick weights for each station.

    :returns: :class: obspy.Catalog with pick weights added
    """
    # First outlier removal
    cat1 = filter_picks(catalog, sigma=2)
    # Second outlier removal
    cat2 = filter_picks(cat1, sigma=2)
    # Calculate std dev and m for remaining residuals
    final_resds = extract_all_residuals(cat2)
    res_std = np.std(final_resds)
    # Quartiles for residuals
    q75, q25 = np.percentile(final_resds, [75, 25])
    print('Lower quartile: %s' % str(q25))
    print('Upper quartile: %s' % str(q75))
    # Ratio of outliers to all arrivals
    m = len([x for x in final_resds if x < q25 or x > q75]) / len(final_resds)
    print('m equals: %s' % str(m))
    # Loop over all arrivals and assign weights
    for event in cat2:
        origin = event.preferred_origin() or event.origins[0]
        for arr in origin.arrivals:
            wt = 1 / (1 + m * np.exp(arr.time_residual ** 2 / res_std ** 2))
            print('weight equals: %s' % str(wt))
            # XXX Putting this here because I don't know where else ATM
            arr.time_weight = wt
    return catalog


def filter_picks(catalog, sigma):
    r"""
    Function to remove single picks with travel time residuals outside of\
    an integer multiple of standard deviation. Works in-place on the catalog\
    so the original catalog cannot be referenced later.

    :type catalog: :class: obspy.Catalog
    :param catalog: catalog from which to remove picks
    :type sigma: int
    :param sigma: Number of standard of deviations outside of which picks are\
        removed

    :returns: :class: obspy.Catalog
    """
    new_cat = catalog.copy()
    residuals = extract_all_residuals(new_cat)
    mean_res = np.mean(residuals)
    res_std = np.std(residuals)
    for event in new_cat:
        origin = event.preferred_origin() or event.origins[0]
        for arrival in origin.arrivals:
            if arrival.time_residual < mean_res - (res_std * sigma) or\
               arrival.time_residual > mean_res + (res_std * sigma):
                pick = arrival.pick_id.get_referred_object()
                event.picks.remove(pick)
            else:
                continue
    return new_cat


def filter_events(catalog, method='avg_residual', plot=False):
    r"""
    Function to remove events with unsatisfactory picks from a Catalog

    :type catalog: :class: 'obspy.Catalog'
    :param catalog: Catalog from which to remove events
    :type method: str
    :param method: Method used to determine which events to remove from\
        catalog. Options are 'avg_residual' and 'single_arrival'. Defaults to\
        'avg_residual' which will remove events with an average arrival time\
        residual outside of one standard deviation accross the catalog.\
        'single_arrival' will remove events which contain one or more arrivals\
        with time residuals outside of one standard deviation for all arrival\
        time residuals in the entire catalog.
    :type plot: bool
    :param plot: If True, will plot the distribution of either average event\
        arrival time residual or all arrival time residuals for the catalog,\
        depending upon which method is used.

    :returns: class: obspy.Catalog
    """
    from obspy import Catalog

    # Extract all arrivals for each preferred origin
    arr_time_resds = extract_all_residuals(catalog)

    # Calculate average arrival time residual for all preferred origins
    avg_arr_res = []
    for event in catalog:
        pref_o = event.preferred_origin() or event.origins[0]
        # Calculate average arrival time residual for origin
        avg_arr_res.append(sum([abs(x.time_residual) for
                           i, x in enumerate(pref_o.arrivals)]) / i)
    # Plot the histograms
    if plot:
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
        ax1.set_title('Arrival-time residuals')
        sns.boxplot(data=arr_time_resds, ax=ax1, orient="h", width=0.1)
        sns.distplot(arr_time_resds, ax=ax2)
        ax3.set_title('Event average arrival-time residuals')
        sns.boxplot(data=avg_arr_res, ax=ax3, orient="h", width=0.1)
        sns.distplot(avg_arr_res, ax=ax4)
        plt.show()
        plt.close()
    # Creat new, filtered catalog
    filtered_cat = Catalog()
    if method == 'avg_residual':
        mean_avg = np.mean(avg_arr_res)
        std_avg = np.std(avg_arr_res)
        for event in catalog:
            pref_o = event.preferred_origin() or event.origins[0]
            avg_arr_res = sum([x.time_residual for
                               x in pref_o.arrivals]) / len(pref_o.arrivals)
            if avg_arr_res < mean_avg + std_avg and\
                    avg_arr_res > mean_avg - std_avg:
                filtered_cat.append(event)
            else:
                continue
    elif method == 'single_arrival':
        mean_res = np.mean(arr_time_resds)
        std_res = np.std(arr_time_resds)
        for event in catalog:
            pref_o = event.preferred_origin() or event.origins[0]
            bad_arrivals = [x for x in pref_o.arrivals
                            if x.time_residual < mean_res - std_res or
                            x.time_residual > mean_res + std_res]
            if bad_arrivals:
                del bad_arrivals
                continue
            else:
                filtered_cat.append(event)
    return filtered_cat


def uncommon_picks(event1, event2):
    r"""
    Helper to find picks in event1 for which there are no matches in event2

    There must be a more elegant way to do this!!

    :returns: List of obspy.event.Pick
    """
    sta_pha1 = [p.waveform_id for p in event1.picks]
    sta_pha2 = [p.waveform_id for p in event2.picks]
    mismatches = [mm for mm in sta_pha1 if mm not in sta_pha2]
    uncom_picks = [p for p in event1.picks if p.waveform_id in mismatches]
    return uncom_picks


def refine_picks(catalog, stream_dict, pre_pick, post_pick, shift_len,
                 cc_thresh, master=None, lowcut=1.0, highcut=20.0,
                 plotvar=False):
    r"""Function to refine picks in a catalog based upon either a pre-chosen\
    master event or the event in the catalog with the highest amplitude.

    :type catalog: class: obspy.Catalog
    :param catalog: Catalog of events which we want to adjust picks for
    :type stream_dict: dict
    :param stream_dict: Dictionary with key:value pairing of event\
        ResourceID:obspy.Stream for each event in catalog.
    :type pre_pick: float
    :param pre_pick: Time before the pick to start the correlation window
    :type post_pick: float
    :param post_pick: Time after the pick to start the correlation window
    :type shift_len: float
    :param shift_len: Time to allow pick to vary
    :type master: bool or str
    :param master: If 'None', master event defaults to the event with the\
        highest SNR. Otherwise, must specify a valid event resoure_id\
        from the catalog.
    :type lowcut: float
    :param lowcut: Lowcut in Hz - default=1.0
    :type highcut: float
    :param highcut: Highcut in Hz - deafult=10.0

    :returns: class: obspy.Catalog
    """

    import obspy
    if int(obspy.__version__.split('.')[0]) > 0:
        from obspy.signal.cross_correlation import xcorr_pick_correction
    else:
        from obspy.signal.cross_correlation import xcorrPickCorrection \
            as xcorr_pick_correction

    # Establish master template if not specified
    if master:
        master_id = obspy.ResourceIdentifier(master)
    else:
        # Find event with highest SNR to be master
        avg_snr = {}
        for event in catalog:
            avg_snr[event.resource_id] =\
                sum([x.snr for x in event.amplitudes]) / len(event.amplitudes)
        master_id = max(avg_snr.iterkeys(), key=(lambda key: avg_snr[key]))
    # Loop back through catalog and extract master event (there better way?)
    master_event = [x for x in catalog if x.resource_id == master_id][0]
    master_stream = stream_dict[master_id]

    new_catalog = obspy.Catalog()
    # Figure total number of picks
    tot_pks = 0
    for event in catalog:
        for cnt_pick in event:
            tot_pks += 1
    refined_num = 0
    # Now loop the master through all events in catalog
    for slave_event in catalog:
        # Copy old slave event and reset the picks (keep the rest of the info)
        # new_event = obspy.core.event.Event()
        new_event = slave_event.copy()
        new_event.picks = []
        slave_stream = stream_dict[slave_event.resource_id]
        # Find UNcommon picks between slave and master
        mismatches = uncommon_picks(slave_event, master_event)
        # Append them to new event (otherwise they get missed)
        for uncom_pick in mismatches:
            new_event.picks.append(uncom_pick)
        for pick in master_event.picks:
            # Find station, phase pairs
            # Added by Carolin
            slave_matches = [p for p in slave_event.picks
                             if p.phase_hint == pick.phase_hint
                             and p.waveform_id.station_code ==
                             pick.waveform_id.station_code]
            if master_stream.select(station=pick.waveform_id.station_code,
                                    channel='*' +
                                    pick.waveform_id.channel_code[-1]):
                mastertr = master_stream.\
                    select(station=pick.waveform_id.station_code,
                           channel='*' +
                           pick.waveform_id.channel_code[-1])[0]
            else:
                print('No waveform data for ' +
                      pick.waveform_id.station_code + '.' +
                      pick.waveform_id.channel_code)
                break
            for slave_pick in slave_matches:
                if slave_stream.select(station=slave_pick.waveform_id.
                                       station_code,
                                       channel='*'+slave_pick.waveform_id.
                                       channel_code[-1]):
                    slavetr = slave_stream.\
                        select(station=slave_pick.waveform_id.station_code,
                               channel='*'+slave_pick.waveform_id.
                               channel_code[-1])[0]
                else:
                    print('No slave data for ' +
                          slave_pick.waveform_id.station_code + '.' +
                          slave_pick.waveform_id.channel_code)
                    break
                try:
                    correction, cc =\
                        xcorr_pick_correction(pick.time, mastertr,
                                              slave_pick.time,
                                              slavetr, pre_pick, post_pick,
                                              shift_len, filter="bandpass",
                                              filter_options={'freqmin':
                                                              lowcut,
                                                              'freqmax':
                                                              highcut},
                                              plot=plotvar)
                    if abs(correction) > shift_len:
                        warnings.warn('Shift correction too large, ' +
                                      'will not use')
                        new_event.picks.append(slave_pick)
                        continue
                    if cc > cc_thresh:
                        print('Threshold exceeded')
                        new_pick_time = slave_pick.time + correction
                        new_pick = slave_pick.copy()
                        new_pick.time = new_pick_time
                        new_pick.creation_info.agency_id = 'VUW'
                        new_pick.creation_info.author = 'eqcorrscan.refine_picks()'
                        new_pick.creation_info.creation_time = obspy.UTCDateTime.now()
                        new_event.picks.append(new_pick)
                        refined_num += 1
                    else:
                        # new_event.picks.append(slave_pick)
                        print('Correlation not good enough to correct pick')
                        new_event.picks.append(slave_pick)
                except:
                    # Should warn here
                    msg = "Couldn't compute correlation correction"
                    warnings.warn(msg)
                    new_event.picks.append(slave_pick)
                    continue
        new_catalog += new_event
    print('Refined %d of %d picks' % (refined_num, tot_pks))
    return new_catalog
