"""
A series of test functions for the core functions in EQcorrscan.
"""

from __future__ import division


def test_lag_calc():
    """
    Function to test the capabilites of lag_calc.  When bugs are found in the\
    code and fixed the cause of the bug should be emulated here to ensure that\
    it is truely fixed.
    """
    # We need to develop a synthetic dataset...
    # This needs to have a list of detection objects, a day of synth data,
    # and a series of templates.
    return True


def test_match_filter(samp_rate=20.0, debug=0):
    """
    Function to test the capabilities of match_filter and just check that it\
    is working!  Uses synthetic templates and seeded, randomised data.

    :type debug: int
    :param debug: Debug level, higher the number the more output.
    """
    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils import EQcorrscan_plotting as plotting
    from eqcorrscan.core import match_filter
    from obspy import UTCDateTime
    import string
    # Generate a random dataset
    templates, data, seeds = generate_synth_data(nsta=5, ntemplates=2,
                                                 nseeds=50,
                                                 samp_rate=samp_rate,
                                                 t_length=6.0, max_amp=5.0,
                                                 debug=debug)
    # Notes to the user: If you use more templates you should ensure they are
    # more different, e.g. set the data to have larger moveouts, otherwise
    # similar templates will detect events seeded by another template.
    # Test the pre_processing functions
    for tr in data:
        pre_processing.dayproc(tr=tr, lowcut=2.0, highcut=8.0, filt_order=3,
                               samp_rate=20.0, debug=0,
                               starttime=UTCDateTime(0))
    if debug > 0:
        data.plot()
    # Filter the data and the templates
    for template in templates:
        pre_processing.shortproc(st=template, lowcut=2.0, highcut=8.0,
                                 filt_order=3, samp_rate=20.0)
        if debug > 0:
            template.plot()
    template_names = list(string.ascii_lowercase)[0:len(templates)]
    detections = match_filter.match_filter(template_names=template_names,
                                           template_list=templates,
                                           st=data, threshold=8.0,
                                           threshold_type='MAD',
                                           trig_int=6.0,
                                           plotvar=False,
                                           plotdir='.',
                                           cores=1)
    # Compare the detections to the seeds
    print('This test made ' + str(len(detections)) + ' detections')
    ktrue = 0
    kfalse = 0
    for detection in detections:
        print(detection.template_name)
        i = template_names.index(detection.template_name)
        t_seeds = seeds[i]
        dtime_samples = int((detection.detect_time - UTCDateTime(0)) *
                            samp_rate)
        if dtime_samples in t_seeds['time']:
            j = list(t_seeds['time']).index(dtime_samples)
            print('Detection at SNR of: ' + str(t_seeds['SNR'][j]))
            ktrue += 1
        else:
            min_diff = min(abs(t_seeds['time'] - dtime_samples))
            if min_diff < 10:
                # If there is a match within ten samples then its good enough
                j = list(abs(t_seeds['time'] - dtime_samples)).index(min_diff)
                print('Detection at SNR of: ' + str(t_seeds['SNR'][j]))
                ktrue += 1
            else:
                print('Detection at sample: ' + str(dtime_samples) +
                      ' does not match anything in seed times:')
                kfalse += 1
            print 'Minimum difference in samples is: ' + str(min_diff)
    # Plot the detections
    if debug > 3:
        for i, template in enumerate(templates):
            times = [d.detect_time.datetime for d in detections
                     if d.template_name == template_names[i]]
            print(times)
            plotting.detection_multiplot(data, template, times)
    # Set an 'acceptable' ratio of positive to false detections
    print(str(ktrue) + ' true detections and ' + str(kfalse) +
          ' false detections')
    if kfalse / ktrue < 0.25:
        return True
    else:
        return False


def generate_synth_data(nsta=5, ntemplates=3, nseeds=100, samp_rate=20.0,
                        t_length=3.0, max_amp=10.0, debug=0):
    """
    Function to generate a synthetic dataset to be used for testing.  This will
    generate both templates and data to scan through.  Templates will be
    generated using the utils.synth_seis functions.  The day of data will be
    random noise, with random signal-to-noise ratio copies of the templates
    randomly seeded throughout the day.  It also returns the seed times and
    signal-to-noise ratios used.

    :type nsta: int
    :param nsta: Number of stations to generate data for < 15.
    :type ntemplates: int
    :param ntemplates: Number of templates to generate, will be generated with\
        random arrival times.
    :type nseeds: int
    :param nseeds: Number of copies of the template to seed within the day of\
        noisey data for each template.
    :type samp_rate: float
    :param samp_rate: Sampling rate to use in Hz
    :type t_length: float
    :param t_length: Length of templates in seconds.
    :type max_amp: float
    :param max_amp: Maximum signal-to-noise ratio of seeds.
    :type debug: int
    :param debug: Debug level, bigger the number, the more plotting/output.

    :returns: Templates: List of obspy.Stream, Data: obspy.Stream of seeded\
        noisy data, Seeds: dictionary of seed SNR and time with time in\
        samples.
    """
    from eqcorrscan.utils import synth_seis
    import numpy as np
    from obspy import UTCDateTime

    # Generate random arrival times
    t_times = np.abs(np.random.random([nsta, ntemplates])) * t_length
    # Generate random node locations - these do not matter as they are only
    # used for naming
    lats = np.random.random(ntemplates) * 90.0
    lons = np.random.random(ntemplates) * 90.0
    depths = np.abs(np.random.random(ntemplates) * 40.0)
    nodes = zip(lats, lons, depths)
    # Generating a 5x3 array to make 3 templates
    stations = ['ALPH', 'BETA', 'GAMM', 'KAPP', 'ZETA', 'BOB', 'MAGG',
                'ALF', 'WALR', 'ALBA', 'PENG', 'BANA', 'WIGG', 'SAUS',
                'MALC']
    if debug > 1:
        print(nodes)
        print(t_times)
        print(stations[0:nsta])
    templates = synth_seis.template_grid(stations=stations[0:nsta],
                                         nodes=nodes,
                                         travel_times=t_times, phase='S',
                                         samp_rate=samp_rate,
                                         flength=int(t_length * samp_rate))
    if debug > 2:
        for template in templates:
            print(template)
            template.plot(size=(800, 600), equal_scale=False)
    # Now we want to create a day of synthetic data
    seeds = []
    data = templates[0].copy()  # Copy a template to get the correct length
    # and stats for data, we will overwrite the data on this copy
    for tr in data:
        tr.data = np.zeros(86400 * int(samp_rate))
        # Set all the traces to have a day of zeros
        tr.stats.starttime = UTCDateTime(0)
    for i, template in enumerate(templates):
        impulses = np.zeros(86400 * int(samp_rate))
        # Generate a series of impulses for seeding
        # Need three seperate impulse traces for each of the three templates,
        # all will be convolved within the data though.
        impulse_times = np.random.randint(86400 * int(samp_rate),
                                          size=nseeds)
        impulse_amplitudes = np.random.randn(nseeds) * max_amp
        # Generate amplitudes up to maximum amplitude in a normal distribution
        seeds.append({'SNR': impulse_amplitudes,
                      'time': impulse_times})
        for j in range(nseeds):
            impulses[impulse_times[j]] = impulse_amplitudes[j]
        # We now have one vector of impulses, we need nsta numbers of them,
        # shifted with the appropriate lags
        mintime = min([template_tr.stats.starttime
                       for template_tr in template])
        for j, template_tr in enumerate(template):
            offset = int((template_tr.stats.starttime - mintime) * samp_rate)
            pad = np.zeros(offset)
            tr_impulses = np.append(pad, impulses)[0:len(impulses)]
            # Convolve this with the template trace to give the daylong seeds
            data[j].data += np.convolve(tr_impulses,
                                        template_tr.data)[0:len(impulses)]
        if debug > 2:
            data.plot(starttime=UTCDateTime(0) +
                      impulse_times[0]/samp_rate - 3,
                      endtime=UTCDateTime(0) +
                      impulse_times[0]/samp_rate + 15)
    # Add the noise
    for tr in data:
        noise = np.random.randn(86400 * int(samp_rate))
        tr.data += noise / max(noise)
        if debug > 2:
            tr.plot()

    return templates, data, seeds

if __name__ == '__main__':
    """
    Run core tests
    """
    test_match_filter()
