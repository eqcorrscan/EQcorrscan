"""
A series of test functions for the core functions in EQcorrscan.
"""


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
    from obspy import Stream, Trace
    from eqcorrscan.utils import synth_seis
    import numpy as np

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
    print(templates)
    if debug > 2:
        for template in templates:
            print(template)
            template.plot(size=(800, 600), equal_scale=False)
    # Now we want to create a day of synthetic data
    seeds = []
    data = templates[0].copy()  # Copy a template to get the correct length
    # and stats for data, we will overwrite the data on this copy
    impulses = np.zeros(86400 * int(samp_rate))
    for tr in data:
        tr.data = impulses  # Set all the traces to have a day of zeros
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
            if debug > 0:
                data[j].plot()
    # Add the noise
    for tr in data:
        tr.data += np.random.randn(86400 * int(samp_rate))
        if debug > 2:
            tr.plot()
    if debug > 0:
        data.plot()

    return templates, data, seeds
