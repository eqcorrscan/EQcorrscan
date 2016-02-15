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


def generate_synth_data():
    """
    Function to generate a synthetic dataset to be used for testing.
    """
    # First generate a synthetic template
    from eqcorrscan.utils import synth_seis
    from obspy import read
    import numpy as np
    import matplotlib.pyplot as plt
    t_times = np.array([np.random.random(5), np.random.random(5),
                        np.random.random(5)]).reshape(5, 3)
    # Generating a 5x3 array to make 3 templates
    templates = synth_seis.template_grid(stations=['ALPH', 'BETA', 'GAMM',
                                                   'KAPP', 'ZETA'],
                                         nodes=[(-42.5, 35.0, 10.0)],
                                         travel_times=t_times, phase='S')
    # Now we want to create a day of synthetic data
    noise = np.random.random(8640000)  # Generate empty arrays to be seeded
    impulses = np.zeros(8640000)  # Generate a series of impulses for seeding
    # Need three seperate impulse traces for each of the three templates, all
    # will be convolved within the data though.
