"""
Test of large-scale match-filter, not to be run on CI.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from eqcorrscan.utils import synth_seis
from eqcorrscan.utils.timer import Timer
from eqcorrscan.core.match_filter import match_filter


@profile
def test_synth_large():
    print('\tGenerating Synthetic data\n\n')
    templates, data, seeds = synth_seis.generate_synth_data(
        nsta=10, ntemplates=20, nseeds=5, samp_rate=100, t_length=6,
        max_amp=5, max_lag=5, debug=0)
    print('\tRunning the parallel detections\n\n')
    with Timer() as t:
        detections = match_filter(
            template_names=[str(i) for i in range(len(templates))],
            template_list=templates, st=data, threshold=8, threshold_type='MAD',
            trig_int=6, plotvar=False, cores=4, output_event=False)
    print('Parallel run took %f seconds' % t.secs)
    print('\tRunning the serial detections\n\n')
    with Timer() as t:
        detections = match_filter(
            template_names=[str(i) for i in range(len(templates))],
            template_list=templates, st=data, threshold=8, threshold_type='MAD',
            trig_int=6, plotvar=False, cores=None, output_event=False)
    print('Serial run took %f seconds' % t.secs)


if __name__ == '__main__':
    """
    Run core tests
    """
    test_synth_large()
