from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import unittest
import numpy as np
import time

from obspy import Trace, Stream

from eqcorrscan.utils.correlate import numpy_normxcorr, fftw_normxcorr, \
    time_multi_normxcorr, multichannel_normxcorr


class CorrelateTests(unittest.TestCase):
    def test_same_various_methods(self):
        templates = np.random.randn(200, 200)
        stream = np.random.randn(10000)
        stream *= stream ** 5
        pads = np.zeros(templates.shape[0], dtype=int)
        tic = time.time()
        numpy_ccc, no_chans = numpy_normxcorr(templates, stream, pads)
        toc = time.time()
        print('Scipy took: %f seconds' % (toc-tic))
        tic = time.time()
        fftw_ccc, no_chans = fftw_normxcorr(
            templates, stream, pads, threaded=False)
        toc = time.time()
        print('FFTW took: %f seconds' % (toc-tic))
        tic = time.time()
        fftw_ccc_2d, no_chans = fftw_normxcorr(
            templates, stream, pads, threaded=True)
        toc = time.time()
        print('FFTW-threaded took: %f seconds' % (toc-tic))
        tic = time.time()
        time_ccc, no_chans = time_multi_normxcorr(templates, stream, pads)
        toc = time.time()
        print('Time-domain took: %f seconds' % (toc-tic))
        self.assertTrue(np.allclose(numpy_ccc, fftw_ccc, atol=0.001))
        self.assertTrue(np.allclose(numpy_ccc, fftw_ccc_2d, atol=0.001))
        self.assertTrue(np.allclose(numpy_ccc, time_ccc, atol=0.12))
        self.assertTrue(np.allclose(fftw_ccc, time_ccc, atol=0.12))

    def test_multi_channel_xcorr(self):
        chans = ['EHZ', 'EHN', 'EHE']
        stas = ['COVA', 'FOZ', 'LARB', 'GOVA', 'MTFO', 'MTBA']
        n_templates = 20
        stream_len = 100000
        template_len = 200
        templates = []
        stream = Stream()
        for station in stas:
            for channel in chans:
                stream += Trace(data=np.random.randn(stream_len))
                stream[-1].stats.channel = channel
                stream[-1].stats.station = station
        for i in range(n_templates):
            template = Stream()
            for station in stas:
                for channel in chans:
                    template += Trace(data=np.random.randn(template_len))
                    template[-1].stats.channel = channel
                    template[-1].stats.station = station
            templates.append(template)
        # print("Running time serial")
        # tic = time.time()
        # cccsums_t_s, no_chans, chans = multichannel_normxcorr(
        #     templates=templates, stream=stream, cores=None, time_domain=True)
        # toc = time.time()
        # print('Time-domain in serial took: %f seconds' % (toc-tic))
        print("Running time parallel")
        tic = time.time()
        cccsums_t_p, no_chans, chans = multichannel_normxcorr(
            templates=templates, stream=stream, cores=4, time_domain=True)
        toc = time.time()
        print('Time-domain in parallel took: %f seconds' % (toc-tic))
        print("Running frequency serial")
        tic = time.time()
        cccsums_f_s, no_chans, chans = multichannel_normxcorr(
            templates=templates, stream=stream, cores=None, time_domain=False)
        toc = time.time()
        print('Frequency-domain in serial took: %f seconds' % (toc-tic))
        print("Running frequency parallel")
        tic = time.time()
        cccsums_f_p, no_chans, chans = multichannel_normxcorr(
            templates=templates, stream=stream, cores=4, time_domain=False)
        toc = time.time()
        print('Frequency-domain in parallel took: %f seconds' % (toc-tic))
        print("Running frequency openmp parallel")
        tic = time.time()
        cccsums_f_op, no_chans, chans = multichannel_normxcorr(
            templates=templates, stream=stream, cores=2, time_domain=False,
            openmp=True)
        toc = time.time()
        print('Frequency-domain in parallel took: %f seconds' % (toc-tic))
        print("Finished")
        # self.assertTrue(np.allclose(cccsums_t_s, cccsums_t_p, atol=0.00001))
        self.assertTrue(np.allclose(cccsums_f_s, cccsums_f_p, atol=0.00001))
        self.assertTrue(np.allclose(cccsums_f_s, cccsums_f_op, atol=0.00001))
        if not np.allclose(cccsums_t_p, cccsums_f_s, atol=0.7):
            import matplotlib.pyplot as plt
            plt.plot(cccsums_t_p[0], 'k', label='Time')
            plt.plot(cccsums_f_s[0], 'r', label='Frequency')
            plt.legend()
            plt.show()
        self.assertTrue(np.allclose(cccsums_t_p, cccsums_f_s, atol=0.7))


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
