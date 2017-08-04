from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import itertools
import time
import unittest
from collections import defaultdict
from functools import partial, wraps

import numpy as np
import pytest
from obspy import Trace, Stream

import eqcorrscan.utils.correlate
from eqcorrscan.utils.correlate import numpy_normxcorr, fftw_normxcorr, \
    time_multi_normxcorr, multichannel_normxcorr, normxcorr

# set seed state for consistent arrays
random = np.random.RandomState(13)


# helper functions for test suite

def time_func(func, name, *args, **kwargs):
    """ call a func with args and kwargs, print name of func and how 
    long it took. """
    print('running %s' % name)
    tic = time.time()
    out = func(*args, **kwargs)
    toc = time.time()
    print('%s took %0.2f seconds' % (name, toc - tic))
    return out


# module fixtures

@pytest.fixture(scope='module')
def templates():
    """ 
    return a set of templates, generated with randomn, for correlation tests.
    """
    # setup search space and templates
    return random.randn(200, 200)


@pytest.fixture(scope='module')
def stream():
    """
    Return a stream genearted with randomn for testing normxcorr functions
    """
    return random.randn(10000) * 5


@pytest.fixture(scope='module')
def pads(templates):
    """
    return an array of zeros for padding, matching templates len.
    :return: 
    """
    return np.zeros(templates.shape[0], dtype=int)


class TestCorrelate:
    # contants for setup
    chans = ['EHZ', 'EHN', 'EHE']
    stas = ['COVA', 'FOZ', 'LARB', 'GOVA', 'MTFO', 'MTBA']
    n_templates = 20
    stream_len = 100000
    template_len = 200

    # fixtures
    @pytest.fixture
    def multichannel_stream(self):
        """ create a multichannel stream for tests """
        stream = Stream()
        for station in self.stas:
            for channel in self.chans:
                stream += Trace(data=random.randn(self.stream_len))
                stream[-1].stats.channel = channel
                stream[-1].stats.station = station
        return stream

    @pytest.fixture
    def multichannel_templates(self):
        """ create multichannel templates """
        templates = []
        for i in range(self.n_templates):
            template = Stream()
            for station in self.stas:
                for channel in self.chans:
                    template += Trace(data=random.randn(self.template_len))
                    template[-1].stats.channel = channel
                    template[-1].stats.station = station
            templates.append(template)
        return templates

    # tests
    def test_same_various_methods(self, templates, stream, pads):
        """ basic test case for each of the correlation methods """
        # functions to test
        cc_funcs = dict(
            scipy=numpy_normxcorr,
            fftw=partial(fftw_normxcorr, threaded=False),
            fftw_threaded=partial(fftw_normxcorr, threaded=True),
            time_domain=time_multi_normxcorr,
        )

        # loop over cc funcs and store results in cc_list
        cc_list = []
        for name, func in cc_funcs.items():
            cc, _ = time_func(func, name, templates, stream, pads)
            cc_list.append(cc)

        # compare each cc against other ccs
        for cc1, cc2 in itertools.combinations(cc_list, 2):
            assert np.allclose(cc1, cc2, atol=0.05)

    def test_multi_channel_xcorr(self, multichannel_templates,
                                 multichannel_stream):
        """ test various correlation methods with multiple channels """
        # collection functions to test
        cc_funcs = dict(
            time_parrallel=partial(multichannel_normxcorr, time_domain=True,
                                   cores=4),
            frequency_serial=partial(multichannel_normxcorr, cores=None,
                                     time_domain=False),
            frequency_parallel=partial(multichannel_normxcorr, cores=4,
                                       time_domain=False),
            freq_open_mp=partial(multichannel_normxcorr, cores=2,
                                 time_domain=False, openmp=True),
        )

        # loop cc funcs and store results
        cc_list = []
        for name, func in cc_funcs.items():
            ccs, _, _ = time_func(func, name,
                                  templates=multichannel_templates,
                                  stream=multichannel_stream)
            cc_list.append(ccs)

        # compare combinations of ccs
        for cc1, cc2 in itertools.combinations(cc_list, 2):
            assert np.allclose(cc1, cc2, atol=0.05)


class TestGenericNormxcorr:
    """ tests for generic normxcorr that acts as facade to xcorr functions """

    counter = defaultdict(lambda: 0)  # a simple counter

    # helper functions
    def measure_counts(self, func):
        """ decorator for counter how often func get called """

        @wraps(func)
        def wrapper(*args, **kwargs):
            self.counter[func.__name__] += 1
            return func(*args, **kwargs)

        return wrapper

    # fixtures
    @pytest.fixture
    def instrumented_default_normxcorr(self, monkeypatch, templates, stream,
                                       pads):
        """ instrument the default xcor """
        default = eqcorrscan.utils.correlate.XCOR_FUNCS['default']
        func_name = default.__name__
        monkeypatch.setitem(eqcorrscan.utils.correlate.XCOR_FUNCS,
                            'default', self.measure_counts(default))
        _ = normxcorr(templates, stream, pads)
        yield func_name
        self.__class__.counter = defaultdict(lambda: 0)  # reset counter

    @pytest.fixture
    def swapped_normxcorr(self, monkeypatch, templates, stream, pads):
        """ ensure the default normxcorr can be changed  """
        new_func = self.measure_counts(time_multi_normxcorr)
        monkeypatch.setitem(eqcorrscan.utils.correlate.XCOR_FUNCS,
                            'time_domain', new_func)
        _ = normxcorr(templates, stream, pads, xcorr_func='time_domain')
        yield new_func.__name__
        self.__class__.counter = defaultdict(lambda: 0)  # reset counter

    @pytest.fixture
    def callable_normxcorr(self, templates, stream, pads):
        def customfunc(templates, stream, pads, *args, **kwargs):
            pass

        func = self.measure_counts(customfunc)

        _ = normxcorr(templates, stream, pads, xcorr_func=func)
        return func.__name__

    # tests
    def test_normxcorr_calls_default(self, instrumented_default_normxcorr):
        """ ensure the default normxcorr function gets called when 
         xcor_func is not used """
        assert self.counter[instrumented_default_normxcorr] == 1

    def test_swap_corr_function(self, swapped_normxcorr):
        """ ensure the xcor_func can be used to change xcor func """
        assert self.counter[swapped_normxcorr] == 1

    def test_xcor_func_callable(self, callable_normxcorr):
        """ ensure a custom function passed to normxcorr gets called """
        assert self.counter[callable_normxcorr] == 1


class TestRegisterNormXcorrs:
    """ Tests for register_normxcorr function, which holds global context
    for which xcorr to use """

    # helper functions
    def name_func_is_registered(self, func_name):
        """ return True if func is registered as a normxcorr func """
        name = func_name.__name__ if callable(func_name) else func_name
        return name in eqcorrscan.utils.correlate.XCOR_FUNCS

    def gen_xcorr_func(self, name):
        """ return an xcorr function with desired name """

        def func(templates, stream, pads, *args, **kwargs):
            pass

        func.__name__ = name
        return func

    # fixtures
    @pytest.fixture(scope='class', autouse=True)
    def swap_registery(self):
        """ get a copy of the current registry, restore it when test finish """
        current = copy.deepcopy(eqcorrscan.utils.correlate.XCOR_FUNCS)
        yield
        eqcorrscan.utils.correlate.XCOR_FUNCS = current

    # tests
    def test_register_as_decorator_no_args(self):
        """ ensure register_normxcorr works as a decorator with no args """

        @eqcorrscan.utils.correlate.register_normxcorr
        def func1(templates, stream, pads, *args, **kwargs):
            pass

        assert self.name_func_is_registered(func1)

    def test_register_as_decorator_with_args(self):
        """ ensure register can be used as a decorator with args """

        @eqcorrscan.utils.correlate.register_normxcorr(name='func2')
        def func(templates, stream, pads, *args, **kwargs):
            pass

        assert self.name_func_is_registered('func2')

    def test_register_as_callable(self):
        """ ensure register can be used as a callable to take a name
        and a normxcorr func """
        func = self.gen_xcorr_func('funky')

        eqcorrscan.utils.correlate.register_normxcorr(name='func3', func=func)
        assert self.name_func_is_registered('func3')

    def test_set_default(self):
        """ ensure the default can be overwritten """
        func = self.gen_xcorr_func('funky')
        eqcorrscan.utils.correlate.register_normxcorr(func, is_default=True)
        assert eqcorrscan.utils.correlate.XCOR_FUNCS['default'] is func


if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
