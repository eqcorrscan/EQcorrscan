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
from eqcorrscan import register_normxcorr
from eqcorrscan.utils.correlate import numpy_normxcorr, fftw_normxcorr, \
    time_multi_normxcorr, multichannel_normxcorr, normxcorr


# set seed state for consistent arrays
random = np.random.RandomState(13)


# -------------------------- helper functions


def gen_xcorr_func(name):
    """ return an xcorr function with desired name """

    def func(templates, stream, pads, *args, **kwargs):
        pass

    func.__name__ = name
    return func


def time_func(func, name, *args, **kwargs):
    """ call a func with args and kwargs, print name of func and how 
    long it took. """
    print('running %s' % name)
    tic = time.time()
    out = func(*args, **kwargs)
    toc = time.time()
    print('%s took %0.2f seconds' % (name, toc - tic))
    return out


def measure_counts(self, func):
    """ decorator for counter how often func get called """

    @wraps(func)
    def wrapper(*args, **kwargs):
        self.counter[func.__name__] += 1
        return func(*args, **kwargs)

    return wrapper


# ----------------------------- module fixtures


# auto run fixtures

@pytest.fixture(scope='class', autouse=True)
def swap_registery():
    """ copy the current registry, restore it when tests finish. This will
     run after every class test suite """
    current = copy.deepcopy(eqcorrscan.utils.correlate.XCOR_FUNCS)
    yield
    eqcorrscan.utils.correlate.XCOR_FUNCS = current

# array fixtures

starting_index = 500


@pytest.fixture(scope='module')
def array_template():
    """ 
    return a set of templates, generated with randomn, for correlation tests.
    """
    return random.randn(200, 200)


@pytest.fixture(scope='module')
def array_stream(array_template):
    """
    Return a stream genearted with randomn for testing normxcorr functions
    """
    stream = random.randn(10000) * 5
    # insert a template into the stream so cc == 1 at some place
    ar = array_template[0]
    stream[starting_index: starting_index + len(ar)] = ar
    return stream


@pytest.fixture(scope='module')
def pads(array_template):
    """
    return an array of zeros for padding, matching templates len.
    """
    return np.zeros(array_template.shape[0], dtype=int)


@pytest.fixture(scope='module')
def array_ccs(array_template, array_stream, pads):
    """  use each function stored in the normcorr cache to correlate the 
     templates and arrays, return a dict with keys as func names and values
     as the cc calculated by said function"""
    out = {}
    for name, func in eqcorrscan.utils.correlate.XCOR_FUNCS.items():
        cc, _ = time_func(func, name, array_template, array_stream, pads)
        out[name] = cc
    return out

# stream fixtures

chans = ['EHZ', 'EHN', 'EHE']
stas = ['COVA', 'FOZ', 'LARB', 'GOVA', 'MTFO', 'MTBA']
n_templates = 20
stream_len = 100000
template_len = 200


@pytest.fixture(scope='module')
def multichannel_templates():
    """ create multichannel templates """
    templates = []
    for i in range(n_templates):
        template = Stream()
        for station in stas:
            for channel in chans:
                template += Trace(data=random.randn(template_len))
                template[-1].stats.channel = channel
                template[-1].stats.station = station
        templates.append(template)
    return templates


@pytest.fixture(scope='module')
def multichannel_stream():
    """ create a multichannel stream for tests """
    stream = Stream()
    for station in stas:
        for channel in chans:
            stream += Trace(data=random.randn(stream_len))
            stream[-1].stats.channel = channel
            stream[-1].stats.station = station
    return stream


# ----------------------------------- tests


class TestCorrelateFunctionsReturnSame:
    """ these tests ensure the various implementations of normxcorr return 
    approximately the same answers """

    # tests
    def test_single_channel_similar(self, array_ccs):
        """ ensure each of the correlation methods return similar answers 
        given the same input data """
        cc_list = list(array_ccs.values())
        for cc1, cc2 in itertools.combinations(cc_list, 2):
            assert np.allclose(cc1, cc2, atol=0.05)

    def test_test_autocorrelation(self, array_ccs):
        """ ensure an auto correlationoccurred in each of ccs where it is 
        expected, defined by starting_index variable """
        for name, cc in array_ccs.items():
            assert np.isclose(cc[0, starting_index], 1., atol=.01)

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

    # fixtures
    @pytest.fixture
    def instrumented_default_normxcorr(self, monkeypatch, array_template,
                                       array_stream, pads):
        """ instrument the default xcor """
        default = eqcorrscan.utils.correlate.XCOR_FUNCS['default']
        func_name = default.__name__
        monkeypatch.setitem(eqcorrscan.utils.correlate.XCOR_FUNCS,
                            'default', measure_counts(self, default))
        _ = normxcorr(array_template, array_stream, pads)
        yield func_name
        self.__class__.counter = defaultdict(lambda: 0)  # reset counter

    @pytest.fixture
    def swapped_normxcorr(self, monkeypatch, array_template, array_stream,
                          pads):
        """ ensure the default normxcorr can be changed  """
        new_func = measure_counts(self, time_multi_normxcorr)
        monkeypatch.setitem(eqcorrscan.utils.correlate.XCOR_FUNCS,
                            'time_domain', new_func)
        _ = normxcorr(array_template, array_stream, pads,
                      xcorr_func='time_domain')
        yield new_func.__name__
        self.__class__.counter = defaultdict(lambda: 0)  # reset counter

    @pytest.fixture
    def callable_normxcorr(self, array_template, array_stream, pads):
        def customfunc(templates, stream, pads, *args, **kwargs):
            pass

        func = measure_counts(self, customfunc)

        _ = normxcorr(array_template, array_stream, pads, xcorr_func=func)
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

    # tests
    def test_register_as_decorator_no_args(self):
        """ ensure register_normxcorr works as a decorator with no args """

        @register_normxcorr
        def func1(templates, stream, pads, *args, **kwargs):
            pass

        assert self.name_func_is_registered(func1)

    def test_register_as_decorator_with_args(self):
        """ ensure register can be used as a decorator with args """

        @register_normxcorr(name='func2')
        def func(templates, stream, pads, *args, **kwargs):
            pass

        assert self.name_func_is_registered('func2')

    def test_register_as_callable(self):
        """ ensure register can be used as a callable to take a name
        and a normxcorr func """
        func = gen_xcorr_func('funky')
        register_normxcorr(name='func3', func=func)
        assert self.name_func_is_registered('func3')

    def test_set_default(self):
        """ ensure the default can be overwritten """
        func = gen_xcorr_func('funky')
        eqcorrscan.utils.correlate.register_normxcorr(func, is_default=True)
        assert eqcorrscan.utils.correlate.XCOR_FUNCS['default'] is func


class TestRegisterAlternativeConcurrency:
    """ Tests for registering alternative concurrency functions """
    counter = defaultdict(lambda: 0)

    # helper functions
    def new_multi(self, template_dict, stream_dict, pad_dict, seed_ids):
        pass

    # fixtures
    @pytest.fixture
    def r_normxcorr(self):
        """ return a registered normxcorr function """
        return register_normxcorr(gen_xcorr_func('normxcorr'))

    @pytest.fixture
    def normxcorr_new_multithread(self, r_normxcorr):
        """ register the new multithread method """
        func = measure_counts(self, self.new_multi)
        r_normxcorr.register('multithread')(func)
        r_normxcorr.multithread(None, None, None, None)
        yield func.__name__
        self.counter.pop(func.__name__)

    # tests
    def test_new_method_was_called(self, normxcorr_new_multithread):
        """ ensure the new method was called """
        assert self.counter[normxcorr_new_multithread]


class TestAlternativeConcurrency:
    """ tests for alternative concurrency functions """

    # fixtures
    @pytest.fixture(scope='class')
    def multichan_normxcorr(self, multichannel_templates,
                            multichannel_stream):
        """ return the result given by default multichannel normxcorr """
        cc, _, _ = multichannel_normxcorr(multichannel_templates,
                                          multichannel_stream)
        return cc

    @pytest.fixture(scope='class')
    def numpy_serial(self, multichannel_templates, multichannel_stream):
        """ return the output of the numpy functions serial interface """
        cc = numpy_normxcorr.serial(multichannel_templates,
                                    multichannel_stream)
        return cc

    # tests
    def test_numpy_serial_is_similar(self, multichan_normxcorr,
                                     numpy_serial):
        """ delete when passes """
        assert multichan_normxcorr.shape == numpy_serial.shape
        assert np.allclose(numpy_serial, multichan_normxcorr)



if __name__ == '__main__':
    """
    Run core tests
    """
    unittest.main()
