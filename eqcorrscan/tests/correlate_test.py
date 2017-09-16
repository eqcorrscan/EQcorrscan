from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import itertools
import time
from collections import defaultdict
from functools import wraps

import numpy as np
import pytest
from obspy import Trace, Stream

import eqcorrscan.utils.correlate as corr
from eqcorrscan.utils.correlate import register_array_xcorr

# set seed state for consistent arrays
random = np.random.RandomState(7)


# -------------------------- helper functions


def gen_xcorr_func(name):
    """ return an xcorr function with desired name """

    def func(templates, stream, pads, *args, **kwargs):
        pass

    func.__name__ = str(name)
    return func


def time_func(func, name, *args, **kwargs):
    """ call a func with args and kwargs, print name of func and how
    long it took. """
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


# generate streams for templates and search space

chans = ['EHZ', 'EHN', 'EHE']
stas = ['COVA', 'FOZ', 'LARB', 'GOVA', 'MTFO', 'MTBA']
n_templates = 20
stream_len = 100000
template_len = 200


def generate_multichannel_stream():
    stream = Stream()
    for station in stas:
        for channel in chans:
            stream += Trace(data=random.randn(stream_len))
            stream[-1].stats.channel = channel
            stream[-1].stats.station = station
    return stream


def generate_multichannel_templates():
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


# ----------------------------- module fixtures


# auto run fixtures

@pytest.fixture(scope='class', autouse=True)
def print_class_name(request):
    """ prints the class name before a class test suite starts """
    try:
        cls_name = request.cls.__name__
    except AttributeError:  # test does not belong to a class
        return
    else:
        dash = '-' * 70
        print('\nstarting tests contained in class %s\n%s' % (cls_name, dash))


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

    for name in list(corr.XCORR_FUNCS_ORIGINAL.keys()):
        func = corr.get_array_xcorr(name)
        print("Running %s" % name)
        cc, _ = time_func(func, name, array_template, array_stream, pads)
        out[name] = cc
    return out


# stream fixtures


@pytest.fixture(scope='module')
def multichannel_templates():
    """ create multichannel templates """
    return generate_multichannel_templates()


@pytest.fixture(scope='module')
def multichannel_stream():
    """ create a multichannel stream for tests """
    return generate_multichannel_stream()


# a dict of all registered stream functions (this is a bit long)
stream_funcs = {fname + '_' + mname: corr.get_stream_xcorr(fname, mname)
                for fname in corr.XCORR_FUNCS_ORIGINAL.keys()
                for mname in corr.XCORR_STREAM_METHODS
                if fname != 'default'}


@pytest.fixture(scope='module')
def stream_cc_output_dict(multichannel_templates, multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    out = {}
    for name, func in stream_funcs.items():
        cc_out = time_func(func, name, multichannel_templates,
                           multichannel_stream, cores=1)
        out[name] = cc_out
    return out


@pytest.fixture(scope='module')
def stream_cc_dict(stream_cc_output_dict):
    """ return just the cc arrays from the stream_cc functions """
    return {name: result[0] for name, result in stream_cc_output_dict.items()}


# ----------------------------------- tests


class TestArrayCorrelateFunctions:
    """ these tests ensure the various implementations of normxcorr return
    approximately the same answers """
    atol = .00001  # how close correlations have to be

    # tests
    def test_single_channel_similar(self, array_ccs):
        """ ensure each of the correlation methods return similar answers
        given the same input data """
        cc_list = list(array_ccs.values())
        for cc1, cc2 in itertools.combinations(cc_list, 2):
            assert np.allclose(cc1, cc2, atol=self.atol)

    def test_test_autocorrelation(self, array_ccs):
        """ ensure an auto correlationoccurred in each of ccs where it is
        expected, defined by starting_index variable """
        for name, cc in array_ccs.items():
            assert np.isclose(cc[0, starting_index], 1., atol=self.atol)


@pytest.mark.serial
class TestStreamCorrelateFunctions:
    """ same thing as TestArrayCorrelateFunction but for stream interface """
    atol = TestArrayCorrelateFunctions.atol

    def test_multi_channel_xcorr(self, stream_cc_dict):
        """ test various correlation methods with multiple channels """
        # get correlation results into a list
        cc_names = list(stream_cc_dict.keys())
        cc_list = [stream_cc_dict[cc_name] for cc_name in cc_names]
        cc_1 = cc_list[0]
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            assert np.allclose(cc_1, cc, atol=self.atol)


class TestXcorrContextManager:
    # fake_cache = copy.deepcopy(corr.XCOR_FUNCS)

    @pytest.fixture
    def cache(self):
        """ this fixtures resets the class level cache after every test """
        yield copy.deepcopy(corr.XCOR_FUNCS)

    @pytest.fixture
    def set_xcorr(self, cache):
        return corr._Context(cache, 'default')

    @pytest.fixture
    def set_value(self, set_xcorr):
        """ set a value in a 'permanent' fashion, return the set function """

        def func(templates, stream, pads):
            pass

        set_xcorr(func)
        return func

    # tests
    def test_value_was_set(self, set_value, cache):
        assert cache['default'] == set_value

    def test_context_manager(self, cache):
        """ ensure the context manager reverts changes """
        context = corr._Context(cache, 'default')
        old_default = cache['default']
        new_val = corr.numpy_normxcorr
        with context(new_val):
            assert cache['default'] == new_val
        assert old_default

    def test_str_accepted(self):
        """ ensure a str of the xcorr function can be passed as well """
        with corr.set_xcorr('numpy'):
            func = corr.get_array_xcorr()
            assert func is corr.numpy_normxcorr


class TestGenericStreamXcorr:
    """ tests for stream_xocrr function """

    # tests
    def test_noargs_returns_default(self):
        """ ensure passing no args to get_stream_xcorr returns default """
        func = corr.get_stream_xcorr()
        default = corr.XCOR_FUNCS['default'].stream_xcorr
        assert func is default

    def test_callable_registered(self, multichannel_templates,
                                 multichannel_stream):
        """ ensure a callable can be registered """
        small_count = {}

        def some_callable(template_array, stream_array, pad_array):
            small_count['name'] = 1
            return corr.numpy_normxcorr(template_array, stream_array,
                                        pad_array)

        func = corr.get_stream_xcorr(some_callable)
        func(multichannel_templates, multichannel_stream)
        assert 'name' in small_count

    def test_bad_concurrency_raises(self):
        """ ensure passing an invalid concurrency argument raises a
        ValueError"""
        with pytest.raises(ValueError):
            corr.get_stream_xcorr(concurrency='node.js')

    def test_loading_unregistered_function_registers(self):
        """ ensure if a function in cache hasn't been decoratored it gets
        decorated when returned """

        def func(templates, streams, pads):
            pass

        corr.XCOR_FUNCS['func_test'] = func
        corr.get_stream_xcorr('func_test')
        assert hasattr(corr.XCOR_FUNCS['func_test'], 'registered')

    def test_using_custom_function_doesnt_change_default(self):
        """ ensure a custom function will not change the default """

        def func(templates, streams, pads):
            pass

        default = corr.get_array_xcorr(None)

        corr.get_array_xcorr(func)

        assert corr.get_array_xcorr(None) is default


class TestRegisterNormXcorrs:
    """ Tests for register_normxcorr function, which holds global context
    for which xcorr to use """

    # helper functions
    def name_func_is_registered(self, func_name):
        """ return True if func is registered as a normxcorr func """
        # Note: don not remove this fixture or bad things will happen
        name = func_name.__name__ if callable(func_name) else func_name
        return name in corr.XCOR_FUNCS

    # fixtures
    @pytest.fixture(scope='class', autouse=True)
    def swap_registery(self):
        """ copy the current registry, restore it when tests finish"""
        current = copy.deepcopy(corr.XCOR_FUNCS)
        yield
        corr.XCOR_FUNCS = current

    # tests
    def test_register_as_decorator_no_args(self):
        """ ensure register_normxcorr works as a decorator with no args """

        @register_array_xcorr
        def func1(templates, stream, pads, *args, **kwargs):
            pass

        assert self.name_func_is_registered(func1)

    def test_register_as_decorator_with_args(self):
        """ ensure register can be used as a decorator with args """

        @register_array_xcorr(name='func2')
        def func(templates, stream, pads, *args, **kwargs):
            pass

        assert self.name_func_is_registered('func2')

    def test_register_as_callable(self):
        """ ensure register can be used as a callable to take a name
        and a normxcorr func """
        func = gen_xcorr_func('funky')
        register_array_xcorr(name='func3', func=func)
        assert self.name_func_is_registered('func3')

    def test_set_default(self):
        """ ensure the default can be overwritten """
        func = gen_xcorr_func('funky')
        corr.register_array_xcorr(func, is_default=True)
        assert corr.XCOR_FUNCS['default'] is func

    def test_register_bad_func_rasies(self):
        """ test trying to register a non-supported function raises """
        func = corr.XCOR_FUNCS['default']

        with pytest.raises(ValueError):
            @func.register('not_supported_value')
            def func():
                pass


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
        return register_array_xcorr(gen_xcorr_func('normxcorr'))

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
