import copy
import itertools
import logging
from collections import defaultdict
from functools import wraps
from os.path import join

import numpy as np
import pytest
from multiprocessing import cpu_count
from obspy import Trace, Stream, read
from obspy.core.util import AttribDict
from scipy.fftpack import next_fast_len

import eqcorrscan.utils.correlate as corr
from eqcorrscan.utils.correlate import register_array_xcorr
from eqcorrscan.utils.timer import time_func
from eqcorrscan.helpers.mock_logger import MockLoggingHandler

# set seed state for consistent arrays
random = np.random.RandomState(7)

log = logging.getLogger(corr.__name__)
_log_handler = MockLoggingHandler(level='DEBUG')
log.addHandler(_log_handler)
log_messages = _log_handler.messages

# -------------------------- helper functions


def gen_xcorr_func(name):
    """ return an xcorr function with desired name """

    def func(templates, stream, pads, *args, **kwargs):
        pass

    func.__name__ = str(name)
    return func


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
unstacked_stream_len = 10000
# Use a reduced length for unstacked to conserve memory
template_len = 200
gap_start = 5000
# fft_len = 2 ** 13
fft_len = next_fast_len(stream_len + template_len + 1)


def generate_multichannel_stream():
    stream = Stream()
    for station in stas:
        for channel in chans:
            stream += Trace(data=random.randn(stream_len))
            stream[-1].stats.channel = channel
            stream[-1].stats.station = station
    return stream


def generate_gappy_multichannel_stream():
    stream = generate_multichannel_stream()
    gappy_station = stream.select(station="COVA").copy()
    for tr in stream.select(station="COVA"):
        stream.remove(tr)
    gappy_station = gappy_station.cutout(
        stream[0].stats.starttime + gap_start,
        stream[0].stats.starttime + gap_start + (template_len * 4)).merge(
        fill_value=0)
    stream += gappy_station
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


def generate_weighted_multichannel_templates():
    random.seed(42)
    templates = generate_multichannel_templates()
    for template in templates:
        for tr in template:
            tr.stats.extra = AttribDict()
            tr.stats.extra.update({"weight": random.random()})
    return templates


def read_gappy_real_template():
    return [read(join(pytest.test_data_path, "DUWZ_template.ms"))]


def read_gappy_real_data():
    """ These data SUCK - gap followed by spike, and long period trend.
    Super fugly"""
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    from eqcorrscan.utils.pre_processing import multi_process

    client = Client("GEONET")
    st = client.get_waveforms(
        network="NZ", station="DUWZ", location="20", channel="BNZ",
        starttime=UTCDateTime(2016, 12, 31, 23, 58, 56),
        endtime=UTCDateTime(2017, 1, 1, 0, 58, 56))
    st = multi_process(
        st=st.merge(), lowcut=2, highcut=20, filt_order=4, samp_rate=50)
    return st


def read_real_multichannel_templates():
    import glob
    tutorial_templates = glob.glob(
        join(pytest.test_data_path, "tutorial_template_0.ms"))

    t = read(tutorial_templates[0])
    return [t.select(station="POWZ") + t.select(station="HOWZ")]


def get_real_multichannel_data():
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    from eqcorrscan.utils.pre_processing import multi_process

    t1 = UTCDateTime("2016-01-04T12:00:00.000000Z")
    t2 = t1 + 600
    bulk = [('NZ', 'POWZ', '*', 'EHZ', t1, t2),
            ('NZ', 'HOWZ', '*', 'EHZ', t1, t2)]
    client = Client("GEONET")
    st = client.get_waveforms_bulk(bulk)
    st = multi_process(
        st.merge(), lowcut=2.0, highcut=9.0, filt_order=4, samp_rate=20.0,
        starttime=t1, endtime=t2)
    return st

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
GLOBAL_WEIGHT = 0.5  # Weight applied to array ccs funcs
# global to allow comparison to unweighted


@pytest.fixture(scope='module')
def array_template():
    """
    return a set of templates, generated with random, for correlation tests.
    """
    return random.randn(200, 200)


@pytest.fixture(scope='module')
def array_stream(array_template):
    """
    Return a stream generated with random for testing normxcorr functions
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
def array_ccs_weighted(array_template, array_stream, pads):
    """ Use each function stored in the normxcorr cache to correlate the
     templates and arrays, return a dict with keys as func names and values
     as the cc calculated by said function"""
    out = {}
    weights = np.ones(array_template.shape[0]) * GLOBAL_WEIGHT
    for name in list(corr.XCORR_FUNCS_ORIGINAL.keys()):
        func = corr.get_array_xcorr(name)
        print("Running %s" % name)
        cc, _ = time_func(
            func, name, array_template, array_stream, pads, weights)
        out[name] = cc
        if "fftw" in name:
            print("Running fixed len fft")
            fft_len = next_fast_len(
                max(len(array_stream) // 4, len(array_template)))
            cc, _ = time_func(func, name, array_template, array_stream, pads,
                              weights, fft_len=fft_len)
            out[name + "_fixed_len"] = cc
    return out


@pytest.fixture(scope='module')
def array_ccs(array_template, array_stream, pads):
    """ Use each function stored in the normxcorr cache to correlate the
     templates and arrays, return a dict with keys as func names and values
     as the cc calculated by said function"""
    out = {}

    for name in list(corr.XCORR_FUNCS_ORIGINAL.keys()):
        func = corr.get_array_xcorr(name)
        print("Running %s" % name)
        cc, _ = time_func(func, name, array_template, array_stream, pads)
        out[name] = cc
        if "fftw" in name:
            print("Running fixed len fft")
            fft_len = next_fast_len(
                max(len(array_stream) // 4, len(array_template)))
            cc, _ = time_func(func, name, array_template, array_stream, pads,
                              fft_len=fft_len)
            out[name + "_fixed_len"] = cc
    return out


@pytest.fixture(scope='module')
def array_ccs_low_amp(array_template, array_stream, pads):
    """ Use each function stored in the normxcorr cache to correlate the
     templates and arrays, return a dict with keys as func names and values
     as the cc calculated by said function.
     This specifically tests low amplitude streams as raised in issue #181."""
    out = {}
    arr_stream = array_stream * 10e-8
    for name in list(corr.XCORR_FUNCS_ORIGINAL.keys()):
        func = corr.get_array_xcorr(name)
        print("Running {0} with low-variance".format(name))
        _log_handler.reset()
        cc, _ = time_func(
            func, name, array_template, arr_stream, pads)
        out[name] = (cc, copy.deepcopy(log_messages['warning']))
        if "fftw" in name:
            print("Running fixed len fft")
            _log_handler.reset()
            fft_len = next_fast_len(
                max(len(array_stream) // 4, len(array_template)))
            cc, _ = time_func(func, name, array_template, array_stream, pads,
                              fft_len=fft_len)
            out[name + "_fixed_len"] = (
                cc, copy.deepcopy(log_messages['warning']))
    return out

# stream fixtures


@pytest.fixture(scope='module')
def multichannel_templates():
    """ create multichannel templates """
    return generate_multichannel_templates()


@pytest.fixture(scope='module')
def weighted_multichannel_templates():
    """ create weighted multichannel templates """
    return generate_weighted_multichannel_templates()


@pytest.fixture(scope='module')
def multichannel_stream():
    """ create a multichannel stream for tests """
    return generate_multichannel_stream()


@pytest.fixture(scope='module')
def gappy_multichannel_stream():
    """ Create a multichannel stream with gaps (padded with zeros). """
    return generate_gappy_multichannel_stream()


@pytest.fixture(scope='module')
def gappy_real_data():
    return read_gappy_real_data()


@pytest.fixture(scope='module')
def gappy_real_data_template():
    return read_gappy_real_template()


@pytest.fixture(scope='module')
def real_templates():
    return read_real_multichannel_templates()


@pytest.fixture(scope='module')
def real_multichannel_stream():
    return get_real_multichannel_data()


# a dict of all registered stream functions (this is a bit long)
stream_funcs = {fname + '_' + mname: corr.get_stream_xcorr(fname, mname)
                for fname in sorted(corr.XCORR_FUNCS_ORIGINAL.keys())
                for mname in corr.XCORR_STREAM_METHODS
                if fname != 'default'}


def _stream_cc_output_dict(multichannel_templates, multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    out = {}
    for name, func in stream_funcs.items():
        for cores in [1, cpu_count()]:
            print("Running {0} with {1} cores".format(name, cores))

            cc_out = time_func(func, name, multichannel_templates,
                               multichannel_stream, cores=cores)
            out["{0}.{1}".format(name, cores)] = cc_out
            if "fftw" in name:
                if cores > 1:
                    print("Running outer core parallel")
                    # Make sure that using both parallel methods gives the same
                    # result
                    cc_out = time_func(
                        func, name, multichannel_templates,
                        multichannel_stream, cores=1, cores_outer=cores)
                    out["{0}.{1}_outer".format(name, cores)] = cc_out
                print("Running fixed fft-len: {0}".format(fft_len))
                # Make sure that running with a pre-defined fft-len works
                cc_out = time_func(
                    func, name, multichannel_templates, multichannel_stream,
                    cores=cores, fft_len=fft_len)
                out["{0}.{1}_fixed_fft".format(name, cores)] = cc_out
    return out


@pytest.fixture(scope='module')
def stream_cc_output_dict(multichannel_templates, multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    return _stream_cc_output_dict(
        multichannel_templates=multichannel_templates,
        multichannel_stream=multichannel_stream
    )


@pytest.fixture(scope='module')
def stream_cc_weighted_output_dict(weighted_multichannel_templates, multichannel_stream):
    return _stream_cc_output_dict(
        multichannel_templates=weighted_multichannel_templates,
        multichannel_stream=multichannel_stream
    )


@pytest.fixture(scope='module')
def stream_cc_dict(stream_cc_output_dict):
    """ return just the cc arrays from the stream_cc functions """
    return {name: result[0] for name, result in stream_cc_output_dict.items()}


@pytest.fixture(scope='module')
def stream_cc_weighted_dict(stream_cc_weighted_output_dict):
    """ return just the cc arrays from the stream_cc functions """
    return {name: result[0] for name, result in stream_cc_weighted_output_dict.items()}


@pytest.fixture(scope='module')
def real_stream_cc_output_dict(real_templates, real_multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    out = {}
    fft_len = next_fast_len(
        real_templates[0][0].stats.npts +
        real_multichannel_stream[0].stats.npts + 1)
    short_fft_len = 2 ** 8
    for name, func in stream_funcs.items():
        for cores in [1, cpu_count()]:
            print("Running {0} with {1} cores".format(name, cores))

            cc_out = time_func(func, name, real_templates,
                               real_multichannel_stream, cores=cores,
                               fft_len=short_fft_len)
            out["{0}.{1}".format(name, cores)] = cc_out
            if "fftw" in name:
                print("Running fixed fft-len: {0}".format(fft_len))
                # Make sure that running with a pre-defined fft-len works
                cc_out = time_func(
                    func, name, real_templates, real_multichannel_stream,
                    cores=cores, fft_len=fft_len)
                out["{0}.{1}_fixed_fft".format(name, cores)] = cc_out
    return out


@pytest.fixture(scope='module')
def real_stream_cc_dict(real_stream_cc_output_dict):
    """ return just the cc arrays from the stream_cc functions """
    return {name: result[0]
            for name, result in real_stream_cc_output_dict.items()}


@pytest.fixture(scope='module')
def gappy_stream_cc_output_dict(
        multichannel_templates, gappy_multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    out = {}
    for name, func in stream_funcs.items():
        for cores in [1, cpu_count()]:
            print("Running {0} with {1} cores".format(name, cores))
            _log_handler.reset()
            cc_out = time_func(func, name, multichannel_templates,
                               gappy_multichannel_stream, cores=cores)
            out["{0}.{1}".format(name, cores)] = (
                cc_out, copy.deepcopy(log_messages['warning']))
            if "fftw" in name:
                if cores > 1:
                    print("Running outer core parallel")
                    _log_handler.reset()
                    # Make sure that using both parallel methods gives the same
                    # result
                    cc_out = time_func(
                        func, name, multichannel_templates,
                        gappy_multichannel_stream, cores=1, cores_outer=cores)
                    out["{0}.{1}_outer".format(name, cores)] = (
                        cc_out, copy.deepcopy(log_messages['warning']))
                print("Running shorter, fixed fft-len")
                # Make sure that running with a pre-defined fft-len works
                cc_out = time_func(
                    func, name, multichannel_templates,
                    gappy_multichannel_stream, cores=cores, fft_len=fft_len)
                out["{0}.{1}_fixed_fft".format(name, cores)] = (
                    cc_out, copy.deepcopy(log_messages['warning']))
    return out


@pytest.fixture(scope='module')
def gappy_stream_cc_dict(gappy_stream_cc_output_dict):
    """ return just the cc arrays from the stream_cc functions """
    return {name: (result[0][0], result[1])
            for name, result in gappy_stream_cc_output_dict.items()}


@pytest.fixture(scope='module')
def gappy_real_cc_output_dict(
        gappy_real_data_template, gappy_real_data):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    out = {}
    for name, func in stream_funcs.items():
        for cores in [1, cpu_count()]:
            _log_handler.reset()
            print("Running {0} with {1} cores".format(name, cores))
            cc_out = time_func(func, name, gappy_real_data_template,
                               gappy_real_data, cores=cores)
            out["{0}.{1}".format(name, cores)] = (
                cc_out, copy.deepcopy(log_messages['warning']))
            if "fftw" in name:
                if cores > 1:
                    _log_handler.reset()
                    print("Running outer core parallel")
                    # Make sure that using both parallel methods gives the same
                    # result
                    cc_out = time_func(
                        func, name, gappy_real_data_template,
                        gappy_real_data, cores=1, cores_outer=cores)
                    out["{0}.{1}_outer".format(name, cores)] = (
                        cc_out, copy.deepcopy(log_messages['warning']))
                _log_handler.reset()
                print("Running shorter, fixed fft-len")
                # Make sure that running with a pre-defined fft-len works
                cc_out = time_func(
                    func, name, gappy_real_data_template, gappy_real_data,
                    cores=cores, fft_len=fft_len)
                out["{0}.{1}_fixed_fft".format(name, cores)] = (
                    cc_out, copy.deepcopy(log_messages['warning']))
    return out


@pytest.fixture(scope='module')
def gappy_real_cc_dict(gappy_real_cc_output_dict):
    """ return just the cc arrays from the stream_cc functions """
    return {name: (result[0][0], result[1])
            for name, result in gappy_real_cc_output_dict.items()}


# ------------------------------------ unstacked setup

@pytest.fixture(scope='module')
def stream_cc_output_dict_unstacked(
        multichannel_templates, multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    for tr in multichannel_stream:
        tr.data = tr.data[0:unstacked_stream_len]
    multichannel_templates = multichannel_templates[0:5]
    out = {}
    for name, func in stream_funcs.items():
        if name.startswith("fmf"):
            print("Skipping fmf - unstacked not implemented")
            continue
        for cores in [1, cpu_count()]:
            print("Running {0} with {1} cores".format(name, cores))

            cc_out = time_func(func, name, multichannel_templates,
                               multichannel_stream, cores=cores, stack=False)
            out["{0}.{1}".format(name, cores)] = cc_out
            if "fftw" in name and cores > 1:
                print("Running outer core parallel")
                # Make sure that using both parallel methods gives the same
                # result
                cc_out = time_func(
                    func, name, multichannel_templates, multichannel_stream,
                    cores=1, cores_outer=cores, stack=False)
                out["{0}.{1}_outer".format(name, cores)] = cc_out
    return out


@pytest.fixture(scope='module')
def stream_cc_dict_unstacked(stream_cc_output_dict_unstacked):
    """ return just the cc arrays from the stream_cc functions """
    return {name: result[0]
            for name, result in stream_cc_output_dict_unstacked.items()}


@pytest.fixture(scope='module')
def gappy_stream_cc_output_dict_unstacked(
        multichannel_templates, gappy_multichannel_stream):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    import warnings

    for tr in gappy_multichannel_stream:
        tr.data = tr.data[0:unstacked_stream_len]
    multichannel_templates = multichannel_templates[0:5]
    out = {}
    for name, func in stream_funcs.items():
        if name.startswith("fmf"):
            print("Skipping fmf - unstacked not implemented")
            continue
        for cores in [1, cpu_count()]:
            # Check for same result both single and multi-threaded
            print("Running {0} with {1} cores".format(name, cores))
            with warnings.catch_warnings(record=True) as w:
                cc_out = time_func(func, name, multichannel_templates,
                                   gappy_multichannel_stream, cores=cores,
                                   stack=False)
                out["{0}.{1}".format(name, cores)] = (cc_out, w)
            if "fftw" in name and cores > 1:
                print("Running outer core parallel")
                # Make sure that using both parallel methods gives the same
                # result
                with warnings.catch_warnings(record=True) as w:
                    cc_out = time_func(
                        func, name, multichannel_templates,
                        gappy_multichannel_stream, cores=1, cores_outer=cores,
                        stack=False)
                out["{0}.{1}_outer".format(name, cores)] = (cc_out, w)
    return out


@pytest.fixture(scope='module')
def gappy_stream_cc_dict_unstacked(gappy_stream_cc_output_dict_unstacked):
    """ return just the cc arrays from the stream_cc functions """
    return {name: (result[0][0], result[1])
            for name, result in gappy_stream_cc_output_dict_unstacked.items()}


@pytest.fixture(scope='module')
def gappy_real_cc_output_dict_unstacked(
        gappy_real_data_template, gappy_real_data):
    """ return a dict of outputs from all stream_xcorr functions """
    # corr._get_array_dicts(multichannel_templates, multichannel_stream)
    import warnings

    for tr in gappy_real_data:
        tr.data = tr.data[0:unstacked_stream_len]
    out = {}
    for name, func in stream_funcs.items():
        if name.startswith("fmf"):
            print("Skipping fmf - unstacked not implemented")
            continue
        for cores in [1, cpu_count()]:
            print("Running {0} with {1} cores".format(name, cores))
            with warnings.catch_warnings(record=True) as w:
                cc_out = time_func(func, name, gappy_real_data_template,
                                   gappy_real_data, cores=cores, stack=False)
                out["{0}.{1}".format(name, cores)] = (cc_out, w)
            if "fftw" in name and cores > 1:
                print("Running outer core parallel")
                # Make sure that using both parallel methods gives the same
                # result
                with warnings.catch_warnings(record=True) as w:
                    cc_out = time_func(
                        func, name, gappy_real_data_template,
                        gappy_real_data, cores=1, cores_outer=cores,
                        stack=False)
                out["{0}.{1}_outer".format(name, cores)] = (cc_out, w)
    return out


@pytest.fixture(scope='module')
def gappy_real_cc_dict_unstacked(gappy_real_cc_output_dict_unstacked):
    """ return just the cc arrays from the stream_cc functions """
    return {name: (result[0][0], result[1])
            for name, result in gappy_real_cc_output_dict_unstacked.items()}

# ----------------------------------- tests


class TestArrayCorrelateFunctions:
    """ these tests ensure the various implementations of normxcorr return
    approximately the same answers """
    # atol = .00001  # how close correlations have to be
    atol = .01

    # tests
    def test_single_channel_similar(self, array_ccs):
        """ ensure each of the correlation methods return similar answers
        given the same input data """
        cc_names = list(array_ccs.keys())
        print(cc_names)
        for key1, key2 in itertools.combinations(cc_names, 2):
            cc1 = array_ccs[key1]
            cc2 = array_ccs[key2]
            if not np.allclose(cc1, cc2, atol=self.atol):
                print("{0} does not match {1}".format(key1, key2))
                np.save("cc1.npy", cc1)
                np.save("cc2.npy", cc2)
            assert np.allclose(cc1, cc2, atol=self.atol)

    def test_known_weight_application(self, array_ccs, array_ccs_weighted):
        cc_names = list(array_ccs.keys())
        print(cc_names)
        for key1, key2 in itertools.combinations(cc_names, 2):
            cc1 = array_ccs_weighted[key1]
            cc2 = array_ccs_weighted[key2]
            print(f"Comparing {key1} and {key2}")
            assert np.allclose(cc1, cc2, atol=self.atol)
        for cc_name in cc_names:
            print(f"Checking for {cc_name}")
            assert np.allclose(
                array_ccs_weighted[cc_name],
                array_ccs[cc_name] * GLOBAL_WEIGHT,
                atol=1e-5)

    def test_autocorrelation(self, array_ccs):
        """ ensure an autocorrelation occurred in each of ccs where it is
        expected, defined by starting_index variable """
        for name, cc in array_ccs.items():
            assert np.isclose(cc[0, starting_index], 1., atol=self.atol)

    def test_non_zero_median(self, array_ccs_low_amp):
        """ Ensure that the median of correlations returned is non-zero,
        this happens with v.0.2.7 when the amplitudes are low."""
        for name, value in array_ccs_low_amp.items():
            cc = value[0]
            warning = value[1]
            assert np.median(cc) != 0.0
            if name == 'fftw':
                assert len(warning) == 1
                assert "Low variance found" in warning[-1]


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

    def test_weighted_multi_channel_xcorr(self, stream_cc_weighted_dict):
        """ test various correlation methods weighted with multiple channels """
        # get correlation results into a list
        cc_names = list(stream_cc_weighted_dict.keys())
        cc_list = [stream_cc_weighted_dict[cc_name] for cc_name in cc_names]
        cc_1 = cc_list[0]
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            assert np.allclose(cc_1, cc, atol=self.atol)

    def test_real_multi_channel_xcorr(self, real_stream_cc_dict):
        """ test various correlation methods with multiple channels """
        # get correlation results into a list
        cc_names = list(real_stream_cc_dict.keys())
        cc_list = [real_stream_cc_dict[cc_name] for cc_name in cc_names]
        cc_1 = cc_list[0]
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            if not np.allclose(cc_1, cc, atol=self.atol):
                print("{0} does not match {1}".format(cc_names[0], cc_name))
                np.save("cc1.npy", cc_1)
                np.save("cc2.npy", cc)
            assert np.allclose(cc_1, cc, atol=self.atol)

    def test_gappy_multi_channel_xcorr(self, gappy_stream_cc_dict):
        """
        test various correlation methods with multiple channels and a gap.
        """
        # get correlation results into a list
        cc_names = list(gappy_stream_cc_dict.keys())
        cc_list = [gappy_stream_cc_dict[cc_name][0] for cc_name in cc_names]
        warning_list = [gappy_stream_cc_dict[cc_name][1]
                        for cc_name in cc_names]
        cc_1 = cc_list[0]
        for cc_name, warning in zip(cc_names, warning_list):
            # fftw_multiprocess doesn't warn?
            if cc_name[0:4] in ['fftw_stream_xcorr', 'fftw_multithread',
                                'fftw_concurrent']:
                assert len(warning) == 1
                assert "are there zeros" in warning[-1]
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            if not np.allclose(cc_1, cc, atol=self.atol * 10):
                print("{0} does not match the master {1}".format(
                    cc_name, cc_names[0]))
                print(np.where((cc_1 - cc) > self.atol * 10))
                np.save("cc.npy", cc)
                np.save("cc_1.npy", cc_1)
                assert np.allclose(cc_1, cc, atol=self.atol * 10)

    def test_gappy_real_multi_channel_xcorr(self, gappy_real_cc_dict):
        """
        test various correlation methods with multiple channels and a gap.

        This test used to fail internally when the variance threshold was too
        low.
        """
        # get correlation results into a list
        cc_names = list(gappy_real_cc_dict.keys())
        cc_list = [gappy_real_cc_dict[cc_name][0] for cc_name in cc_names]
        warning_list = [gappy_real_cc_dict[cc_name][1]
                        for cc_name in cc_names]
        cc_1 = cc_list[0]
        for cc_name, warning in zip(cc_names, warning_list):
            # fftw_multiprocess doesn't warn?
            if cc_name[0:4] in ['fftw_stream_xcorr', 'fftw_multithread',
                                'fftw_concurrent']:
                assert len(warning) == 1
                assert issubclass(warning[-1].category, UserWarning)
                assert "are there zeros" in str(warning[-1].message)
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            if not np.allclose(cc_1, cc, atol=self.atol * 100):
                print("{0} does not match the master {1}".format(
                    cc_name, cc_names[0]))
                print(np.where((cc_1 - cc) > self.atol * 100))
                np.save("cc.npy", cc)
                np.save("cc_1.npy", cc_1)
                assert np.allclose(cc_1, cc, atol=self.atol * 100)


@pytest.mark.serial
class TestStreamCorrelateFunctionsUnstacked:
    """ same thing as TestArrayCorrelateFunction but for stream interface """
    atol = TestArrayCorrelateFunctions.atol

    def test_multi_channel_xcorr(self, stream_cc_dict_unstacked):
        """ test various correlation methods with multiple channels """
        # get correlation results into a list
        cc_names = list(stream_cc_dict_unstacked.keys())
        cc_list = [stream_cc_dict_unstacked[cc_name] for cc_name in cc_names]
        cc_1 = cc_list[0]
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            assert np.allclose(cc_1, cc, atol=self.atol)

    def test_gappy_multi_channel_xcorr(self, gappy_stream_cc_dict_unstacked):
        """
        test various correlation methods with multiple channels and a gap.
        """
        # get correlation results into a list
        cc_names = list(gappy_stream_cc_dict_unstacked.keys())
        cc_list = [gappy_stream_cc_dict_unstacked[cc_name][0]
                   for cc_name in cc_names]
        warning_list = [gappy_stream_cc_dict_unstacked[cc_name][1]
                        for cc_name in cc_names]
        cc_1 = cc_list[0]
        for cc_name, warning in zip(cc_names, warning_list):
            # fftw_multiprocess doesn't warn?
            if cc_name[0:4] in ['fftw_stream_xcorr', 'fftw_multithread',
                                'fftw_concurrent']:
                assert len(warning) == 1
                assert issubclass(warning[-1].category, UserWarning)
                assert "are there zeros" in str(warning[-1].message)
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            if not np.allclose(cc_1, cc, atol=self.atol * 10):
                print("{0} does not match the master {1}".format(
                    cc_name, cc_names[0]))
                print(np.where((cc_1 - cc) > self.atol * 10))
                np.save("cc.npy", cc)
                np.save("cc_1.npy", cc_1)
                assert np.allclose(cc_1, cc, atol=self.atol * 10)

    def test_gappy_real_multi_channel_xcorr(
            self, gappy_real_cc_dict_unstacked):
        """
        test various correlation methods with multiple channels and a gap.

        This test used to fail internally when the variance threshold was too
        low.
        """
        # get correlation results into a list
        cc_names = list(gappy_real_cc_dict_unstacked.keys())
        cc_list = [gappy_real_cc_dict_unstacked[cc_name][0]
                   for cc_name in cc_names]
        warning_list = [gappy_real_cc_dict_unstacked[cc_name][1]
                        for cc_name in cc_names]
        cc_1 = cc_list[0]
        for cc_name, warning in zip(cc_names, warning_list):
            # fftw_multiprocess doesn't warn?
            if cc_name[0:4] in ['fftw_stream_xcorr', 'fftw_multithread',
                                'fftw_concurrent']:
                assert len(warning) == 1
                assert issubclass(warning[-1].category, UserWarning)
                assert "are there zeros" in str(warning[-1].message)
        # loop over correlations and compare each with the first in the list
        # this will ensure all cc are "close enough"
        for cc_name, cc in zip(cc_names[2:], cc_list[2:]):
            if not np.allclose(cc_1, cc, atol=self.atol * 100):
                print("{0} does not match the master {1}".format(
                    cc_name, cc_names[0]))
                print(np.where((cc_1 - cc) > self.atol * 100))
                np.save("cc.npy", cc)
                np.save("cc_1.npy", cc_1)
                assert np.allclose(cc_1, cc, atol=self.atol * 100)


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
        old_default = corr.get_array_xcorr()
        old_default_stream = corr.get_stream_xcorr()
        with corr.set_xcorr('numpy'):
            func = corr.get_array_xcorr()
            assert func is corr.numpy_normxcorr
        assert corr.get_array_xcorr() == old_default
        assert corr.get_stream_xcorr() == old_default_stream


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

        def some_callable(template_array, stream_array, pad_array,
                          weight_array):
            small_count['name'] = 1
            return corr.numpy_normxcorr(template_array, stream_array,
                                        pad_array, weight_array)

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
        # Note: do not remove this fixture or bad things will happen
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
