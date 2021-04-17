import os
import shutil
from os.path import join, dirname

import pytest

# ------------------------- Paths constants
# these are added to pytest namespace for convenience if finding things
# eg test data

PKG_PATH = join(dirname(__file__), 'eqcorrscan')
TEST_PATH = join(PKG_PATH, 'tests')
TEST_DATA_PATH = join(TEST_PATH, 'test_data')


# --------------- pytest configuration


def pytest_addoption(parser):
    try:
        parser.addoption("--runslow", action="store_true",
                         help="run slow tests")
    except ValueError as e:
        print(e)
        pass
    try:
        parser.addoption("--runsuperslow", action="store_true",
                         help="run super-slow tests")
    except ValueError as e:
        print(e)
        pass


# ------------------ session fixtures


@pytest.fixture(scope='session', autouse=True)
def clean_up_test_files():
    """ cleanup the file droppings produced by doc tests.
     This fixture will run after all other tests have finished """
    files_to_kill = [
        'test_csv_write.csv',
        'test_family.tgz',
        'test_quakeml.xml',
        'test_tar_write.tgz',
        'test_template.tgz',
        'test_template_read.tgz',
        'test_tribe.tgz',
        'test_waveforms.ms',
        'mag_calc.out',
        'station.dat',
        'station_elev.dat',
        'test_waveform.ms',
        '01-0410-35L.S201309',
        '04-0007-55R.S201601',
        '04-0045-52L.S201601',
        'dt.cc',
        'dt.cc2',
        'dt.ct',
        'dt.ct2',
        'phase.dat'
    ]

    yield

    # remove files
    for fi in files_to_kill:
        if os.path.isfile(fi):
            try:
                os.remove(fi)
            except Exception as e:
                print("File not found, assuming already cleaned")
                print(e)


@pytest.fixture(scope='session', autouse=True)
def clean_up_test_directories():
    """ cleanup the file droppings produced by doc tests.
     This fixture will run after all other tests have finished """
    directories_to_kill = [
        'temp1',
        'temp2',
        'test_family',
        'test_party_out',
        'test_tar_write',
        'tmp1',
        'cc_exported,'
    ]

    yield

    # remove files
    for directory in directories_to_kill:
        if os.path.isdir(directory):
            try:
                shutil.rmtree(directory)
            except Exception as e:
                print("Could not find directory, already cleaned?")
                print(e)


# ------------- add objects to global pytest scope


def append_name(list_like):
    """
    Decorator to append a function name to an object with an append method.

    Useful for making meta-fixtures using the following recipe:
    https://github.com/pytest-dev/pytest/issues/349#issuecomment-189370273

    :param list_like:  Any object with append method
    :return: func
    """

    def _append_func_name(func):
        list_like.append(func.__name__)
        return func

    return _append_func_name


# Stop-gap as per
# https://docs.pytest.org/en/latest/deprecations.html#pytest-namespace
def pytest_configure():
    pytest.append_name = append_name
    pytest.test_path = TEST_PATH
    pytest.test_data_path = TEST_DATA_PATH
    pytest.pkg_path = PKG_PATH


# Over-ride the -n auto in travis as this goes of the wall
# import sys
# import os
#
#
# def pytest_cmdline_preparse(config, args):
#     is_travis = 'TRAVIS' in os.environ
#     if 'xdist' in sys.modules and not is_travis:
#         if "-n" in args:
#             args[args.index("-n") + 1] = "2"
#         else:
#             args[:] = args + ["-n2"]
#     print(args)
