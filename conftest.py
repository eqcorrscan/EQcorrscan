import os
import shutil

import pytest


# --------------- pytest configuration


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",
                     help="run slow tests")


# ------------------ session fixtures


@pytest.fixture(scope='session', autouse=True)
def clean_up_doctests():
    """ cleanup the file droppings produced by doc tests """
    files_to_kill = [
        'test_csv_write.csv',
        'test_family.tgz',
        'test_quakeml.ml',
        'test_tar_write.tgz',
        'test_template.tgz',
        'test_template_read.tgz',
        'test_tribe.tgz',
        'test_waveforms.ms',
    ]
    directories_to_kill = [

    ]
    yield

    # remove files
    for fi in files_to_kill:
        if os.path.isfile(fi):
            os.remove(fi)

    # remove directories
    for directory in directories_to_kill:
        if os.path.isdir(directory):
            shutil.rmtree(directory)

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
