#!/usr/bin/env python
"""
Collect and run the EQcorrscan tests.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import os
import logging
import sys
import pytest
import requests
import glob
import zipfile
import io
import shutil

import eqcorrscan
from eqcorrscan import tests

logging.basicConfig(
    level=logging.INFO, stream=sys.stdout,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")
Logger = logging.getLogger(__name__)

VERSION = eqcorrscan.__version__
TEST_PATH = os.path.dirname(tests.__file__)
PKG_PATH = os.path.dirname(os.path.dirname(eqcorrscan.__file__))
TAG_URL = "https://github.com/eqcorrscan/EQcorrscan/archive/{0}.zip".format(
    VERSION)
TEST_DATA_PATH = os.path.join(TEST_PATH, 'test_data')


class cd:
    """
    Context manager for changing the current working directory.

    From: https://stackoverflow.com/questions/431684/\
    how-do-i-change-directory-cd-in-python/13197763#13197763
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def download_test_data():
    """Check if test data are installed, and if not, download them"""
    if os.path.isdir(TEST_DATA_PATH):
        if len(glob.glob(TEST_DATA_PATH)) > 0:
            Logger.info("Found test data at: {0}".format(TEST_DATA_PATH))
            return
    Logger.info("Downloading test data from github")
    os.makedirs(TEST_DATA_PATH)
    with cd(TEST_DATA_PATH):
        # Get the whole zip for this release
        Logger.info("Downloading from {0}".format(TAG_URL))
        r = requests.get(TAG_URL)
        assert r.ok
        contents = zipfile.ZipFile(io.BytesIO(r.content))
        # Files that we want, test-data, conftest.py and pytest.ini
        test_data_path = "EQcorrscan-{0}/eqcorrscan/tests/test_data/".format(
            VERSION)
        conftest = "EQcorrscan-{0}/conftest.py".format(VERSION)
        contents.extract(conftest, '.')
        pytestini = "EQcorrscan-{0}/pytest.ini".format(VERSION)
        contents.extract(pytestini, '.')
        for file in contents.namelist():
            if file.startswith(test_data_path):
                contents.extract(file, '.')
        test_data = glob.glob(
            "EQcorrscan-{0}/eqcorrscan/tests/test_data/*".format(VERSION))
        for test_file in test_data:
            Logger.debug("Moving {0}".format(test_file))
            shutil.move(test_file, ".")
        Logger.debug("Moving {0} to {1}".format(conftest, PKG_PATH))
        shutil.move(conftest, PKG_PATH)
        Logger.debug("Moving {0} to {1}".format(pytestini, PKG_PATH))
        shutil.move(pytestini, PKG_PATH)
        shutil.rmtree("EQcorrscan-{0}".format(VERSION))
    return


def run_tests(arg_list):
    """
    Run the tests.
    """
    arg_list.extend(["--doctest-modules", "--runslow"])
    arg_list.extend(
        ["--ignore", "EGG-INFO", "--ignore", "eqcorrscan/utils/lib"])
    with cd(PKG_PATH):
        pytest.main(args=arg_list)


if __name__ == "__main__":
    arg_list = sys.argv[1:]
    download_test_data()
    run_tests(arg_list=arg_list)
