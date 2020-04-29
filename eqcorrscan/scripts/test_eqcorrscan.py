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
PKG_PATH = os.path.dirname(eqcorrscan.__file__)
WORKING_DIR = os.path.join(os.path.expanduser("~"), ".eqcorrscan")
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


def setup_ci():
    if not os.path.isdir(WORKING_DIR):
        os.makedirs(WORKING_DIR)
    for f in ["pytest.ini", "conftest.py"]:
        shutil.copy(
            os.path.join(os.getcwd(), f), WORKING_DIR)
    check_path_conftest(os.path.join(WORKING_DIR, "conftest.py"))
    rewrite_coveragerc(os.path.join(os.getcwd(), ".coveragerc"),
                       os.path.join(WORKING_DIR, ".coveragerc"))
    if os.path.isdir(TEST_DATA_PATH) and \
       len(glob.glob(os.path.join(TEST_DATA_PATH, "*"))) == 0:
        shutil.rmtree(TEST_DATA_PATH)
    if not os.path.isdir(TEST_DATA_PATH):
        test_data_path = os.path.join(
            os.getcwd(), "eqcorrscan", "tests", "test_data")
        shutil.copytree(test_data_path, TEST_DATA_PATH)
    if os.path.isdir(os.path.join(PKG_PATH, "doc")) and \
       len(glob.glob(os.path.join(PKG_PATH, "doc", "*"))) == 0:
        shutil.rmtree(os.path.join(PKG_PATH, "doc"))
    if not os.path.isdir(os.path.join(PKG_PATH, "doc")):
        doc_path = os.path.join(
            os.getcwd(), "eqcorrscan", "doc")
        shutil.copytree(doc_path, os.path.join(PKG_PATH, "doc"))


def download_test_data():
    """Check if test data are installed, and if not, download them"""
    test_data_downloaded = False
    control_files = [
        {"name": "pytest.ini", "downloaded": False},
        {"name": "conftest.py", "downloaded": False},
        {"name": ".coveragerc", "downloaded": False}]
    if os.path.isdir(TEST_DATA_PATH):
        if len(glob.glob(os.path.join(TEST_DATA_PATH, "*"))) > 0:
            Logger.info("Found test data at: {0}".format(TEST_DATA_PATH))
            test_data_downloaded = True
        if os.path.isdir(WORKING_DIR):
            for control_file in control_files:
                working_path = os.path.join(WORKING_DIR, control_file["name"])
                if os.path.isfile(working_path):
                    Logger.info("{0} already here: {1}".format(
                        control_file["name"], working_path))
                    control_file["downloaded"] = True
    else:
        os.makedirs(TEST_DATA_PATH)
    if not os.path.isdir(WORKING_DIR):
        os.makedirs(WORKING_DIR)

    control_files_downloaded = all([c["downloaded"] for c in control_files])
    if test_data_downloaded and control_files_downloaded:
        return

    Logger.info("Downloading necessary files from github")
    with cd(TEST_DATA_PATH):
        # Get the whole zip for this release
        Logger.info("Downloading from {0}".format(TAG_URL))
        r = requests.get(TAG_URL)
        assert r.ok
        contents = zipfile.ZipFile(io.BytesIO(r.content))
        # Files that we want, test-data, conftest.py and pytest.ini
        test_data_path = "EQcorrscan-{0}/eqcorrscan/tests/test_data/".format(
            VERSION)

        for control_file in control_files:
            if not control_file["downloaded"]:
                extract_path = "EQcorrscan-{0}/{1}".format(
                    VERSION, control_file["name"])
                contents.extract(extract_path, ".")
                Logger.debug("Moving {0} to {1}".format(
                    extract_path, WORKING_DIR))
                if control_file["name"] == "conftest.py":
                    check_path_conftest(extract_path)
                if control_file["name"] == ".coveragerc":
                    rewrite_coveragerc(
                        extract_path, os.path.join(WORKING_DIR, ".coveragerc"))
                else:
                    shutil.move(extract_path, WORKING_DIR)

        if not test_data_downloaded:
            for file in contents.namelist():
                if file.startswith(test_data_path):
                    contents.extract(file, '.')
            test_data = glob.glob(
                "EQcorrscan-{0}/eqcorrscan/tests/test_data/*".format(VERSION))
            for test_file in test_data:
                Logger.debug("Moving {0}".format(test_file))
                shutil.move(test_file, ".")
        shutil.rmtree("EQcorrscan-{0}".format(VERSION))
    return


def rewrite_coveragerc(infile, outfile):
    with open(infile, "r") as f:
        contents = f.read().split("\n")
    with open(outfile, "w") as f:
        for line in contents:
            f.write(line.replace("eqcorrscan", PKG_PATH) + "\n")


def check_path_conftest(conftest):
    lines_out = []
    with open(conftest, "r") as f:
        for line in f:
            line.rstrip("\n")
            if line.startswith("PKG_PATH"):
                lines_out.append("PKG_PATH = r'{0}'".format(PKG_PATH))
            else:
                lines_out.append(line)
    with open(conftest, "w") as f:
        for line in lines_out:
            f.write(line + "\n")


def run_tests(arg_list):
    """
    Run the tests.
    """
    # Convert arguments to real paths
    if "--doc" in arg_list:
        doc_files = glob.glob(
            os.path.join(PKG_PATH, "doc", "tutorials", "*.rst"))
        doc_files.extend(glob.glob(
            os.path.join(PKG_PATH, "doc", "submodules", "*.rst")))
        arg_list.extend(doc_files)
        arg_list.remove("--doc")
    elif "--runsuperslow" in arg_list:
        arg_list.append(os.path.join(PKG_PATH, "tests", "tutorials_test.py"))
    else:
        arg_list.append(PKG_PATH)
    if "--runsuperslow" not in arg_list and "--runslow" not in arg_list:
        arg_list.append("--runslow")
    arg_list.extend(
        ["--ignore", "EGG-INFO", "--ignore", PKG_PATH + "/utils/lib",
         "--doctest-modules", "--cov", "--cov-config",
         os.path.join(WORKING_DIR, ".coveragerc"), "--ignore=setup.py"])
    # arg_list.append(PKG_PATH)
    with cd(WORKING_DIR):
        # Copy files to eqcorrscan tree so that pytest can find them
        shutil.copy("pytest.ini", os.path.join(PKG_PATH, "pytest.ini"))
        shutil.copy(".coveragerc", os.path.join(PKG_PATH, ".coveragerc"))
        shutil.copy("conftest.py", os.path.join(PKG_PATH, "conftest.py"))

        Logger.info("Working in {0}".format(WORKING_DIR))
        Logger.info("Running tests from {0}".format(PKG_PATH))
        Logger.info("pytest {0}".format(' '.join(arg_list)))
        # Run the tests!
        ret = pytest.main(args=arg_list)
        # Remove the config files
        if os.path.isfile(os.path.join(PKG_PATH, "conftest.py")):
            os.remove(os.path.join(PKG_PATH, "pytest.ini"))
            os.remove(os.path.join(PKG_PATH, ".coveragerc"))
            os.remove(os.path.join(PKG_PATH, "conftest.py"))
        if ret != 0:
            raise SystemExit("Failed tests")


if __name__ == "__main__":
    arg_list = sys.argv[1:]
    if "--ci" in arg_list:
        setup_ci()
        arg_list.remove("--ci")
    else:
        download_test_data()
    run_tests(arg_list=arg_list)
