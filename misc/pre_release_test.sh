#!/usr/bin/env bash

# Pre-release testing script to create a clean virtual environment and try
# running the install and tests in it - designed to work on a linux machine
# with cv2 installed globally (will link to it)

# Needs virtualenv installed as well as cv2

virtualenv -p /usr/bin/python2.7 ~/clean_test
source ~/clean_test/bin/activate
pip install numpy cython
ln -s /usr/local/lib/python2.7/dist-packages/cv2.so ~/clean_test/local/lib/python2.7/.

# Install test runners
pip install pytest pytest-cov pytest-pep8

cd ..
pip install .
python setup.py test

deactivate