#!/usr/bin/env bash

cd ~

git clone --depth=50 --branch=develop https://github.com/calum-chamberlain/EQcorrscan.git calum-chamberlain/EQcorrscan
cd calum-chamberlain/EQcorrscan

export OBSPY_VERSION=1.0.1
source ~/virtualenv/python3.5/bin/activate

python --version
pip --version

export OS="Linux"; export py="3.5"

wget https://repo.continuum.io/miniconda/Miniconda3-latest-${OS}-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda

export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda

conda info -a

if [[ "$py" == "3.3" ]]; then
    NUMPY_VERSION=1.9.2
    SCIPY_VERSION=0.16.0
    MPL_VERSION=1.4.3
    BASEMAP_VERSION=1.0.7
    PYPROJ=""
    PYTEST="pytest>=2.8"
    PYFLAKES="pyflakes=0.9.0"
elif [[ "$py" == "3.5" ]]; then
    NUMPY_VERSION=1.9.2
    SCIPY_VERSION=0.16.0
    MPL_VERSION=1.4.3
    BASEMAP_VERSION=1.0.7
    PYPROJ=""
    PYTEST="pytest>=2.8"
    PYFLAKES="pyflakes=1.0.0"
elif [[ "${py:0:1}" == '3' ]]; then
    NUMPY_VERSION=1.10.4
    SCIPY_VERSION=0.17.0
    MPL_VERSION=1.5.1
    BASEMAP_VERSION=1.0.7
    PYPROJ="pyproj"
    PYTEST=""
    PYFLAKES="pyflakes=0.9.0"
else
    NUMPY_VERSION=1.10.4
    SCIPY_VERSION=0.17.0
    MPL_VERSION=1.5.1
    BASEMAP_VERSION=1.0.7
    PYPROJ="pyproj"
    PYTEST=""
    PYFLAKES="pyflakes=0.9.0"
fi

if [[ "${py:0:1}" == '3' ]]; then
  OPENCV="-c menpo opencv3=3.1.0"
else
  OPENCV="opencv"
fi

conda create -n test-environment python=$py numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION matplotlib=$MPL_VERSION basemap=$BASEMAP_VERSION $PYPROJ flake8 future lxml decorator sqlalchemy mock nose gdal docopt coverage requests
source activate test-environment

conda install $PYFLAKES
conda install $OPENCV
pip install coveralls
pip install geographiclib
pip install https://github.com/megies/PyImgur/archive/py3.zip
pip install pep8-naming
pip install pytest
pip install pytest-cov
pip install obspy==$OBSPY_VERSION
pip freeze
conda list
git version
pip install .

python setup.py test

source deactivate
conda remove --name test-environment --all

rm /home/calumch/miniconda/bin/conda

cd ~
rm -r calum-chamberlain
rm -r $HOME/miniconda