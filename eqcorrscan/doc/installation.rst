EQcorrscan installation
=======================

EQcorrscan is a Python package with C extensions. The C extensions in EQcorrscan
have their own dependencies on compiled libraries. We heavily recommend installing
EQcorrscan using conda because this will:

* make your life easier;
* separate your EQcorrscan install from your system Python, meaning you can
  experiment to your hearts-content without breaking your operating system (yay);
* ensure that compiled modules are compiled using the correct C-compiler against
  the correct libraries


If you do not have either a miniconda or anaconda installation you can follow
the |conda-install| instructions.

If you do not already have a conda environment we recommend creating one
with the following:

.. code-block:: bash

    conda create -n eqcorrscan -c conda-forge colorama numpy scipy matplotlib obspy bottleneck pyproj python=3.8
    source activate eqcorrscan

To then install EQcorrscan you can simply run:

.. code-block:: bash

    conda install -c conda-forge eqcorrscan

Installation without conda
--------------------------

Installing EQcorrscan without conda involves two steps:

1. Installing fftw3 libraries;
2. Installing python dependancies and EQcorrscan.


How you undertake the first step depends on your operating system and system
package manager.

Non-Python dependencies--Ubuntu:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prior to installing the python routines you will need to install the fftw
library.  On linux use apt (or your default package manager - note you may need
sudo access):

.. code-block:: bash

    apt-get install libfftw3-dev

Note that you will need to ensure you have the single-precision libraries of
fftw3 (files named fftw3f...). On CentOS you can install the `fftw-libs` package.

Non-Python dependencies--OSX:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For MacOS/OS-X systems we have tested using homebrew and macports (fink options
are available, but we haven't tested them).

Homebrew
........

You will need a recent version of gcc (the homebrew gcc-4.9 port has issues with openMP).
We have tested the following and found it to work (note that you may need to prepend
sudo depending on your configuration):

.. code-block:: bash

    brew install gcc6
    brew install fftw

Then run the following to install EQcorrscan (note the need to select CC=gcc, you can
install using clang, but you will need additional libraries for openmp support):

.. code-block:: bash

    CC=gcc pip install eqcorrscan


MacPorts
........

The following has been tested and found to work (note that you may need to prepend
sudo depending on your configuration):

1. Install an up-to-date gcc (gcc is needed for openmp compatibility) - any gcc should work (>4), here we use gcc6 for example:

.. code-block:: bash

    port install gcc6

2. Install python from macports (tested for python35, but its up to you)

.. code-block:: bash

    port install python35`
    # optional: select python35 as default python for terminal:
    port select --set python python35

3. Install numpy and pip from macports:

.. code-block:: bash

    port install py35-numpy py35-pip
    # optional, select pip35 as default pip
    port select --set pip pip35

4. Install fftw3 from source:

    a. |fftw-3.3.7| - link to fftw 3.3.7, most recent as of 10/01/2018
    b. unzip/untar
    c. Run the following from within the expanded directory:

    .. code-block:: bash

        ./configure --enable-threads --enable-float && make
        make install
        ./configure --enable-threads && make # Need both double and float precision files
        make install

5. Run: (if you didn't run the `port select --set pip pip35` command you will need to replace `pip` with `pip35`)

.. code-block:: bash

    CC=gcc pip install eqcorrscan


Non-Python dependencies--Windows:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Windows systems you should follow the instructions on the |fftw-windows|
page and use the pre-compiled dynamic libraries. These should be installed
somewhere on your system path, or the install location added to your path.
The correlation routines use openMP for parallel workflows, however, some aspects
of this run into issues with version of MSVC < 10.0 (due to old C standards being
used), as such, by default, the correlation routines are compiled as serial
workflows on windows.  If you have a need for this threading in windows please
get in touch with the developers.

EQcorrscan install via pip:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have installed fftw the EQcorrscan install should be as simple as:

.. code-block:: bash

    pip install eqcorrscan

.. |conda-install| raw:: html

    <a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">conda-install</a>


.. |fftw-install| raw:: html

    <a href="http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix" target="_blank">fftw installation</a>

.. |fftw-3.3.7| raw:: html

    <a href="http://www.fftw.org/fftw-3.3.7.tar.gz" target="_blank">Download</a>

.. |fftw-windows| raw:: html

    <a href="http://www.fftw.org/install/windows.html" target="_blank">fftw-windows install</a>

Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

pip pulls the package from the |pypi| package repository and runs the `setup.py` file.
If instead you wish to install from source, download the package (either by cloning
the git repository, or by downloading the source code) from |eqcorrscan-github|,
change directory to the `EQcorrscan` directory and run:

.. code-block:: bash

    python setup.py install

If this fails because the default compiler is `clang` you can run:

.. code-block:: bash

    CC=gcc python setup.py install

Note though that this will compile EQcorrscan using a different compiler than
used to build your Python, which may have unwanted effects, if you do this you
MUST test you install using the instructions here: :ref:`RunningTests`.


.. |pypi| raw:: html

    <a href="https://pypi.org/project/EQcorrscan/" target="_blank">PyPi</a>

.. |eqcorrscan-github| raw:: html

    <a href="https://github.com/eqcorrscan/EQcorrscan" target="_blank">github</a>

Using Intel's MKL
~~~~~~~~~~~~~~~~~

For versions >= 0.3.0 EQcorrscan supports compilation against the Intel Math Kernel
Libraries (MKL). This has shown |speed-ups| compared to the standard FFTW library.
To enable this you must install MKL before compiling EQcorrscan.  MKL is available from
most package managers (including conda). Once you have MKL installed you can
follow the `Installation from source`_ section.  Check that near the top of the
install that the MKL libraries are found.


Notes
-----

You may have issues with these installs if you don't have numpy installed: but if
you don't have numpy installed then you have bigger issues...

If you plan to generate a grid of synthetic templates you will need to have
grid csv files, which the authors have previously used NonLinLoc to generate.
This is not provided here and should be sourced from |NLLoc_link|. This will provide
the Grid2Time routine which is required to set-up a lag-time grid for your
velocity model.  You should read the NonLinLoc documentation for more
information regarding how this process works and the input files you are
required to give.

.. |NLLoc_link| raw:: html

  <a href="http://alomax.free.fr/nlloc/" target="_blank">NonLinLoc</a>

.. |speed-ups| raw:: html

  <a href="https://github.com/eqcorrscan/EQcorrscan/pull/168" target="_blank">speed ups</a>
