Fast Matched Filter
===================

For those of you looking to combine the best of time-domain correlations, GPUs and
EQcorrscan, the developers of Fast Matched Filter have written some very fast, very
parallel correlation routines.  These can be plugged into EQcorrscan to allow you
to use the EQcorrscan front-end and the Fast Matched Filter backend.

For information on the Fast Matched Filter project see their |FMFgithub| page.  If
you do use their codes, please cite the paper they recommend (at the time of writing,
this is |FMFpaper|).

Fast Matched Filter is most useful on massively multi-threaded systems, such as those with more
than 40 CPU cores, or NVIDIA GPU cards (other vendors GPU cards are not supported by Fast
Matched Filter as of 12/06/2018). It is fast, and very light on memory compared to
EQcorrscan's frequency domain calculations. For smaller, CPU-only machines, the frequency domain
wins out though.

Install
-------

Fast Matched Milter is currently only distributed on github.  You can either, follow
the install instructions on the Fast Matched Filter |FMFgithub| page, or, if you
have the pre-requisite compilers (gcc and a recent CUDA C library), install via pip:

.. code-block:: bash

    pip install git+https://github.com/beridel/fast_matched_filter

This should compile the libraries, and install them on your python path along with
the python wrappers.  Fast Matched Filter also comes with Matlab wrappers, but... Matlab.
The install process checks to see if it can compile the CUDA (GPU) version of the
code as well as the CPU version.  If both can be compiled, they will both be available
in the next step.

Usage
-----

To use the Fast Matched Filter correlation routines in EQcorrscan you need to
wrap their functions to provide the expected input and output.  Note that 
Fast Matched Filter allows weighting of the correlations by channel, which
is (as of EQcorrscan v.0.3.0) not supported for native correlation functions. In
this example we set the weights to one for all traces. To wrap Fast Matched Filter
write a function like:

.. code-block:: python

    import numpy as np

    from fastmatchedfilter import fastmatchedfilter

    def fmf_time(templates, stream, pads, args, **kwargs):
        weights = np.ones()
        used_chans = []
        ccc = fastmatchedfilter()
        return ccc, used_chans

This function can then either be passed to any of the matched_filter_ functions
and methods, or set as a default correlation routine as shown in set_correlation_.

.. _matched_filter: core.match_filter.html
.. _set_correlation: utils.correlate.html#switching-which-correlation-function-is-used

.. |FMFgithub| raw:: html

    <a href="https://github.com/beridel/fast_matched_filter" target="_blank">github</a>

.. |FMFpaper| raw:: html

    <a href="https://doi.org/10.1785/0220170181" target="_blank">Beaucé, Eric, W. B. Frank, and Alexey Romanenko (2017). Fast matched-filter (FMF):
     an efficient seismic matched-filter search for both CPU and GPU architectures. Seismological
     Research Letters, doi: 10.1785/0220170181</a>