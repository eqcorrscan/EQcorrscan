correlate
---------

.. currentmodule:: eqcorrscan.utils.correlate
.. automodule:: eqcorrscan.utils.correlate

    .. comment to end block

    Classes & Functions
    -------------------
    .. autosummary::
       :toctree: autogen
       :nosignatures:

       fftw_multi_normxcorr
       fftw_normxcorr
       numpy_normxcorr
       time_multi_normxcorr
       get_array_xcorr
       get_stream_xcorr
       register_array_xcorr


    .. comment to end block

All functions within this module (and therefore all correlation functions used
in EQcorrscan) are normalised cross-correlations, which follow the definition
that Matlab uses for |normxcorr2|.

In the time-domain this is as follows
    .. math::
        c(t) = \frac{\sum_{y}(a_{y} - \overline{a})(b_{t+y} - \overline{b_t})}{\sqrt{\sum_{y}(a_{y} - \overline{a})(a_{y} - \overline{a})\sum_{y}(b_{t+y} - \overline{b_t})(b_{t+y} - \overline{b_t})}}

where :math:`c(t)` is the normalised cross-correlation at at time :math:`t`,
:math:`a` is the template window which has length :math:`y`, :math:`b` is the
continuous data, and :math:`\overline{a}` is the mean of the template and
:math:`\overline{b_t}` is the local mean of the continuous data.

In practice the correlations only remove the mean of the template once, but we do
remove the mean from the continuous data at every time-step.  Without doing this
(just removing the global mean), correlations are affected by jumps in average
amplitude, which is most noticeable during aftershock sequences, or when using
short templates.

In the frequency domain (functions :func:`eqcorrscan.utils.correlate.numpy_normxcorr`,
:func:`eqcorrscan.utils.correlate.fftw_normxcorr`, :func:`fftw_multi_normxcorr`),
correlations are computed using an equivalent method.  All methods are tested
against one-another to check for internal consistency, and checked against results
computed using Matlab's |normxcorr2| function. These routines also give the same
(within 0.00001) of openCV cross-correlations.

.. |normxcorr2| raw:: html

    <a href="https://au.mathworks.com/help/images/ref/normxcorr2.html?s_tid=gn_loc_drop" target="_blank">normxcorr2</a>

Selecting a correlation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EQcorrscan strives to use sensible default algorithms for calculating
correlation values, however, you may want to change how correlations are
calculated to be more advantageous to your specific needs.


There are currently 3 different correlations functions currently included in EQcorrscan:

    1. :func:`eqcorrscan.utils.correlate.numpy_normxcorr` known as "numpy"

    2. :func:`eqcorrscan.utils.correlate.time_multi_normxcorr` known as "time_domain"

    3. :func:`eqcorrscan.utils.correlate.fftw_normxcorr` known as "fftw"

Number 3 is the default.

A further time-domain correlation backend is available, see note below on using
Fast Matched Filter within EQcorrscan.

Setting FFT length
~~~~~~~~~~~~~~~~~~
For version 0.4.0 onwards, the "fftw" backend allows the user to pass an `fft_len`
keyword to it.  This will set the length of the FFT internally.  Smaller FFT's use
less memory, but can be faster. Larger FFTs use more memory and can be slower.
By default this is set to the minimum of 2 ** 13, or the next fast length of the
sum of the template length and the stream length. |#285| showed that 2 ** 13 was
consistently fastest over a range of data shapes on an intel i7 with 8-threads.
Powers of two are generally fastest.

.. |#285| raw:: html

    <a href="https://github.com/eqcorrscan/EQcorrscan/pull/285" target="_blank">#285</a>


Using Fast Matched Filter within EQcorrscan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|FMF| provides fast time-domain correlations for both CPU and GPU
architectures.  For massively multi-threaded environment this runs faster than
the frequency-domain routines native to EQcorrscan (when more than 40 CPU cores,
or an NVIDIA GPU card is available) with less memory consumption. |FMF| is now
supported natively, and can be used as a backend by setting `xcorr_func="fmf"`.

.. |FMF| raw:: html

    <a href="https://github.com/beridel/fast_matched_filter" target="_blank">Fast Matched Filter</a>

Using weights
~~~~~~~~~~~~~

New for version 0.5.0: support for weighting cross-correlation stacks. All array
correlation functions have an optional `weights` argument to allow you to weight
different templates. This makes most sense when using the stream correlation functions
where different template channels might have different quality. For the stream
correlation functions the weight for each channel of each template is required to
be stored in `trace.stats.extra.weight`. Individual channel cross-correlations
are multiplied by this weight before stacking. Weights can be any floating point
number, although weights between 0 and 1 make the most sense.


Switching which correlation function is used
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can switch between correlation functions using the `xcorr_func` parameter
included in:
- :func:`eqcorrscan.core.match_filter.match_filter`,
- :meth:`eqcorrscan.core.match_filter.Tribe.detect`,
- :meth:`eqcorrscan.core.match_filter.Template.detect`

by:

1. passing a string (eg "numpy", "time_domain", or "fftw") or;
2. passing a function

for example:

.. code-block:: python

    >>> import obspy
    >>> import numpy as np
    >>> from eqcorrscan.utils.correlate import numpy_normxcorr, set_xcorr
    >>> from eqcorrscan.core.match_filter import match_filter


    >>> # generate some toy templates and stream
    >>> random = np.random.RandomState(42)
    >>> template = obspy.read()
    >>> stream = obspy.read()
    >>> for num, tr in enumerate(stream):  # iter stream and embed templates
    ...     data = tr.data
    ...     tr.data = random.randn(6000) * 5
    ...     tr.data[100: 100 + len(data)] = data

    >>> # do correlation using numpy rather than fftw
    >>> detections = match_filter(['1'], [template], stream, .5, 'absolute',
    ...                           1, False, xcorr_func='numpy')

    >>> # do correlation using a custom function
    >>> def custom_normxcorr(templates, stream, pads, weights *args, **kwargs):
    ...     # Just to keep example short call other xcorr function
    ...     print('calling custom xcorr function')
    ...     return numpy_normxcorr(templates, stream, pads, weights, *args, **kwargs)

    >>> detections = match_filter(
    ...     ['1'], [template], stream, .5, 'absolute', 1, False,
    ...     xcorr_func=custom_normxcorr) # doctest:+ELLIPSIS
    calling custom xcorr function...


You can also use the set_xcorr object (eqcorrscan.utils.correlate.set_xcorr)
to change which correlation function is used. This can be done permanently
or within the scope of a context manager:

.. code-block:: python

    >>> # change the default xcorr function for all code in the with block
    >>> with set_xcorr(custom_normxcorr):
    ...     detections = match_filter(['1'], [template], stream, .5,
    ...                               'absolute', 1, False) # doctest:+ELLIPSIS
    calling custom xcorr function...

    >>> # permanently set the xcorr function (until the python kernel restarts)
    >>> set_xcorr(custom_normxcorr)
    default changed to custom_normxcorr
    >>> detections = match_filter(['1'], [template], stream, .5, 'absolute',
    ...                           1, False) # doctest:+ELLIPSIS
    calling custom xcorr function...
    >>> set_xcorr.revert()  # change it back to the previous state


Notes on accuracy
~~~~~~~~~~~~~~~~~
To cope with floating-point rounding errors, correlations may not be
calculated for data with low variance.  If you see a warning saying:
*"Some correlations not computed, are there zeros in the data? If not
consider increasing gain"*, check whether your data have zeros, and if not,
but have low, but real amplitudes, multiply your data by some value.

Warning
~~~~~~~
If data are padded with zeros prior to filtering then correlations may be
incorrectly calculated where there are zeros. You should always pad after
filtering.  If you see warnings saying *"Low variance, possible zeros, check
result"*, you should check that you have padded correctly and check that
correlations haven't been calculated when you don't expect them.