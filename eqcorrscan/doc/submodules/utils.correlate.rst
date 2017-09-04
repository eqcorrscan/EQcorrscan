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

Selecting a correlation function
--------------------------------

EQcorrscan strives to use sensible default algorithms for calculating
correlation values, however, you may want to change how correlations are
caclulated to be more advantageous to your specific needs.


There are currently 3 different correlations functions currently included in EQcorrscanq:

    1. :func:`eqcorrscan.utils.correlate.numpy_normxcorr` known as "numpy"

    2. :func:`eqcorrscan.utils.correlate.time_multi_normxcorr` known as "time_domain"

    3. :func:`eqcorrscan.utils.correlate.fftw_normxcorr` known as "fftw"

Number 3 is the default.

Switching which correlation function is used
--------------------------------------------

You can switch between correlation functions using the `xcorr_func` paramater
included in:

:func:`eqcorrscan.core.match_filter.match_filter`,
:meth:`eqcorrscan.core.Tribe.detect`,
:meth:`eqcorrscan.core.Template.detect`

by

1. passing a string (eg "numpy", "time_domain", or "fftw")

or

2. passing a function

for example:

.. code-block:: python

    import obspy
    from eqcorrscan.utils.correlate import numpy_normxcorr, set_xcorr
    from eqcorrscan.core.match_filter import match_filter


    # generate some toy templates and stream
    random = np.random.RandomState(42)
    template = obspy.read()
    stream = obspy.read()
    for num, tr in enumerate(stream):  # iter stream and embed templates
        data = tr.data
        tr.data = random.randn(6000) * 5
        tr.data[100: 100 + len(data)] = data


    # do correlation using numpy rather than fftw
    match_filter(['1'], [template], stream, .5, 'absolute', 1, False,
                 xcorr_func='numpy)


    # do correlation using a custom function
    def custom_normxcorr(templates, stream, pads, *args, **kwargs):
        # Just to keep example short call other xcorr function
        print('calling custom xcorr function')
        return numpy_normxcorr(templates, stream, pads, *args, **kwargs)


    match_filter(['1'], [template], stream, .5, 'absolute', 1, False,
                 xcorr_func=custom_normxcorr)
    # prints "calling custom xcorr function


You can also use the set_xcorr object (eqcorrscan.utils.correlate.set_xcorr)
to change which correlation function is used:

.. code-block:: python

    # this can also be done using the set_xcorr function / context manager
    # to change the default xcorr function for all functions that use it
    with set_xcorr(custom_normxcorr):
        match_filter(['1'], [template], stream, .5, 'absolute', 1, False)
        #  prints "calling custom xcorr function"

    # you can also permanently change the xcorr function (until your python
    # kernel is restarted) by calling set_xcorr
    set_xcorr(custom_normxcorr)
    match_filter(['1'], [template], stream, .5, 'absolute', 1, False)
    # prints "calling custom xcorr function
    set_xcorr.revert()  # change it back
