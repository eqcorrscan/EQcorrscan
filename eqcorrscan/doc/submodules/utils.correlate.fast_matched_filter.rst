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

    from eqcorrscan.utils.correlate import (
        register_array_xcorr, numpy_normxcorr)

    def flatten_data(templates, stream):
        """
        Helper function to convert to fast-matched-filter's
        require input shapes.
        """
        # Use channel_weights to define relative weighting
        # for specific channel ID's
        # channel_weights = {'NZ.FOZ.10.HHZ': 0.5}
        channel_weights = {}
        # Set a default weight for everything else.
        default_weight = 1.0

        # Ensure stream is zero mean
        stream.detrend()
        # Do some reshaping
        stream.sort(['network', 'station', 'location', 'channel'])
        t_starts = []
        for template in templates:
            template.sort(['network', 'station', 'location', 'channel'])
            t_starts.append(min([tr.stats.starttime for tr in template]))
        stations = list(set([tr.stats.station for tr in template
                             for template in templates]))
        channels = {}
        n_components = 0
        for station in stations:
            n_station_components = 0
            for template in templates:
                unique_station_channels = set(
                    [tr.stats.channel for tr in template.select(station=station)])
                # Need to do this to allow for multiple template waveforms from 
                # the same chanel
                station_channels = []
                for channel in unique_station_channels:
                    _st = template.select(station=station, channel=channel)
                    for i, _tr in enumerate(_st):
                        station_channels.append((_tr.stats.channel, i))
                if len(station_channels) > n_station_components:
                    channels.update({station: station_channels})
                if len(station_channels) > n_components:
                    n_components = len(station_channels)
        template_array = np.zeros((
            len(templates), 
            len(stations),
            n_components, 
            templates[0][0].stats.npts))
        stream_array = np.empty((
            len(stations),
            n_components,
            stream[0].stats.npts))
        pad_array = np.zeros((
            len(templates),
            len(stations),
            n_components))
        weight_array = np.ones_like(pad_array)
        for i, template in enumerate(templates):
            t_start = template.sort(['starttime'])[0].stats.starttime
            for j, station in enumerate(stations):
                for k, channel in enumerate(channels[station]):
                    chan = template.select(
                        station=station, channel=channel[0])[channel[1]]
                    template_array[i][j][k] = chan.data - chan.data.mean()
                    pad_array[i][j][k] = int(
                        (chan.stats.starttime - t_start) *
                        chan.stats.sampling_rate)
                    try:
                        weight = channel_weights[chan.id]
                    except KeyError:
                        weight = default_weight
                    weight_array[i][j][k] = weight
        for j, station in enumerate(stations):
            for k, channel in enumerate(channels[station]):
                chan = stream.select(
                    station=station, channel=channel[0])[0]
                stream_array[j][k] = chan.data - chan.data.mean()
        template_array = np.ascontiguousarray(
            template_array, dtype=np.float32)
        stream_array = np.ascontiguousarray(stream_array, dtype=np.float32)
        pad_array = np.ascontiguousarray(pad_array, dtype=np.int32)
        weight_array = np.ascontiguousarray(weight_array, dtype=np.float32)
        return template_array, stream_array, pad_array, weight_array


    @register_array_xcorr("fmf")
    def fmf_xcorr(templates, stream, pads, *args, **kwargs):
        print("This function is just here as a mapper and does nothing.")
        return numpy_normxcorr(templates, stream, pads, *args, **kwargs)


    @fmf_xcorr.register("stream_xcorr")
    @fmf_xcorr.register("concurrent")
    @fmf_xcorr.register("multiprocess")
    @fmf_xcorr.register("multithread")
    def fmf_multi_xcorr(templates, stream, *args, **kwargs):
        """
        Apply FastMatchedFilter routine concurrently.

        :type templates: list
        :param templates:
            A list of templates, where each one should be an obspy.Stream object
            containing multiple traces of seismic data and the relevant header
            information.
        :type stream: obspy.core.stream.Stream
        :param stream:
            A single Stream object to be correlated with the templates.

        :returns:
            New list of :class:`numpy.ndarray` objects.  These will contain
            the correlation sums for each template for this day of data.
        :rtype: list
        :returns:
            list of ints as number of channels used for each cross-correlation.
        :rtype: list
        :returns:
            list of list of tuples of station, channel for all cross-correlations.
        :rtype: list
        """
        try:
            from fast_matched_filter import matched_filter as fmf
        except ImportError:
            raise ImportError("FastMatchedFilter is not available")

        Logger.info("Flattening data")
        t_arr, d_arr, pads, weights = flatten_data(templates, stream)
        Logger.info("Running cross-correlations")
        cccsums = fmf(
            templates=t_arr, weights=weights, moveouts=pads,
            data=d_arr, step=1, arch='gpu')
        Logger.info("Correlations finished")
        # set arch='gpu' if you want to use the gpu and it
        # is available.
        no_chans = []
        chans = []
        for template in templates:
            no_chans.append(len(template))
            chans.append([tr.id for tr in template])
        return cccsums, no_chans, chans

This function-name ("fmf") can then either be passed to any of the matched_filter_ functions
and methods, or set as a default correlation routine as shown in set_correlation_.

.. _matched_filter: core.match_filter.html
.. _set_correlation: utils.correlate.html#switching-which-correlation-function-is-used

.. |FMFgithub| raw:: html

    <a href="https://github.com/beridel/fast_matched_filter" target="_blank">github</a>

.. |FMFpaper| raw:: html

    <a href="https://doi.org/10.1785/0220170181" target="_blank">Beaucé, Eric, W. B. Frank, and Alexey Romanenko (2017). Fast matched-filter (FMF):
     an efficient seismic matched-filter search for both CPU and GPU architectures. Seismological
     Research Letters, doi: 10.1785/0220170181</a>
