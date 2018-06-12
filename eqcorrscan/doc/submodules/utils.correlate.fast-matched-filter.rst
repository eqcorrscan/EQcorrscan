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

    from eqcorrscan.utils.correlate import register_array_xcorr
    from fast_matched_filter import matched_filter as fmf


    def flatten_data(templates, stream):
        """
        Helper function to convert to fast-matched-filter's
        require input shapes.
        """
        # Do some reshaping
        stream.sort(['network', 'station', 'location', 'channel'])
        t_starts = []
        for template in templates:
            template.sort(['network', 'station', 'location', 'channel'])
            t_starts.append(min([tr.stats.starttime for tr in template]))
        stations = list(set([tr.stats.station for tr in template
                             for template in templates]))
        channels = {}
        for station in stations:
            channels.update({station: list(set(
                [tr.stats.channel for tr in template.select(station=station)
                 for template in templates]))})
        template_array = []
        stream_array = []
        pad_array = []
        for template in templates:
            sta_array = []
            pads = []
            t_start = template.sort(['starttime'])[0].stats.starttime
            for station in stations:
                chan_array = []
                chan_pad_array = []
                for channel in channels[station]:
                    chan = template.select(station=station, channel=channel)
                    if len(chan) == 0:
                        print("Padding template {0}.{1} with zeros".format(
                              station, channel))
                        chan = template.select(station=station)
                        if len(chan) == 0:
                            chan = template[0].copy()
                        else:
                            chan = chan[0].copy()
                        chan.data = np.zeros(len(chan.data))
                    else:
                        chan = chan[0]
                    chan_array.append(chan.data - np.mean(chan.data))
                    chan_pad_array.append(
                        int((chan.stats.starttime - t_start) *
                            chan.stats.sampling_rate))
                pads.append(chan_pad_array)
                sta_array.append(chan_array)
            template_array.append(sta_array)
            pad_array.append(pads)
        for station in stations:
            chan_array = []
            for channel in channels[station]:
                chan = stream.select(station=station, channel=channel)
                if len(chan) == 0:
                    print("Padding continuous data {0}.{1} with zeros".format(
                        station, channel))
                    chan = stream[0].copy()
                    chan.data = np.zeros(len(chan.data))
                else:
                    chan = chan[0]
                    chan = chan.data
                    chan -= np.mean(chan)
                    chan_array.append(chan)
            stream_array.append(chan_array)
        template_array = np.ascontiguousarray(
            template_array, dtype=np.float32)
        stream_array = np.ascontiguousarray(stream_array, dtype=np.float32)
        pad_array = np.ascontiguousarray(pad_array, dtype=np.float32)
        return template_array, stream_array, pad_array


    def get_weights(templates, stream):
        """
        Get the weighting array - note that your could define
        default weights here for some stations.
        """
        # Use channel_weights to define relative weighting
        # for specific channel ID's
        channel_weights = {'NZ.FOZ.10.HHZ': 0.5}
        # Set a default weight for everything else.
        default_weight = 1.0

        weights = []
        for template in templates:
            template_weights = []
            for station in set([tr.stats.station for tr in stream]):
                station_weights = []
                for tr in stream.select(station=station):
                    try:
                        station_weights.append(channel_weights[tr.id])
                    except KeyError:
                        station_weights.append(default_weight)
                template_weights.append(station_weights)
            weights.append(template_weights)
        return np.ascontiguousarray(weights, dtype=np.float32)


    @register_array_xcorr("fmf")
    def fmf_xcorr(templates, stream, pads, *args, **kwargs):
        """ Reshape arrays and call fast-matched-filter. """
        t_arr, d_arr, pads = flatten_data(templates, stream)
        weights = get_weights(templates, stream)
        cccsums = fmf(
            templates=t_arr, weights=weights, moveouts=pads,
            data=d_arr, step=1, arch='cpu')
        # Use the cpu architecture for single-channels - feel
        # free to change this!
        return cccsums[0]

    @fmf_xcorr.register("concurrent")
    def fmf_multi_xcorr(templates, stream, pads, *args, **kwargs):
        t_arr, d_arr, pads = flatten_data(templates, stream)
        weights = get_weights(templates, stream)
        cccsums = fmf(
            templates=t_arr, weights=weights, moveouts=pads,
            data=d_arr, step=1, arch='cpu')
        # set arch='gpu' if you want to use the gpu and it
        # is available.
        no_chans = []
        chans = []
        for template in templates:
            no_chans.append(len(template))
            chans.append([tr.id for tr in template])
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