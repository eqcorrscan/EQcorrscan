r"""Functions to locate earthquakes detected by EQcorrscan - designed first to\
locate stacks of detections to give family locations.  Extensions may later be\
written.

Copyright 2015 Calum Chamberlain

This file is part of EQcorrscan.

    EQcorrscan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EQcorrscan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EQcorrscan.  If not, see <http://www.gnu.org/licenses/>.
"""


def synth_compare(stream, stream_list, cores=4, debug=0):
    r"""Compare a specific stream to a list of synthetic templates, or \
    earthquakes of known source and find the best matching event.

    :type stream: :class:obspy.Stream
    :param stream: Stream to be compared to streams with known locations.
    :type stream_list: list
    :param stream_list: List of streams with known locations
    :type cores: int
    :param cores: Number of cores to parallel over
    :type debug: int
    :param debug: Debug level, high is more debug

    :returns: int, float: index of best match and cross-correlation sum
    """

    from eqcorrscan.core.match_filter import _channel_loop
    import numpy as np
    import copy
    from obspy import Trace

    stream_copy = stream.copy()
    templates = copy.deepcopy(stream_list)
    # Need to fill the stream_list - template - channels
    template_stachan = []
    for template in templates:
        for tr in template:
            template_stachan += [(tr.stats.station, tr.stats.channel)]
    template_stachan = list(set(template_stachan))

    for stachan in template_stachan:
        if not stream_copy.select(station=stachan[0], channel=stachan[1]):
            # Remove template traces rather than adding NaN data
            for template in templates:
                if template.select(station=stachan[0], channel=stachan[1]):
                    for tr in template.select(station=stachan[0],
                                              channel=stachan[1]):
                        template.remove(tr)
    # Remove un-needed channels
    for tr in stream_copy:
        if not (tr.stats.station, tr.stats.channel) in template_stachan:
            stream_copy.remove(tr)
    # Also pad out templates to have all channels
    for template in templates:
        for stachan in template_stachan:
            if not template.select(station=stachan[0], channel=stachan[1]):
                nulltrace = Trace()
                nulltrace.stats.station = stachan[0]
                nulltrace.stats.channel = stachan[1]
                nulltrace.stats.sampling_rate = template[0].stats.sampling_rate
                nulltrace.stats.starttime = template[0].stats.starttime
                nulltrace.data = np.array([np.NaN] * len(template[0].data),
                                          dtype=np.float32)
                template += nulltrace
    # Hand off  cross-correaltion to _channel_loop, which runs in parallel
    [cccsums, no_chans] = _channel_loop(templates, stream_copy, cores, debug)
    cccsums = [np.max(cccsum) for cccsum in cccsums]
    # Find the maximum cccsum and index thereof
    index = np.argmax(cccsums)
    cccsum = cccsums[index]
    return index, cccsum
