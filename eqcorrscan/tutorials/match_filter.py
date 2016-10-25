"""
Simple tutorial to demonstrate some of the basic capabilities of the EQcorrscan
matched-filter detection routine.  This builds on the template generation
tutorial and uses those templates.  If you haven't run that tutorial script
then you will need to before you can run this script.
"""

def run_tutorial(plot=False):
    """Main function to run the tutorial dataset."""

    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils import plotting
    from eqcorrscan.core import match_filter
    import glob
    from multiprocessing import cpu_count

    # This import section copes with namespace changes between obspy versions
    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client
    from obspy import UTCDateTime, Stream, read

    # First we want to load our templates
    template_names = glob.glob('tutorial_template_*.ms')

    if len(template_names) == 0:
        raise IOError('Template files not found, have you run the template ' +
                      'creation tutorial?')

    templates = [read(template_name) for template_name in template_names]

    # Work out what stations we have and get the data for them
    stations = []
    for template in templates:
        for tr in template:
            stations.append((tr.stats.station, tr.stats.channel))
    # Get a unique list of stations
    stations = list(set(stations))

    # We will loop through the data chunks at a time, these chunks can be any
    # size, in general we have used 1 day as our standard, but this can be
    # as short as five minutes (for MAD thresholds) or shorter for other
    # threshold metrics. However the chunk size should be the same as your
    # template process_len.

    # You should test different parameters!!!
    start_time = UTCDateTime(2016, 1, 4)
    end_time = UTCDateTime(2016, 1, 5)
    process_len = 3600
    chunks = []
    chunk_start = start_time
    while chunk_start < end_time:
        chunk_end = chunk_start + process_len
        if chunk_end > end_time:
            chunk_end = end_time
        chunks.append((chunk_start, chunk_end))
        chunk_start += process_len

    unique_detections = []
    detections = []

    # Set up a client to access the GeoNet database
    client = Client("GEONET")

    # Note that these chunks do not rely on each other, and could be paralleled
    # on multiple nodes of a distributed cluster, see the SLURM tutorial for
    # an example of this.
    for t1, t2 in chunks:
        # Generate the bulk information to query the GeoNet database
        bulk_info = []
        for station in stations:
            bulk_info.append(('NZ', station[0], '*',
                              station[1][0] + 'H' + station[1][-1], t1, t2))

        # Note this will take a little while.
        print('Downloading seismic data, this may take a while')
        st = client.get_waveforms_bulk(bulk_info)
        # Merge the stream, it will be downloaded in chunks
        st.merge(fill_value='interpolate')

        # Set how many cores we want to parallel across, we will set this to four
        # as this is the number of templates, if your machine has fewer than four
        # cores/CPUs the multiprocessing will wait until there is a free core.
        # Setting this to be higher than the number of templates will have no
        # increase in speed as only detections for each template are computed in
        # parallel.  It may also slow your processing by using more memory than
        # needed, to the extent that swap may be filled.
        if cpu_count() < 4:
            ncores = cpu_count()
        else:
            ncores = 4

        # Pre-process the data to set frequency band and sampling rate
        # Note that this is, and MUST BE the same as the parameters used for the
        # template creation.
        print('Processing the seismic data')
        st = pre_processing.shortproc(st, lowcut=2.0, highcut=9.0,
                                      filt_order=4, samp_rate=20.0,
                                      debug=2, num_cores=ncores, starttime=t1,
                                      endtime=t2)
        # Convert from list to stream
        st = Stream(st)

        # Now we can conduct the matched-filter detection
        detections += match_filter.match_filter(template_names=template_names,
                                                template_list=templates,
                                                st=st, threshold=8.0,
                                                threshold_type='MAD',
                                                trig_int=6.0, plotvar=plot,
                                                plotdir='.', cores=ncores,
                                                debug=1, plot_format='jpg')

    # Now lets try and work out how many unique events we have just to compare
    # with the GeoNet catalog of 20 events on this day in this sequence
    for master in detections:
        keep = True
        for slave in detections:
            if not master == slave and\
               abs(master.detect_time - slave.detect_time) <= 1.0:
                # If the events are within 1s of each other then test which
                # was the 'best' match, strongest detection
                if not master.detect_val > slave.detect_val:
                    keep = False
                    break
        if keep:
            unique_detections.append(master)

    print('We made a total of ' + str(len(unique_detections)) + ' detections')

    for detection in unique_detections:
        print('Detection at :' + str(detection.detect_time) +
              ' for template ' + detection.template_name +
              ' with a cross-correlation sum of: ' +
              str(detection.detect_val))
        # We can plot these too
        if plot:
            stplot = st.copy()
            template = templates[template_names.index(detection.template_name)]
            lags = sorted([tr.stats.starttime for tr in template])
            maxlag = lags[-1] - lags[0]
            stplot.trim(starttime=detection.detect_time - 10,
                        endtime=detection.detect_time + maxlag + 10)
            plotting.detection_multiplot(stplot, template,
                                         [detection.detect_time.datetime])
    return unique_detections

if __name__ == '__main__':
    run_tutorial()
