"""
Basic script for running the matched-filter routines in eqcorrscan.
You, the user, are likely to need to add and change things within this file
for your specific use-case, think of this file as providing a basic beginning
for your project.
"""


def run():
    """Internal run function so that this can be called from interactive \
    python session for debugging."""
    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils.archive_read import read_data
    from eqcorrscan.core.match_filter import match_filter
    from obspy import UTCDateTime, Stream
    from eqcorrscan.utils.parameters import read_parameters
    import warnings
    import os
    import datetime as dt
    from obspy import read
    import copy

    # Read parameter files
    par = read_parameters('../parameters/VSP_parameters.txt')
    # Log the input parameters
    log_name = ('EQcorrscan_detection_log_' +
                dt.datetime.now().strftime('%Y.%j.%H:%M:%S') + '.log')
    f = open(os.path.join('..', 'detections', log_name), 'w')
    for parameter in par.__dict__.keys():
        f.write(parameter + ': ' + str(par.__dict__.get(parameter)) + '\n')
    f.write('\n###################################\n')
    f.write('template, detect-time, cccsum, threshold, number of channels\n')

    days = (par.enddate.date - par.startdate.date).days
    dates = [par.startdate + (i * 86400)
             for i in range(days)]

    # Read in templates
    templates = [read(os.path.join('..', 'templates', template))
                 for template in par.template_names]
    # We don't need the full file path in the match-filter routine, just the
    # final 'name'
    template_names_short = [t_name.split(os.sep)[-1]
                            for t_name in par.template_names]
    warnings.warn('Unable to check whether filters are correct in templates')
    # Check that the sampling rate is correct...
    for st in templates:
        for tr in st:
            if not tr.stats.sampling_rate == par.samp_rate:
                msg = 'Template sampling rate is not correct: ' + tr.__str__()
                raise IOError(msg)
    # Work out which stations and channels we will be using
    stachans = [(tr.stats.station, tr.stats.channel)
                for st in templates
                for tr in st]
    stachans = list(set(stachans))
    # Loop through days
    for date in dates:
        # Read in the data
        st = read_data(par.archive, par.arc_type, date.date, stachans)
        # Process the data
        st.merge(fill_value='interpolate')
        st = pre_processing.dayproc(st, lowcut=par.lowcut, highcut=par.highcut,
                                    filt_order=par.filt_order,
                                    samp_rate=par.samp_rate, debug=par.debug,
                                    starttime=UTCDateTime(date.date))
        # Will remove templates if they are deemed useless
        # (eg no matching channels)
        template_names_short_copy = copy.deepcopy(template_names_short)
        templates_copy = copy.deepcopy(templates)
        # Now conduct matched-filter
        detections = match_filter(template_names=template_names_short_copy,
                                  template_list=templates_copy,
                                  st=st, threshold=par.threshold,
                                  threshold_type=par.threshold_type,
                                  trig_int=par.trigger_interval,
                                  plotvar=par.plotvar,
                                  plotdir=par.plotdir,
                                  cores=par.cores,
                                  tempdir=par.tempdir,
                                  debug=par.debug,
                                  plot_format=par.plot_format)
        # Log the output
        for detection in detections:
            f.write(', '.join([detection.template_name,
                               str(detection.detect_time),
                               str(detection.detect_val),
                               str(detection.threshold),
                               str(detection.no_chans)+'\n']))
    f.close()


if __name__ == '__main__':
    run()
