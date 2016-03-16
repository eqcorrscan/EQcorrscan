"""
Basic script for running the matched-filter routines in eqcorrscan.
You, the user, are likely to need to add and change things within this file
for your specific use-case, think of this file as providing a basic beginning
for your project.
"""

if __name__ == '__main__':
    from eqcorrscan.utils import pre_processing
    from eqcorrscan.utils.archive_read import read_data
    from eqcorrscan.core.match_filter import match_filter
    from obspy import UTCDateTime, Stream
    from eqcorrscan.utils.parameters import read_parameters
    import warnings

    # Read parameter files
    par = read_parameters('../parameters/EQcorrscan_parameters.txt')
    # Log the input parameters
    f = open(os.path.join('..', 'detections', log_name), 'w')
    for parameter in par.__dict__.keys():
        f.write(parameter + ': ' + par.__dict__.get(parameter) + '\n')
    f.write('\n###################################\n')
    f.write('template, detect-time, cccsum, threshold, number of channels\n')

    days = (par.enddate.date - par.startdate.date).days
    dates = [par.startdate + (i * 86400)
             for i in range(days)]

    # Read in templates
    templates = [read(os.path.join('..', 'templates', template))
                 for template in par.template_names]
    warnings.warn('Unable to check whether filters are correct in templates')
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
        st = Stream(st)

        # Now conduct matched-filter
        detections = match_filter(template_names=par.template_names,
                                  template_list=templates,
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
        f.write(','.join([detection.template_name, str(detection.detect_time),
                          str(detection.detect_val), str(detection.threshold),
                          str(detection.no_chans)+'\n']))
    f.close()
