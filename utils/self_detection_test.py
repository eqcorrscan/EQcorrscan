#!/usr/bin/python
"""
Functions to ape the LFEsearch functionality but for a single template and\
        single day.  This is intended to ensure that the match-filter routine\
        is working as it should be and it can detect the template that has been\
        made!
"""

def self_test(template):
    """
    :type template: :class: obspy.Stream
    :param template: Template to check for self-detectability
    """
    import sys
    sys.path.append('..')
    import datetime as dt
    import pre_processing
    from EQcorrscan.core import match_filter
    from obspy import read
    from EQcorrscan.par import match_filter_par as matchdef
    from EQcorrscan.par import template_gen_par as templatedef
    # Work out the date of the template
    date=template[0].stats.starttime.datetime.date()
    # Read in the appropriate data
    sta_chans=[(tr.stats.station,tr.stats.channel,tr.stats.network) for tr in template]
    for sta_chan in sta_chans:
        base=matchdef.contbase[[i for i in xrange(len(matchdef.contbase)) \
                                if matchdef.contbase[i][2]==sta_chan[2]][0]]
        if base[1]=='yyyymmdd':
            daydir=date.strftime('%Y%m%d')
            staform='*'+sta_chan[0]+'.'+sta_chan[1][0]+'*'+sta_chan[1][1]+'.*'
        elif base[1]=='Yyyyy/Rjjj.01':
            daydir=date.strftime('Y%Y/R%j.01')
            staform=sta_chan[0]+'.*.'+sta_chan[1][0]+'*'+sta_chan[1][1]+\
                    '.'+date.strftime('%Y.%j')
        else:
            raise IOError('Not in the correct form'+base[1])
        if not 'image' in locals():
            image=read(base[0]+'/'+daydir+'/'+staform)
        else:
            image+=read(base[0]+'/'+daydir+'/'+staform)
    # Process the data using pre-processing
    for tr in image:
        tr=pre_processing.dayproc(tr, templatedef.lowcut, templatedef.highcut,\
                                templatedef.filter_order, templatedef.samp_rate,\
                                matchdef.debug, date)
    # image.plot(size=(800,600), equal_scale=False)
    # Apply the detection routine with plot on
    detections=match_filter.match_filter(str(template[0].stats.starttime), \
                                         [template], image, matchdef.threshold,\
                                         matchdef.threshtype, matchdef.trig_int,\
                                         True)
    for detection in detections:
        print 'Detection using template: '+detection.template_name+' at '+\
                str(detection.detect_time)+' with a cccsum of: '+\
                str(detection.detect_val)


if __name__ == '__main__':
    import sys
    if not len(sys.argv)==2:
        raise IOError('Needs one input, the name of the template to use')
    template=str(sys.argv[1])
    from obspy import read
    template=read(template)
    self_test(template)
