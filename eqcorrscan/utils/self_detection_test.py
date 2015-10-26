#!/usr/bin/python
"""
Functions to ape the LFEsearch functionality but for a single template and\
        single day.  This is intended to ensure that the match-filter routine\
        is working as it should be and it can detect the template that has been\
        made!

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

def self_test(template, low_cut, high_cut, filt_order, samp_rate,\
                threshold, thresh_type, trig_int ,debug=0):
    """
    :type template: :class: obspy.Stream
    :param template: Template to check for self-detectability
    :type highcut: float
    :param highcut: High cut in Hz for bandpass
    :type lowcut: float
    :type lowcut: Low cut in Hz for bandpass
    :type filt_order: int
    :param filt_order: Corners for bandpass
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type threshold: float
    :param threshold: A threshold value set based on the threshold_type\
    :type threshold_type: str
    :param threshold_type: The type of threshold to be used, can be MAD,\
        absolute or av_chan_corr.    MAD threshold is calculated as the\
        threshold*(median(abs(cccsum))) where cccsum is the cross-correlation\
        sum for a given template. absolute threhsold is a true absolute\
        threshold based on the cccsum value av_chan_corr is based on the mean\
        values of single-channel cross-correlations assuming all data are\
        present as required for the template, \
        e.g. av_chan_corr_thresh=threshold*(cccsum/len(template)) where\
        template is a single template from the input and the length is the\
        number of channels within this template.
    :type trig_int: float
    :param trig_int: Minimum gap between detections in seconds.
    :type debug: int
    :param debug: Debug output level, higher=more output.
    """
    import sys
    sys.path.append('..')
    import datetime as dt
    import pre_processing
    from eqcorrscan.core import match_filter
    from obspy import read
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
        tr=pre_processing.dayproc(tr, lowcut, highcut, filt_order, samp_rate,\
                                matchdef.debug, date)
    # image.plot(size=(800,600), equal_scale=False)
    # Apply the detection routine with plot on
    detections=match_filter.match_filter(str(template[0].stats.starttime), \
                                         [template], image, threshold,\
                                         threshtype, trig_int,\
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
