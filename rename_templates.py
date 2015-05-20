#!/usr/bin/python
def chan_rename(tr):
    if tr.stats.channel=='SHZ':
        tr.stats.channel='SH1'
    elif tr.stats.channel=='SHN':
        tr.stats.channel='SH2'
    elif tr.stats.channel=='SHE':
        tr.stats.channel='SH3'
    return tr

def chan_rename_2(tr):
    if tr.stats.channel=='SHZ':
        tr.stats.channel='SHZ'
    elif tr.stats.channel=='SH1':
        tr.stats.channel='SHN'
    elif tr.stats.channel=='SH2':
        tr.stats.channel='SHE'
    return tr

sfiles=['63877_1_2_matlab',
 '17208_1_2_3_matlab',
 '60070_1_2_matlab',
 '37575_1_2_matlab',
 '55115_1_1_matlab',
 '60905_1_matlab',
 '61044_1_2_3_matlab',
 '54905_1_2_3_4_5_matlab',
 '61220_1_2_3_matlab',
 '30442_1_2_matlab',
 '55432_1_2_3_4_matlab',
 '55200_1_2_3_4_5_matlab',
 '61100_1_matlab',
 '59966_1_2_matlab']
templates=[]
from obspy import read, Stream
for sfile in sfiles:
    templates+=[read('templates/'+sfile+'_template.ms')]
i=0
for template in templates:
    template_new=Stream([])
    for tr in template:
        if tr.stats.station=='POCR':
            tr.stats.station='POCR2'
            tr=chan_rename_2(tr)
            tr=chan_rename(tr)
        elif tr.stats.station=='WHAT':
            tr.stats.station='WHAT2'
            tr=chan_rename_2(tr)
            tr=chan_rename(tr)
        elif tr.stats.station=='FRAN':
            tr=chan_rename_2(tr)
            tr=chan_rename(tr)
        elif tr.stats.station=='FOZ':
            tr.stats.network='NZ'
    template.write('templates/'+sfiles[i]+'_renamed_template.ms', format='MSEED')
    i+=1
