#!/usr/bin/python
"""
Utility module of the EQcorrscan package to allow for different methods of
stacking of seismic signal in one place.

In alpha stages and only with linear stacking implimented thusfar

Calum Chamberlain 24/06/2015
"""

def linstack(streams):
    """
    Function to compute the linear stack of a series of seismic streams of
    multiplexed data

    :type streams: List of Streams

    :returns: stack - Stream
    """
    import matplotlib.pyplot as plt
    stack=streams[0]
    for tr in stack:
        tr.data=tr.data/max(tr.data)
    for i in xrange(1,len(streams)):
        print "Stacking stream "+str(i)
        j=0
        for tr in stack:
            print tr.stats.station+'.'+tr.stats.channel
            matchtr=streams[i].select(station=tr.stats.station,\
                                       channel=tr.stats.channel)[0]
            tr.data+=matchtr.data/max(matchtr.data)
            if j==0:
                plt.plot(tr.data)
            j+=1
    plt.show()
    return stack
