#!/usr/bin/python

from multiprocessing import Pool, Array

import matplotlib.pyplot as plt

cccs_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    cccs_list.append(result)

if __name__ == '__main__':
    from obspy import read
    from utils.timer import Timer
    import numpy as np
    from core.match_filter import normxcorr2
    iterations=100
    step=200
    runs=1
    st=read('test_data/HUVA-H1-2015-05-04-processed.ms')
    st.decimate(10)
    # st.decimate(10)
    st=st[0]
    templates=[]
    for i in xrange(iterations):
        stcopy=st.copy()
        templates+=[stcopy.trim(stcopy.stats.starttime+(i*step)+1200,\
                               stcopy.stats.starttime+(i*step)+1800)]
    # for template in templates:
        # print template.stats.starttime
    childnos=np.arange(2,15,1, dtype=int)
    childnos=[20]
    looptimes=[]
    looptime=0
    for j in xrange(runs):
        with Timer() as t:
            # Linear run
            for i in xrange(iterations):
                cccs_list+=[normxcorr2(templates[i].data,st.data, i)[1]]
        looptime+=t.secs
    looptimes+=[looptime/runs]
    for childno in childnos:
        print "Computing for "+str(childno)+" processes"
        looptime=0
        # Run process multiple times and average
        for j in xrange(runs):
            with Timer() as t:
                pool=Pool(processes=childno, maxtasksperchild=None)
                results=[pool.apply_async(normxcorr2, args=(templates[i].data,\
                                                            st.data, i))
                             for i in xrange(iterations)]
                cccs_list=[p.get() for p in results]
                pool.close()
                pool.join()
                cccs_list.sort(key=lambda tup: tup[0]) # Sort by placeholder returned from function
                cccs_list = [ccc[1] for ccc in cccs_list]
            looptime+=t.secs
            print t.secs
        looptimes+=[looptime/runs]
    print len(cccs_list)
    print len(cccs_list[0][0])
    # for ccc in cccs_list:
        # plt.plot(ccc[0])
    # plt.show()
    plt.scatter([0]+childnos,looptimes)
    plt.xlabel("Number of processes in parallel")
    plt.ylabel("Average time after "+str(runs)+" runs")
    plt.show()
