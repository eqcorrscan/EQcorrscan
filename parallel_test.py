#!/usr/bin/python
"""
Test of parallel speed ups
"""
from math import sqrt
from joblib import Parallel, delayed
import time

import numpy as np
from obspy import read
import pylab as plt

print "First i will run a simple example using sqrt from math, \
        this will be slower in parallel"
# Non parallel process
tic=time.clock()
n=[sqrt(i ** 2) for i in range(1000)]
toc=time.clock()
lintime=toc-tic

# Parallel equivilant
tic=time.clock()
n=Parallel(n_jobs=10, verbose=0)(delayed(sqrt)(i) for i in range(1000))
toc=time.clock()
partime=toc-tic

print 'Time taken in linear run: '
print lintime
print 'Time taken in parallel run: '
print partime

print "Now I will run numpy's correlate function both linearly and in parallel,\
        this should be faster in parallel"

st=read('/media/Emilys_Backup/COSA_processing/day_volumes_S/Y2014/R297.01/TEPE.CO..HHZ.2014.297')
image=st[0].data
image=(image-np.mean(image))/(np.std(image)*len(image))
templates=[]
for i in range(20):
    template=st.slice(st[0].stats.starttime+100*i,\
                        st[0].stats.starttime+(100*i)+10)[0].data
    templates.append((template-np.mean(template))/np.std(template))

# Linear run
ccclin=[]
tic=time.clock()
for i in range(len(templates)):
    ccclin.append(np.correlate(image,templates[i],mode='full'))
toc=time.clock()
lintime=toc-tic

# for i in range(len(ccclin)):
    # plt.plot(ccclin[i])
    # plt.show()


# Parallel run
tic=time.clock()
cccpar=Parallel(n_jobs=10, verbose=0)(delayed(np.correlate)\
                                      (image,templates[i],mode='full') \
                                      for i in range(len(templates)))
toc=time.clock()
partime=toc-tic

# for i in range(len(cccpar)):
    # plt.plot(cccpar[i])
    # plt.show()

print 'Time taken in linear run: '
print lintime
print 'Time taken in parallel run: '
print partime
