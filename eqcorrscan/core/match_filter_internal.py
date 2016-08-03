"""
Functions written to be compilled by Cython as the inner loops of \
the match_filter.py routine.

:copyright:
    Calum Chamberlain, Chet Hopp.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np


def _channel_loop(templates, stream, delays, ktemplates, savedir=False,
                  cores=10):
    """
    Loop to generate cross channel correaltion sums for a series of templates
    hands off the actual correlations to a sister function which can be run in
    parallel.

    :type templates: .ndarray
    :param templates: A series of templates, organized into one np.ndarray \
        to the extent that the returned ccc[i][:] will correspond to \
        template[i][:]
    :type stream: np.ndarray
    :param stream: Image stream - np.ndarray, should be daylong
    :type delays: np.ndarray
    :param delays: Delays for each template in the templates array - must be \
            in samples, these will be the length of the pads applied to the \
            stream arrays
    :type ktemplates: int
    :param ktemplates: The number of templates given, should be the actual \
            number of templates, with each template potentially containing \
            multiple traces.
    :type savedir: Str or bool
    :param savedir: If false, data will be kept in memory, otherwise, data \
        will be stored on disk if memory is tight.
    :type cores: int
    :param cores: Number of cores to use.

    :return: :class: 'numpy.ndarray' objects.  These will contain the \
        correlation sums for each template for this day of data.
    :return: list of ints as number of channels used for each cross-correlation

    .. rubric: Note
    templates must be arranged into a numpy array of numpy arrays.  The inner \
    numpy arrays must be shaped (len(template_trace),) - the outer np.ndarray \
    must be shaped (number of channels in stream*number of templates,), such \
    that if there are 5 traces in the stream (e.g. stream is shaped (5,n)) \
    and there are 10 templates each of length 20, the templates ndarray is \
    shaped (50,20).
    """
    # cimport numpy as np
    import os
    import multiprocessing as mp
    DTYPE = np.float32
    # ctypedef np.float32_t DTYPE_t
    num_cores = cores
    if len(templates) < num_cores:
        num_cores = len(templates)

    # from cython.parallel import parallel, prange
    # from libc.stdlib cimport abort, malloc, free

    # Make some lovely ctypes static declarations
    # cdef int image_ind, template_ind, i
    # cdef np.ndarray ccc = np.zeros(len(stream[0])-len(templates[0])+1)
    # Initialize ndarray for cccsums, which should be one cccsum for
    # each of the ktemplates, with each cccsum as long as the
    # correlation overlap window.
    # cdef np.ndarray cccsums=np.array([np.array([0.0]*len(stream[0])-\
    cccsums = np.array([np.array([0.0]*(len(stream[0])-len(templates[0])+1),
                        dtype=DTYPE)]*ktemplates)
    # Initialize an empty array for the return of the number of channels
    # cdef np.ndarray nchans=np.array([0]*ktemplates, dtype=int)
    nchans = np.array([0] * ktemplates, dtype=int)
    # Loop through the templates, using some kind of clever indexing
    # Note, this is where we could parallelise this!
    image_ind = 0
    template_ind = 0
    j_ind = np.concatenate([np.arange(0, len(stream))] * ktemplates)
    # Array of indexes for stream!
    pool = mp.Pool(processes=num_cores)
    if savedir:
        results = [pool.apply_async(_template_loop, args=(templates[i],
                                                          stream[j_ind[i]],
                                                          delays[i],
                                                          savedir+'/'+str(i),
                                                          i))
                   for i in range(len(templates))]
    else:
        results = [pool.apply_async(_template_loop, args=(templates[i],
                                                          stream[j_ind[i]],
                                                          delays[i], False, i))
                   for i in range(len(templates))]
    pool.close()
    if not savedir:
        ccc_list = [p.get() for p in results]
        ccc_list.sort(key=lambda tup: tup[0])
        ccc_list = [ccc[1] for ccc in ccc_list]
    else:
        # order_list = [p.get() for p in results]
        del order_list
    pool.join()
    print("Finished parallel run")
    for i in range(len(templates)):
        # if i in range(0,len(templates),len(templates)/100):
            # print(str(i/len(templates))+' % read back in')
        # Check if there was data for that station for both the
        if not (np.all(np.isnan(stream[image_ind])) or
                np.all(np.isnan(templates[i]))):
            nchans[template_ind] += 1
        if not savedir:
            cccsums[template_ind] = np.sum([cccsums[template_ind],
                                            ccc_list[i]], axis=0)
        else:
            cccsums[template_ind] = np.sum([cccsums[template_ind],
                                            np.load(savedir+'/'+str(i) +
                                            '.npy')],
                                           axis=0)
            os.remove(savedir+'/'+str(i)+'.npy')
        if image_ind < len(stream) - 1:
            image_ind += 1
        else:
            # Move on to the next template
            image_ind = 0
            template_ind += 1
    # Reshape the array to give what we want
    for i in range(len(cccsums)):
        cccsums[i] = cccsums[i].reshape(len(cccsums[i],))
    return cccsums, nchans


def _template_loop(template, stream, delay, savefile=False, i=0):
    """
    Helper loop for parallelisation
    """
    import cv2
    image = np.append(stream,
                      np.array([0] * int(round(delay))))[int(round(delay)):]
    # Compute the cross correlation
    ccc = cv2.matchTemplate(image.astype(np.float32),
                            template.astype(np.float32),
                            cv2.TM_CCOEFF_NORMED)
    ccc = ccc.T.reshape(len(ccc),)
    if savefile:
        np.save(savefile, ccc)
        del ccc
        return i
    else:
        return(i, ccc)
