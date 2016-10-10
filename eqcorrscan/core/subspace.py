r"""This module contains functions relevant to executing subspace detection \
for earthquake catalogs.

We recommend that you read Harris' detailed report on subspace detection \
theory which can be found here: https://e-reports-ext.llnl.gov/pdf/335299.pdf

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import warnings
import time
import h5py
import getpass
import eqcorrscan
import copy
from obspy import Trace, UTCDateTime, Stream
from obspy.core.event import Event, CreationInfo, ResourceIdentifier, Comment,\
    WaveformStreamID, Pick
from eqcorrscan.utils.clustering import svd
from eqcorrscan.utils import findpeaks, pre_processing, stacking, plotting
from eqcorrscan.core.match_filter import DETECTION, extract_from_stream
import matplotlib.pyplot as plt


class Detector(object):
    """
    Class to serve as the base for subspace detections.

    :type name: str
    :param name: Name of subspace detector, used for book-keeping
    :type sampling_rate: float
    :param sampling_rate: Sampling rate in Hz of original waveforms
    :type multiplex: bool
    :param multiplex: Is this detector multiplexed.
    :type stachans: list
    :param stachans:
        List of tuples of (station, channel) used in detector.
        If multiplexed, these must be in the order that multiplexing was done.
    :type lowcut: float
    :param lowcut: Lowcut filter in Hz
    :type highcut: float
    :param highcut: Highcut filter in Hz
    :type filt_order: int
    :param filt_order: Number of corners for filtering
    :type data: numpy.ndarray
    :param data: The actual detector
    :type u: numpy.ndarray
    :param u: Full rank U matrix of left (input) singular vectors.
    :type sigma: numpy.ndarray
    :param sigma: Full rank vector of singular values.
    :type v: numpy.ndarray
    :param v: Full rank right (output) singular vectors.
    :type dimension: int
    :param dimension: Dimension of data.
    """
    def __init__(self, name=None, sampling_rate=None, multiplex=None,
                 stachans=None, lowcut=None, highcut=None,
                 filt_order=None, data=None, u=None, sigma=None, v=None,
                 dimension=None):
        self.name = name
        self.sampling_rate = sampling_rate
        self.multiplex = multiplex
        self.stachans = stachans
        # self.delays = delays
        self.lowcut = lowcut
        self.highcut = highcut
        self.filt_order = filt_order
        self.data = data
        self.u = u
        self.sigma = sigma
        self.v = v
        self.dimension = dimension

    def __repr__(self):
        if self.name:
            out = 'Detector: ' + self.name
        else:
            out = 'Empty Detector object'
        return out

    def __str__(self):
        out = 'Detector object: \n'
        for key in ['name', 'sampling_rate', 'multiplex', 'lowcut', 'highcut',
                    'filt_order', 'dimension']:
            if self.__getattribute__(key):
                out += ('\t' + key + ': ' + str(self.__getattribute__(key)) +
                        '\n')
        return out

    def __eq__(self, other):
        if not isinstance(other, Detector):
            return False
        for key in ['name', 'sampling_rate', 'multiplex', 'lowcut', 'highcut',
                    'filt_order', 'dimension', 'stachans']:
            if not self.__getattribute__(key) == other.__getattribute__(key):
                return False
        for key in ['data', 'u', 'v', 'sigma']:
            list_item = self.__getattribute__(key)
            other_list = other.__getattribute__(key)
            if not len(list_item) == len(other_list):
                return False
            for item, other_item in zip(list_item, other_list):
                if not np.allclose(item, other_item):
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.data)

    def construct(self, streams, lowcut, highcut, filt_order,
                  sampling_rate, multiplex, name, align, shift_len=0,
                  reject=0.3, no_missed=True, plot=False):
        """
        Construct a subspace detector from a list of streams, full rank.

        Subspace detector will be full-rank, further functions can be used \
        to select the desired dimensions.

        :type streams: list
        :param streams:
            List of :class:`obspy.core.stream.Stream` to be used to generate
            the subspace detector.  These should be pre-clustered and aligned.
        :type lowcut: float
        :param lowcut: Lowcut in Hz, can be None to not apply filter
        :type highcut: float
        :param highcut: Highcut in Hz, can be None to not apply filter
        :type filt_order: int
        :param filt_order: Number of corners for filter.
        :type sampling_rate: float
        :param sampling_rate: Desired sampling rate in Hz
        :type multiplex: bool
        :param multiplex:
            Whether to multiplex the data or not.  Data are multiplexed
            according to the method of Harris, see the multi function for
            details.
        :type name: str
        :param name: Name of the detector, used for book-keeping.
        :type align: bool
        :param align:
            Whether to align the data or not - needs to be done at some point
        :type shift_len: float
        :param shift_len: Maximum shift allowed for alignment in seconds.
        :type reject: float
        :param reject:
            Minimum correlation to include traces - only used if align=True.
        :type no_missed: bool
        :param no_missed:
            Reject streams with missed traces, defaults to True. A missing
            trace from lots of events will reduce the quality of the subspace
            detector if multiplexed.  Only used when multi is set to True.
        :type plot: bool
        :param plot: Whether to plot the alignment stage or not.

        .. note::
            The detector will be normalized such that the data, before
            computing the singular-value decomposition, will have unit energy.
            e.g. We divide the amplitudes of the data by the L1 norm of the
            data.

        .. warning::
            EQcorrscan's alignment will attempt to align over the whole data
            window given.  For long (more than 2s) chunks of data this can give
            poor results and you might be better off using the
            :func:`eqcorrscan.utils.stacking.align_traces` function externally,
            focusing on a smaller window of data.  To do this you would align
            the data prior to running construct.
        """
        self.lowcut = lowcut
        self.highcut = highcut
        self.filt_order = filt_order
        self.sampling_rate = sampling_rate
        self.name = name
        self.multiplex = multiplex
        # Pre-process data
        p_streams, stachans = \
            _subspace_process(streams=copy.deepcopy(streams),
                              lowcut=lowcut, highcut=highcut,
                              filt_order=filt_order,
                              sampling_rate=sampling_rate, multiplex=multiplex,
                              align=align, shift_len=shift_len, reject=reject,
                              plot=plot, no_missed=no_missed)
        # Compute the SVD, use the cluster.SVD function
        v, sigma, u, svd_stachans = svd(stream_list=p_streams, full=True)
        if not multi:
            stachans = [tuple(stachan.split('.')) for stachan in svd_stachans]
        self.stachans = stachans
        # self.delays = delays
        self.u = u
        self.v = v
        self.sigma = sigma
        self.data = copy.deepcopy(u)  # Set the data matrix to be full rank U.
        self.dimension = np.inf
        return self

    def partition(self, dimension):
        """
        Partition subspace into desired dimension.

        :type dimension: int
        :param dimension: Maximum dimension to use.
        """
        # Take leftmost 'dimension' input basis vectors
        for i, channel in enumerate(self.u):
            if channel.shape[1] < dimension:
                raise IndexError('Channel is max dimension %s'
                                 % channel.shape[1])
            self.data[i] = channel[:, 0:dimension]
        self.dimension = dimension
        return self

    def energy_capture(self):
        """
        Calculate the average percentage energy capture for this subspace.

        :return: Percentage energy capture
        :rtype: float
        """
        percent_captrue = 0
        if np.isinf(self.dimension):
            return 100
        for channel in self.sigma:
            fc = np.sum(channel[0:self.dimension]) / np.sum(channel)
            percent_captrue += fc
        return 100 * (percent_captrue / len(self.sigma))

    def detect(self, st, threshold, trig_int, moveout=0, min_trig=0,
               process=True, extract_detections=False, debug=0):
        """
        Detect within continuous data using the subspace method.

        :type st: obspy.core.stream.Stream
        :param st:
            Un-processed stream to detect within using the subspace detector.
        :type threshold: float
        :param threshold: Threshold value for detections between 0-1
        :type trig_int: float
        :param trig_int: Minimum trigger interval in seconds.
        :type moveout: float
        :param moveout:
            Maximum allowable moveout window for non-multiplexed, network
            detection.  See note.
        :type min_trig: int
        :param min_trig:
            Minimum number of stations exceeding threshold for non-multiplexed,
            network detection. See note.
        :type process: bool
        :param process:
            Whether or not to process the stream according to the parameters
            defined by the detector.  Default is True, which will process the
            data.
        :type extract_detections: bool
        :param extract_detections:
            Whether to extract waveforms for each detection or not, if True
            will return detections and streams.
        :type debug: int
        :param debug: Debug output level from 0-5.

        :return: list of :class:`eqcorrscan.core.match_filter.DETECTION`
        :rtype: list

        .. warning::
            Subspace is currently in beta, see note in the subspace tutorial
            for information.

        .. note::
            If running in bulk with detectors that all have the same
            parameters then you can pre-process the data and set process to
            False.  This will speed up this detect function dramatically.

        .. warning::
            If the detector and stream are multiplexed then they must
            contain the same channels and multiplexed in the same order. This
            is handled internally when process=True, but if running in bulk
            you must take care.

        .. note::
            Non-multiplexed, network detection.  When the detector is
            not multiplexed, but there are multiple channels within the
            detector, we do not stack the single-channel detection statistics
            because we do not have a one-size-fits-all solution for computing
            delays for a subspace detector (if you want to implement one, then
            please contribute it!).  Therefore, these parameters provide a
            means for declaring a network coincidence trigger using
            single-channel detection statistics, in a similar fashion to the
            commonly used network-coincidence trigger with energy detection
            statistics.
        """
        return _detect(detector=self, st=st, threshold=threshold,
                       trig_int=trig_int, moveout=moveout, min_trig=min_trig,
                       process=process, extract_detections=extract_detections,
                       debug=debug)

    def write(self, filename):
        """
        Write detector to a file - uses HDF5 file format.

        Meta-data are stored alongside numpy data arrays. See h5py.org for \
        details of the methods.

        :type filename: str
        :param filename: Filename to save the detector to.
        """
        f = h5py.File(filename, "w")
        # Must store eqcorrscan version number, username would be useful too.
        data_group = f.create_group(name="data")
        for i, data in enumerate(self.data):
            dset = data_group.create_dataset(name="data_" + str(i),
                                             shape=data.shape,
                                             dtype=data.dtype)
            dset[...] = data
        data_group.attrs['length'] = len(self.data)
        data_group.attrs['name'] = self.name.encode("ascii", "ignore")
        data_group.attrs['sampling_rate'] = self.sampling_rate
        data_group.attrs['multiplex'] = self.multiplex
        data_group.attrs['lowcut'] = self.lowcut
        data_group.attrs['highcut'] = self.highcut
        data_group.attrs['filt_order'] = self.filt_order
        data_group.attrs['dimension'] = self.dimension
        data_group.attrs['user'] = getpass.getuser()
        data_group.attrs['eqcorrscan_version'] = str(eqcorrscan.__version__)
        # Convert station-channel list to something writable
        ascii_stachans = ['.'.join(stachan).encode("ascii", "ignore")
                          for stachan in self.stachans]
        stachans = f.create_dataset(name="stachans",
                                    shape=(len(ascii_stachans),),
                                    dtype='S10')
        stachans[...] = ascii_stachans
        u_group = f.create_group("u")
        for i, u in enumerate(self.u):
            uset = u_group.create_dataset(name="u_" + str(i),
                                          shape=u.shape, dtype=u.dtype)
            uset[...] = u
        u_group.attrs['length'] = len(self.u)
        sigma_group = f.create_group("sigma")
        for i, sigma in enumerate(self.sigma):
            sigmaset = sigma_group.create_dataset(name="sigma_" + str(i),
                                                  shape=sigma.shape,
                                                  dtype=sigma.dtype)
            sigmaset[...] = sigma
        sigma_group.attrs['length'] = len(self.sigma)
        v_group = f.create_group("v")
        for i, v in enumerate(self.v):
            vset = v_group.create_dataset(name="v_" + str(i),
                                          shape=v.shape, dtype=v.dtype)
            vset[...] = v
        v_group.attrs['length'] = len(self.v)
        f.flush()
        f.close()
        return self

    def read(self, filename):
        """
        Read detector from a file, must be HDF5 format.

        Reads a Detector object from an HDF5 file, usually created by \
        eqcorrscan.

        :type filename: str
        :param filename: Filename to save the detector to.
        """
        f = h5py.File(filename, "r")
        self.data = []
        for i in range(f['data'].attrs['length']):
            self.data.append(f['data']['data_' + str(i)].value)
        self.u = []
        for i in range(f['u'].attrs['length']):
            self.u.append(f['u']['u_' + str(i)].value)
        self.sigma = []
        for i in range(f['sigma'].attrs['length']):
            self.sigma.append(f['sigma']['sigma_' + str(i)].value)
        self.v = []
        for i in range(f['v'].attrs['length']):
            self.v.append(f['v']['v_' + str(i)].value)
        self.stachans = [tuple(stachan.decode('ascii').split('.'))
                         for stachan in f['stachans'].value]
        self.dimension = f['data'].attrs['dimension']
        self.filt_order = f['data'].attrs['filt_order']
        self.highcut = f['data'].attrs['highcut']
        self.lowcut = f['data'].attrs['lowcut']
        self.multiplex = bool(f['data'].attrs['multiplex'])
        self.sampling_rate = f['data'].attrs['sampling_rate']
        if isinstance(f['data'].attrs['name'], str):
            self.name = f['data'].attrs['name']
        else:
            self.name = f['data'].attrs['name'].decode('ascii')
        return self

    def plot(self, stachans='all', size=(10, 7), show=True):
        """
        Plot the output basis vectors for the detector at the given dimension.

        Corresponds to the first n horizontal vectors of the V matrix.

        :type stachans: list
        :param stachans: list of tuples of station, channel pairs to plot.
        :type stachans: list
        :param stachans: List of tuples of (station, channel) to use.  Can set\
            to 'all' to use all the station-channel pairs available. If \
            detector is multiplexed, will just plot that.
        :type size: tuple
        :param size: Figure size.
        :type show: bool
        :param show: Whether or not to show the figure.

        :returns: Figure
        :rtype: matplotlib.pyplot.Figure
        """
        if stachans == 'all' and not self.multiplex:
            stachans = self.stachans
        elif self.multiplex:
            stachans = [('multi', ' ')]
        if np.isinf(self.dimension):
            nrows = self.data[0].shape[1]
        fig, axes = plt.subplots(nrows=nrows, ncols=len(stachans),
                                 sharex=True, sharey=True, figsize=size)
        x = np.arange(len(self.v[0]), dtype=np.float32)
        if self.multiplex:
            x /= len(self.stachans) * self.sampling_rate
        else:
            x /= self.sampling_rate
        for column, stachan in enumerate(stachans):
            channel = self.v[column]
            for row, vector in enumerate(channel.T[0:nrows]):
                if len(stachans) == 1:
                    if nrows == 1:
                        axis = axes
                    else:
                        axis = axes[row]
                else:
                    axis = axes[row, column]
                if row == 0:
                    axis.set_title('.'.join(stachan))
                axis.plot(x, vector, 'k', linewidth=1.1)
                if column == 0:
                    axis.set_ylabel('Basis %s' % (row + 1))
                if row == nrows - 1:
                    axis.set_xlabel('Time (s)')
        plt.subplots_adjust(hspace=0.05)
        plt.subplots_adjust(wspace=0.05)
        if show:
            plt.show()
        return fig


def _detect(detector, st, threshold, trig_int, moveout=0, min_trig=0,
            process=True, extract_detections=False, debug=0):
    """
    Detect within continuous data using the subspace method.

    Not to be called directly, use the detector.detect method.

    :type detector: eqcorrscan.core.subspace.Detector
    :param detector: Detector to use.
    :type st: obspy.core.stream.Stream
    :param st: Un-processed stream to detect within using the subspace \
        detector
    :type threshold: float
    :param threshold: Threshold value for detections between 0-1
    :type trig_int: float
    :param trig_int: Minimum trigger interval in seconds.
    :type moveout: float
    :param moveout: Maximum allowable moveout window for non-multiplexed,
        network detection.  See note.
    :type min_trig: int
    :param min_trig: Minimum number of stations exceeding threshold for \
        non-multiplexed, network detection. See note.
    :type process: bool
    :param process: Whether or not to process the stream according to the \
        parameters defined by the detector.  Default is to process the \
        data (True).
    :type extract_detections: bool
    :param extract_detections: Whether to extract waveforms for each \
        detection or not, if true will return detections and streams.
    :type debug: int
    :param debug: Debug output level from 0-5.

    :return: list of detections
    :rtype: list of eqcorrscan.core.match_filter.DETECTION
    """
    from eqcorrscan.core import subspace_statistic
    detections = []
    # First process the stream
    if process:
        if debug > 0:
            print('Processing Stream')
        stream, stachans = _subspace_process(streams=[st.copy()],
                                             lowcut=detector.lowcut,
                                             highcut=detector.highcut,
                                             filt_order=detector.filt_order,
                                             sampling_rate=detector.
                                             sampling_rate,
                                             multiplex=detector.multiplex,
                                             stachans=detector.stachans,
                                             parallel=True,
                                             align=False,
                                             shift_len=None,
                                             reject=False)
    else:
        # Check the sampling rate at the very least
        for tr in st:
            if not tr.stats.sampling_rate == detector.sampling_rate:
                raise ValueError('Sampling rates do not match.')
        stream = [st]
        stachans = detector.stachans
    outtic = time.clock()
    if debug > 0:
        print('Computing detection statistics')
    stats = np.zeros((len(stream[0]),
                      len(stream[0][0]) - len(detector.data[0][0]) + 1),
                     dtype=np.float32)
    for det_channel, in_channel, i in zip(detector.data, stream[0],
                                          np.arange(len(stream[0]))):
        stats[i] = subspace_statistic.\
            det_statistic(detector=det_channel.astype(np.float32),
                          data=in_channel.data.astype(np.float32))
        if debug > 0:
            print(stats[i].shape)
        if debug > 3:
            plt.plot(stats[i])
            plt.show()
        # Hard typing in Cython loop requires float32 type.
    # statistics
    if detector.multiplex:
        trig_int_samples = (len(detector.stachans) *
                            detector.sampling_rate * trig_int)
    else:
        trig_int_samples = detector.sampling_rate * trig_int
    if debug > 0:
        print('Finding peaks')
    peaks = []
    for i in range(len(stream[0])):
        peaks.append(findpeaks.find_peaks2_short(arr=stats[i],
                                                 thresh=threshold,
                                                 trig_int=trig_int_samples,
                                                 debug=debug))
    if not detector.multiplex:
        # Conduct network coincidence triggering
        peaks = findpeaks.coin_trig(peaks=peaks,
                                    samp_rate=detector.sampling_rate,
                                    moveout=moveout, min_trig=min_trig,
                                    stachans=stachans, trig_int=trig_int)
    else:
        peaks = peaks[0]
    if len(peaks) > 0:
        for peak in peaks:
            if detector.multiplex:
                detecttime = st[0].stats.starttime + (peak[1] /
                                                      (detector.sampling_rate *
                                                       len(detector.stachans)))
            else:
                detecttime = st[0].stats.starttime + (peak[1] /
                                                      detector.sampling_rate)
            rid = ResourceIdentifier(id=detector.name + '_' +
                                     str(detecttime),
                                     prefix='smi:local')
            ev = Event(resource_id=rid)
            cr_i = CreationInfo(author='EQcorrscan',
                                creation_time=UTCDateTime())
            ev.creation_info = cr_i
            # All detection info in Comments for lack of a better idea
            thresh_str = 'threshold=' + str(threshold)
            ccc_str = 'detect_val=' + str(peak[0])
            used_chans = 'channels used: ' +\
                ' '.join([str(pair) for pair in detector.stachans])
            ev.comments.append(Comment(text=thresh_str))
            ev.comments.append(Comment(text=ccc_str))
            ev.comments.append(Comment(text=used_chans))
            for stachan in detector.stachans:
                tr = st.select(station=stachan[0], channel=stachan[1])
                if tr:
                    net_code = tr[0].stats.network
                else:
                    net_code = ''
                pick_tm = detecttime
                wv_id = WaveformStreamID(network_code=net_code,
                                         station_code=stachan[0],
                                         channel_code=stachan[1])
                ev.picks.append(Pick(time=pick_tm, waveform_id=wv_id))
            detections.append(DETECTION(detector.name,
                                        detecttime,
                                        len(detector.stachans),
                                        peak[0],
                                        threshold,
                                        'subspace', detector.stachans,
                                        event=ev))
    outtoc = time.clock()
    print('Detection took %s seconds' % str(outtoc - outtic))
    if extract_detections:
        detection_streams = extract_from_stream(st, detections)
        return detections, detection_streams
    return detections


def _subspace_process(streams, lowcut, highcut, filt_order, sampling_rate,
                      multiplex, align, shift_len, reject, no_missed=True,
                      stachans=None, parallel=False, plot=False):
    """
    Process stream data, internal function.

    :type streams: list
    :param streams: List of obspy.core.stream.Stream to be used to \
        generate the subspace detector.  These should be pre-clustered \
        and aligned.
    :type lowcut: float
    :param lowcut: Lowcut in Hz, can be None to not apply filter
    :type highcut: float
    :param highcut: Highcut in Hz, can be None to not apply filter
    :type filt_order: int
    :param filt_order: Number of corners for filter.
    :type sampling_rate: float
    :param sampling_rate: Desired sampling rate in Hz
    :type multiplex: bool
    :param multiplex: Whether to multiplex the data or not.  Data are \
        multiplexed according to the method of Harris, see the multi \
        function for details.
    :type stachans: list of tuple
    :param stachans: list of tuples of (station, channel) to use.
    :type align: bool
    :param align: Whether to align the data or not - needs to be done \
        at some point
    :type shift_len: float
    :param shift_len: Maximum shift allowed for alignment in seconds.
    :type reject: float
    :param reject: Minimum correlation for traces, only used if align=True.
    :type no_missed: bool
    :param: no_missed: Reject streams with missed traces, defaults to True. \
        A missing trace from lots of events will reduce the quality of the \
        subspace detector if multiplexed.  Only used when multi is set to True.
    :type plot: bool
    :param plot: Passed down to align traces - used to check alignment process.

    :return: Processed streams
    :rtype: list
    :return: Station, channel pairs in order
    :rtype: list of tuple
    :return: List of delays
    :rtype: list
    """
    from multiprocessing import Pool, cpu_count
    processed_streams = []
    if not stachans:
        input_stachans = list(set([(tr.stats.station, tr.stats.channel)
                                   for st in streams for tr in st.sort()]))
    else:
        input_stachans = stachans
    input_stachans.sort()  # Make sure stations and channels are in order
    # Check that all channels are the same length in seconds
    first_length = len(streams[0][0].data) /\
        streams[0][0].stats.sampling_rate
    for st in streams:
        for tr in st:
            if not len(tr) / tr.stats.sampling_rate == first_length:
                msg = 'All channels of all streams must be the same length'
                raise IOError(msg)
    for st in streams:
        if not parallel:
            processed_stream = Stream()
            for stachan in input_stachans:
                dummy, tr = _internal_process(st=st, lowcut=lowcut,
                                              highcut=highcut,
                                              filt_order=filt_order,
                                              sampling_rate=sampling_rate,
                                              first_length=first_length,
                                              stachan=stachan, debug=0)
                processed_stream += tr
            processed_streams.append(processed_stream)
        else:
            pool = Pool(processes=cpu_count())
            results = [pool.apply_async(_internal_process, (st,),
                                        {'lowcut': lowcut,
                                         'highcut': highcut,
                                         'filt_order': filt_order,
                                         'sampling_rate': sampling_rate,
                                         'first_length': first_length,
                                         'stachan': stachan,
                                         'debug': 0,
                                         'i': i})
                       for i, stachan in enumerate(input_stachans)]
            pool.close()
            processed_stream = [p.get() for p in results]
            pool.join()
            processed_stream.sort(key=lambda tup: tup[0])
            processed_stream = Stream([p[1] for p in processed_stream])
            processed_streams.append(processed_stream)
        if no_missed and multiplex:
            for tr in processed_stream:
                if np.count_nonzero(tr.data) == 0:
                    processed_streams.remove(processed_stream)
                    print('Removed stream with empty trace')
                    break
    if align:
        processed_streams = align_design(design_set=processed_streams,
                                         shift_len=shift_len,
                                         reject=reject, multiplex=multiplex,
                                         plot=plot, no_missed=no_missed)
    output_streams = []
    for processed_stream in processed_streams:
        if len(processed_stream) == 0:
            # If we have removed all of the traces then onwards!
            continue
        # Need to order the stream according to input_stachans
        _st = Stream()
        for stachan in input_stachans:
            tr = processed_stream.select(station=stachan[0],
                                         channel=stachan[1])
            if len(tr) >= 1:
                _st += tr[0]
            elif multiplex and len(tr) == 0:
                raise IndexError('Missing data for %s.%s' %
                                 (stachan[0], stachan[1]))
        if multiplex:
            st = multi(stream=_st)
            st = Stream(Trace(st))
            st[0].stats.station = 'Multi'
            st[0].stats.sampling_rate = sampling_rate
        else:
            st = _st
        for tr in st:
            # Normalize the data
            norm = np.linalg.norm(tr.data)
            if not norm == 0:
                tr.data /= norm
        output_streams.append(st)
    return output_streams, input_stachans


def _internal_process(st, lowcut, highcut, filt_order, sampling_rate,
                      first_length, stachan, debug, i=0):
    tr = st.select(station=stachan[0], channel=stachan[1])
    if len(tr) == 0:
        tr = Trace(np.zeros(int(first_length * sampling_rate)))
        tr.stats.station = stachan[0]
        tr.stats.channel = stachan[1]
        tr.stats.sampling_rate = sampling_rate
        tr.stats.starttime = st[0].stats.starttime  # Do this to make more
        # sensible plots
        warnings.warn('Padding stream with zero trace for ' +
                      'station ' + stachan[0] + '.' + stachan[1])
    elif len(tr) == 1:
        tr = tr[0]
        tr.detrend('simple')
        tr = pre_processing.process(tr=tr, lowcut=lowcut, highcut=highcut,
                                    filt_order=filt_order,
                                    samp_rate=sampling_rate, debug=debug,
                                    seisan_chan_names=False)
    else:
        msg = ('Multiple channels for ' + stachan[0] + '.' +
               stachan[1] + ' in a single design stream.')
        raise IOError(msg)
    return i, tr


def read_detector(filename):
    """
    Read detector from a filename.

    :type filename: str
    :param filename: Filename to save the detector to.

    :return: Detector object
    :rtype: eqcorrscan.core.subspace.Detector
    """
    detector = Detector()
    detector.read(filename=filename)
    return detector


def multi(stream):
    """
    Internal multiplexer for multiplex_detect.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream to multiplex

    :return: trace of multiplexed data
    :rtype: obspy.core.trace.Trace

    .. Note: Requires all channels to be the same length.

    Maps a standard multiplexed stream of seismic data to a single traces of \
    multiplexed data as follows:

    Input:
    x = [x1, x2, x3, ...]
    y = [y1, y2, y3, ...]
    z = [z1, z2, z3, ...]

    Output:
    xyz = [x1, y1, z1, x2, y2, z2, x3, y3, z3, ...]
    """
    stack = stream[0].data
    for tr in stream[1:]:
        stack = np.dstack(np.array([stack, tr.data]))
    multiplex = stack.reshape(stack.size, )
    return multiplex


def align_design(design_set, shift_len, reject, multiplex, no_missed=True,
                 plot=False):
    """
    Align individual traces within streams of the design set.

    Perform before Detector.construct to align traces before computing the \
    singular value decomposition.

    :type design_set: list
    :param design_set: List of obspy.core.stream.Stream to be aligned
    :type shift_len: float
    :param shift_len: Maximum shift (plus/minus) in seconds.
    :type reject: float
    :param reject: Minimum correlation for traces, only used if align=True.
    :type multiplex: bool
    :param multiplex: If you are going to multiplex the data, then there has \
        to be data for all channels, so we will pad with zeros, otherwise \
        there is no need.
    :type no_missed: bool
    :param: no_missed: Reject streams with missed traces, defaults to True. \
        A missing trace from lots of events will reduce the quality of the \
        subspace detector if multiplexed.  Only used when multi is set to True.
    :type plot: bool
    :param plot: Whether to plot the aligned traces as we go or not.

    :rtype: list
    :return: List of obspy.core.stream.Stream of aligned streams

    .. Note:: Assumes only one trace for each channel for each stream in the \
        design_set. If more are present will only use the first one.

    .. Note:: Will cut all traces to be the same length as required for the \
        svd, this length will be the shortest trace length - 2 * shift_len
    """
    trace_lengths = [tr.stats.endtime - tr.stats.starttime for st in design_set
                     for tr in st]
    clip_len = min(trace_lengths) - (2 * shift_len)
    stachans = list(set([(tr.stats.station, tr.stats.channel)
                         for st in design_set for tr in st]))
    remove_set = []
    for stachan in stachans:
        trace_list = []
        trace_ids = []
        for i, st in enumerate(design_set):
            tr = st.select(station=stachan[0], channel=stachan[1])
            if len(tr) > 0:
                trace_list.append(tr[0])
                trace_ids.append(i)
            if len(tr) > 1:
                warnings.warn('Too many matches for %s %s' % (stachan[0],
                                                              stachan[1]))
        shift_len_samples = int(shift_len * trace_list[0].stats.sampling_rate)
        shifts, cccs = stacking.align_traces(trace_list=trace_list,
                                             shift_len=shift_len_samples,
                                             positive=True)
        for i, shift in enumerate(shifts):
            st = design_set[trace_ids[i]]
            start_t = st.select(station=stachan[0],
                                channel=stachan[1])[0].stats.starttime
            start_t += shift_len
            start_t -= shift
            st.select(station=stachan[0],
                      channel=stachan[1])[0].trim(start_t,
                                                  start_t + clip_len)
            if cccs[i] < reject:
                if multiplex and not no_missed:
                    st.select(station=stachan[0],
                              channel=stachan[1])[0].data =\
                        np.zeros(int(clip_len *
                                     (st.select(station=stachan[0],
                                                channel=stachan[1])[0].
                                      stats.sampling_rate) + 1))
                    warnings.warn('Padding stream with zero trace for ' +
                                  'station ' + stachan[0] + '.' + stachan[1])
                elif multiplex and no_missed:
                    remove_set.append(st)
                    warnings.warn('Will remove stream due to low-correlation')
                    continue
                else:
                    st.remove(st.select(station=stachan[0],
                              channel=stachan[1])[0])
                    print('Removed channel with correlation at %s' % cccs[i])
                    continue
    if no_missed:
        for st in remove_set:
            if st in design_set:
                design_set.remove(st)
    if plot:
        for stachan in stachans:
            trace_list = []
            for st in design_set:
                tr = st.select(station=stachan[0], channel=stachan[1])
                if len(tr) > 0:
                    trace_list.append(tr[0])
            if len(trace_list) > 1:
                plotting.multi_trace_plot(traces=trace_list, corr=True,
                                          stack=None, title='.'.join(stachan))
            else:
                print('No plot for you, only one trace left after rejection')
    return design_set


def subspace_detect(detectors, stream, threshold, trig_int, moveout=0,
                    min_trig=1, parallel=True, num_cores=None):
    """
    Conduct subspace detection with chosen detectors.

    :type detectors: list
    :param detectors:
        list of :class:`eqcorrscan.core.subspace.Detector` to be used
        for detection.
    :type stream: obspy.core.stream.Stream
    :param stream: Stream to detect within.
    :type threshold: float
    :param threshold:
        Threshold between 0 and 1 for detection, see :func:`Detector.detect`
    :type trig_int: float
    :param trig_int: Minimum trigger interval in seconds.
    :type moveout: float
    :param moveout:
        Maximum allowable moveout window for non-multiplexed, network
        detection.  See note.
    :type min_trig: int
    :param min_trig:
        Minimum number of stations exceeding threshold for non-multiplexed,
        network detection. See note in :func:`Detector.detect`.
    :type parallel: bool
    :param parallel: Whether to run detectors in parallel in groups.
    :type num_cores: int
    :param num_cores:
        How many cpu cores to use if parallel==True. If set to None (default),
        will use all available cores.

    :rtype: list
    :return:
        List of :class:`eqcorrscan.core.match_filter.DETECTION` detections.

    .. Note::
        This will loop through your detectors using their detect method.
        If the detectors are multiplexed it will run groups of detectors with
        the same channels at the same time.
    """
    from multiprocessing import Pool, cpu_count
    # First check that detector parameters are the same
    parameters = []
    detections = []
    for detector in detectors:
        parameter = (detector.lowcut, detector.highcut,
                     detector.filt_order, detector.sampling_rate,
                     detector.multiplex, detector.stachans)
        if parameter not in parameters:
            parameters.append(parameter)
    for parameter_set in parameters:
        parameter_detectors = []
        for detector in detectors:
            det_par = (detector.lowcut, detector.highcut, detector.filt_order,
                       detector.sampling_rate, detector.multiplex,
                       detector.stachans)
            if det_par == parameter_set:
                parameter_detectors.append(detector)
        stream, stachans = \
            _subspace_process(streams=[stream.copy()],
                              lowcut=parameter_set[0],
                              highcut=parameter_set[1],
                              filt_order=parameter_set[2],
                              sampling_rate=parameter_set[3],
                              multiplex=parameter_set[4],
                              stachans=parameter_set[5],
                              parallel=True, align=False, shift_len=None,
                              reject=False)
        if not parallel:
            for detector in parameter_detectors:
                detections += _detect(detector=detector, st=stream[0],
                                      threshold=threshold, trig_int=trig_int,
                                      moveout=moveout, min_trig=min_trig,
                                      process=False, extract_detections=False,
                                      debug=0)
        else:
            if num_cores:
                ncores = num_cores
            else:
                ncores = cpu_count()
            pool = Pool(processes=ncores)
            results = [pool.apply_async(_detect,
                                        args=(detector, stream[0], threshold,
                                              trig_int, moveout, min_trig,
                                              False, False, 0))
                       for detector in parameter_detectors]
            pool.close()
            _detections = [p.get() for p in results]
            pool.join()
            for d in _detections:
                if isinstance(d, list):
                    detections += d
                else:
                    detections.append(d)
    return detections
