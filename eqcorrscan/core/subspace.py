r"""This module contains functions relevant to executing subspace detection \
for earthquake catalogs.

We recommend that you read Harris' detailed report on subspace detection \
theory which can be found here: https://e-reports-ext.llnl.gov/pdf/335299.pdf

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
import warnings
import time
import h5py
import getpass
import eqcorrscan
from obspy import Trace, UTCDateTime, Stream
from obspy.core.event import Event, CreationInfo, ResourceIdentifier, Comment,\
    WaveformStreamID, Pick
from eqcorrscan.utils.clustering import svd
from eqcorrscan.utils import findpeaks, pre_processing
from eqcorrscan.core.match_filter import DETECTION, extract_from_stream



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
    :param stachans: List of tuples of (station, channel) used in detector. \
        If multiplexed, these must be in the order that multiplexing was done.
    :type lowcut: float
    :param lowcut: Lowcut filter in Hz
    :type highcut: float
    :param highcut: Highcut filter in Hz
    :type filt_order: int
    :param filt_order: Number of corners for filtering
    :type data: np.ndarray
    :param data: The actual detector
    :type u: np.ndarray
    :param u: Full rank U matrix of left (input) singular vectors.
    :type sigma: np.ndarray
    :param sigma: Full rank vector of singular values.
    :type v: np.ndarray
    :param v: Full rank right (output) singular vectors.
    :type dimension: int
    :param dimension: Dimension of data.
    """
    def __init__(self, name=None, sampling_rate=None, multiplex=None,
                 stachans=None, lowcut=None, highcut=None, filt_order=None,
                 data=None, u=None, sigma=None, v=None, dimension=None):
        self.name = name
        self.sampling_rate = sampling_rate
        self.multiplex = multiplex
        self.stachans = stachans
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
            if not np.allclose(self.__getattribute__(key),
                               other.__getattribute__(key)):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.data)

    def construct(self, streams, lowcut, highcut, filt_order,
                  sampling_rate, multiplex, name):
        """
        Construct a subspace detector from a list of streams, full rank.

        Subspace detector will be full-rank, further functions can be used \
        to select the desired dimensions.

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
        :type name: str
        :param name: Name of the detector, used for book-keeping.

        .. note:: The detector will be normalized such that the data, before \
            computing the singular-value decomposition, will have unit energy. \
            e.g. We divide the amplitudes of the data by the L1 norm of the \
            data.
        """
        self.lowcut = lowcut
        self.highcut = highcut
        self.filt_order = filt_order
        self.sampling_rate = sampling_rate
        self.name = name
        self.multiplex = multiplex
        # Pre-process data
        p_streams, stachans = _subspace_process(streams=streams,
                                                lowcut=lowcut,
                                                highcut=highcut,
                                                filt_order=filt_order,
                                                sampling_rate=sampling_rate,
                                                multiplex=multiplex)
        # Compute the SVD, use the cluster.SVD function
        v, sigma, u, dummy = svd(stream_list=p_streams)
        self.stachans = stachans
        self.u = u
        self.v = v
        self.sigma = sigma
        self.data = u  # Set the data matrix to be full rank U.
        self.dimension = np.inf
        return self

    def partition(self, dimension):
        """
        Partition subspace into desired dimension.

        :type dimension: int
        :param dimension: Maximum dimension to use.
        """
        # Take leftmost 'dimension' input basis vectors
        for i, channel in enumerate(self.data):
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

    def detect(self, st, threshold, trig_int, process=True,
               extract_detections=False, debug=0):
        """
        Detect within continuous data using the subspace method.

        :type st: obspy.core.stream.Stream
        :param st: Un-processed stream to detect within using the subspace \
            detector
        :type threshold: float
        :param threshold: Threshold value for detections between 0-1
        :type trig_int: float
        :param trig_int: Minimum trigger interval in seconds.
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

        .. note:: If running in bulk with detectors that all have the same \
            parameters then you can pre-process the data and set process to \
            False.  This will speed up this detect function dramatically.

        .. warning:: If the detector and stream are multiplexed then they must \
            contain the same channels and multiplexed in the same order. This \
            is handled internally when process=True, but if running in bulk \
            you must take care.
        """
        from eqcorrscan.core import subspace_statistic
        detections = []
        # First process the stream
        if process:
            if debug > 0:
                print('Processing Stream')
            stream, stachans = _subspace_process(streams=[st.copy()],
                                                 lowcut=self.lowcut,
                                                 highcut=self.highcut,
                                                 filt_order=self.filt_order,
                                                 sampling_rate=self.
                                                 sampling_rate,
                                                 multiplex=self.multiplex,
                                                 stachans=self.stachans,
                                                 parallel=True)
        else:
            # Check the sampling rate at the very least
            stream = [st]
            for tr in stream[0]:
                if not tr.stats.sampling_rate == self.sampling_rate:
                    raise ValueError('Sampling rates do not match.')
        outtic = time.clock()
        if debug > 0:
            print('Computing detection statistics')
        stats = np.zeros(len(stream[0][0]) - len(self.data[0][0]) + 1,
                         dtype=np.float32)
        for det_channel, in_channel in zip(self.data, stream[0]):
            stats += subspace_statistic.\
                det_statistic(detector=det_channel.astype(np.float32),
                              data=in_channel.data.astype(np.float32))
            # Hard typing in Cython loop requires float32 type.
        stats /= len(stream[0])  # Average the non-multiplexed detection
        # statistics
        if self.multiplex:
            trig_int_samples = (len(self.stachans) *
                                self.sampling_rate * trig_int)
        else:
            trig_int_samples = self.sampling_rate * trig_int
        if debug > 0:
            print('Finding peaks')
        peaks = findpeaks.find_peaks2_short(arr=stats, thresh=threshold,
                                            trig_int=trig_int_samples,
                                            debug=debug)
        if peaks:
            for peak in peaks:
                if self.multiplex:
                    detecttime = st[0].stats.starttime + (peak[1] /
                                                          (self.sampling_rate *
                                                           len(self.stachans)))
                else:
                    detecttime = st[0].stats.starttime + (peak[1] /
                                                          self.sampling_rate)
                rid = ResourceIdentifier(id=self.name + '_' +
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
                    ' '.join([str(pair) for pair in self.stachans])
                ev.comments.append(Comment(text=thresh_str))
                ev.comments.append(Comment(text=ccc_str))
                ev.comments.append(Comment(text=used_chans))
                for stachan in self.stachans:
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
                detections.append(DETECTION(self.name,
                                            detecttime,
                                            len(self.stachans),
                                            peak[0],
                                            threshold,
                                            'subspace', self.stachans,
                                            event=ev))
        outtoc = time.clock()
        print('Detection took %s seconds' % str(outtoc - outtic))
        if extract_detections:
            detection_streams = extract_from_stream(st, detections)
            return detections, detection_streams
        return detections

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
        # Reshape data into a numpy array
        data_write = np.array(self.data)
        u_array = np.array(self.u)
        sigma_array = np.array(self.sigma)
        v_array = np.array(self.v)
        dset = f.create_dataset(name="data", shape=data_write.shape,
                                dtype=data_write.dtype)
        dset[...] = data_write
        dset.attrs['name'] = self.name
        dset.attrs['sampling_rate'] = self.sampling_rate
        dset.attrs['multiplex'] = self.multiplex
        dset.attrs['lowcut'] = self.lowcut
        dset.attrs['highcut'] = self.highcut
        dset.attrs['filt_order'] = self.filt_order
        dset.attrs['dimension'] = self.dimension
        dset.attrs['user'] = getpass.getuser()
        dset.attrs['eqcorrscan_version'] = str(eqcorrscan.__version__)
        # Convert station-channel list to something writable
        ascii_stachans = ['.'.join(stachan).encode("ascii", "ignore")
                          for stachan in self.stachans]
        stachans = f.create_dataset(name="stachans",
                                    shape=(len(ascii_stachans),),
                                    dtype='S10')
        stachans[...] = ascii_stachans
        uset = f.create_dataset(name="u", shape=u_array.shape,
                                dtype=u_array.dtype)
        uset[...] = u_array
        sigmaset = f.create_dataset(name="sigma", shape=sigma_array.shape,
                                    dtype=sigma_array.dtype)
        sigmaset[...] = sigma_array
        vset = f.create_dataset(name="v", shape=v_array.shape,
                                dtype=v_array.dtype)
        vset[...] = v_array
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
        self.data = list(f['data'].value)
        self.u = list(f['u'].value)
        self.v = list(f['v'].value)
        self.sigma = list(f['sigma'].value)
        self.stachans = [tuple(stachan.decode('ascii').split('.'))
                         for stachan in f['stachans'].value]
        self.dimension = f['data'].attrs['dimension']
        self.filt_order = f['data'].attrs['filt_order']
        self.highcut = f['data'].attrs['highcut']
        self.lowcut = f['data'].attrs['lowcut']
        self.multiplex = bool(f['data'].attrs['multiplex'])
        self.sampling_rate = f['data'].attrs['sampling_rate']
        self.name = f['data'].attrs['name'].decode('ascii')
        return self


def _subspace_process(streams, lowcut, highcut, filt_order, sampling_rate,
                     multiplex, stachans=None, parallel=False):
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

    :return: Processed streams
    :rtype: list
    :return: Station, channel pairs in order
    :rtype: list of tuple
    """
    from multiprocessing import Pool, cpu_count
    processed_streams = []
    if not stachans:
        input_stachans = list(set([(tr.stats.station, tr.stats.channel)
                                   for st in streams for tr in st]))
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
        if multiplex:
            st = multi(stream=processed_stream)
            st = Stream(Trace(st))
            st[0].stats.station = 'Multi'
            st[0].stats.sampling_rate = sampling_rate
        for tr in st:
            # Normalize the data
            tr.data /= np.linalg.norm(tr.data)
        processed_streams.append(st)
    return processed_streams, input_stachans


def _internal_process(st, lowcut, highcut, filt_order, sampling_rate,
                      first_length, stachan, debug, i=0):
    tr = st.select(station=stachan[0], channel=stachan[1])
    if len(tr) == 0:
        tr = Trace(np.zeros(first_length * sampling_rate))
        tr.stats.station = stachan[0]
        tr.stats.channel = stachan[1]
        warnings.warn('Padding stream with zero trace for ' +
                      'station ' + stachan[0] + '.' + stachan[1])
    elif len(tr) == 1:
        tr = tr[0]
        pre_processing.process(tr=tr, lowcut=lowcut, highcut=highcut,
                               filt_order=filt_order, samp_rate=sampling_rate,
                               debug=debug)
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
    :type normalize: bool
    :param normalize: Normalize data before multiplexing, will normalize to \
        rms amplitude.

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


def subspace_detect(detectors, stream, threshold, trig_int):
    """
    Conduct subspace detection with chosen detectors.

    :type detectors: list
    :param detectors: list of eqcorrscan.core.subspace.Detector to be used for \
        detection
    :param stream:
    :param threshold:
    :param trig_int:
    :return:
    """
    # First check that detector parameters are the same
    parameters = []
    for detector in detectors:
        parameters.append(detector.lowcut, detector.highcut,
                          detector.filt_order, detector.sampling_rate,
                          detector.multiplex)
    parameters = list(set(parameters))
    if not len(parameters) == 1:
        msg = ('Multiple parameters used for detectors, group your ' +
               'detectors first')
        raise IOError(msg)
    # Now process the stream according to those parameters

