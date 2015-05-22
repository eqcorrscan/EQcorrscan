def resample(tr, sampling_rate, window='hanning', no_filter=True,
                 strict_length=False):
        """
        Resample trace data using Fourier method. Spectra are linearly
        interpolated if required.

        :type sampling_rate: float
        :param sampling_rate: The sampling rate of the resampled signal.
        :type window: array_like, callable, str, float, or tuple, optional
        :param window: Specifies the window applied to the signal in the
            Fourier domain. Defaults to ``'hanning'`` window. See
            :func:`scipy.signal.resample` for details.
        :type no_filter: bool, optional
        :param no_filter: Deactivates automatic filtering if set to ``True``.
            Defaults to ``True``.
        :type strict_length: bool, optional
        :param strict_length: Leave traces unchanged for which end time of
            trace would change. Defaults to ``False``.

        .. note::

            The :class:`~Trace` object has three different methods to change
            the sampling rate of its data: :meth:`~.resample`,
            :meth:`~.decimate`, and :meth:`~.interpolate`

            Make sure to choose the most appropriate one for the problem at
            hand.

        .. note::

            This operation is performed in place on the actual data arrays. The
            raw data is not accessible anymore afterwards. To keep your
            original data, use :meth:`~obspy.core.trace.Trace.copy` to create
            a copy of your trace object.
            This also makes an entry with information on the applied processing
            in ``stats.processing`` of this trace.

        Uses :func:`scipy.signal.resample`. Because a Fourier method is used,
        the signal is assumed to be periodic.

        .. rubric:: Example

        >>> tr = Trace(data=np.array([0.5, 0, 0.5, 1, 0.5, 0, 0.5, 1]))
        >>> len(tr)
        8
        >>> tr.stats.sampling_rate
        1.0
        >>> tr.resample(4.0)  # doctest: +ELLIPSIS
        <...Trace object at 0x...>
        >>> len(tr)
        32
        >>> tr.stats.sampling_rate
        4.0
        >>> tr.data  # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        array([ 0.5       ,  0.40432914,  0.3232233 ,  0.26903012,  0.25 ...
        """
        from scipy.signal import get_window
        from scipy.fftpack import rfft, irfft
        import numpy as np
        from future.utils import native_str
        factor = tr.stats.sampling_rate / float(sampling_rate)
        # check if end time changes and this is not explicitly allowed
        if strict_length:
            if len(tr.data) % factor != 0.0:
                msg = "End time of trace would change and strict_length=True."
                raise ValueError(msg)
        # do automatic lowpass filtering
        print "Conducting automatic low-pass filtering"
        if not no_filter:
            # be sure filter still behaves good
            if factor > 16:
                msg = "Automatic filter design is unstable for resampling " + \
                      "factors (current sampling rate/new sampling rate) " + \
                      "above 16. Manual resampling is necessary."
                raise ArithmeticError(msg)
            freq = tr.stats.sampling_rate * 0.5 / float(factor)
            tr.filter('lowpassCheby2', freq=freq, maxorder=12)

        orig_dtype = tr.data.dtype
        new_dtype = np.float32 if orig_dtype.itemsize == 4 else np.float64

        print "Resampling in the frequency domain"
        # resample in the frequency domain
        X = rfft(np.require(tr.data, dtype=new_dtype))
        X = np.insert(X, 1, 0)
        if tr.stats.npts % 2 == 0:
            X = np.append(X, [0])
        Xr = X[::2]
        Xi = X[1::2]

        if window is not None:
            print "Windowing"
            if callable(window):
                W = window(np.fft.fftfreq(tr.stats.npts))
            elif isinstance(window, np.ndarray):
                if window.shape != (tr.stats.npts,):
                    msg = "Window has the wrong shape. Window length must " + \
                          "equal the number of points."
                    raise ValueError(msg)
                W = window
            else:
                W = np.fft.ifftshift(get_window(native_str(window),
                                                tr.stats.npts))
            Xr *= W[:tr.stats.npts//2+1]
            Xi *= W[:tr.stats.npts//2+1]

        # interpolate
        print "Interpolating"
        num = int(tr.stats.npts / factor)
        df = 1.0 / (tr.stats.npts * tr.stats.delta)
        dF = 1.0 / num * sampling_rate
        f = df * np.arange(0, tr.stats.npts // 2 + 1, dtype=np.int32)
        nF = num // 2 + 1
        F = dF * np.arange(0, nF, dtype=np.int32)
        Y = np.zeros((2*nF))
        Y[::2] = np.interp(F, f, Xr)
        Y[1::2] = np.interp(F, f, Xi)

        Y = np.delete(Y, 1)
        if num % 2 == 0:
            Y = np.delete(Y, -1)
        tr.data = irfft(Y) * (float(num) / float(tr.stats.npts))
        tr.data = np.require(tr.data, dtype=orig_dtype)
        tr.stats.sampling_rate = sampling_rate

        return tr

