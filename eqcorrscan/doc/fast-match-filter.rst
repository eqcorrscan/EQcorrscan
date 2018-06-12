Fast Matched Filter
===================

For those of you looking to combine the best of time-domain correlations, GPUs and
EQcorrscan, the developers of Fast Matched Filter have written some very fast, very
parallel correlation routines.  These can be plugged into EQcorrscan to allow you
to use the EQcorrscan front-end and the Fast Matched Filter backend.

For information on the Fast Matched Filter project see their |FMFgithub| page.  If
you do use their codes, please cite the paper they recommend (at the time of writing,
this is |FMFpaper|).

Fast Matched Filter is most useful on massively multi-threaded systems, such as those with more
than 40 CPU cores, or NVIDIA GPU cards (other vendors GPU cards are not supported by Fast
Matched Filter as of 12/06/2018). It is wicked fast, and very light on memory compared to
EQcorrscan's frequency domain calculations. For smaller, CPU-only machines, the frequency domain
wins out though.

To install Fast Matched Filter:


To register it as a correlation routine and use it:



.. |FMFgithub| raw:: html

    <a href="https://github.com/beridel/fast_matched_filter" target="_blank">github</a>

.. |FMFpaper| raw:: html

    <a href="https://doi.org/10.1785/0220170181" target="_blank">Beaucé, Eric, W. B. Frank, and Alexey Romanenko (2017). Fast matched-filter (FMF):
     an efficient seismic matched-filter search for both CPU and GPU architectures. Seismological
     Research Letters, doi: 10.1785/0220170181</a>