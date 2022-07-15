EQcorrscan FAQs
===============

This is a developing list of frequently asked questions. If you find yourself
experiencing a problem similar to one of these, then try the solution here first!

If your problem/question isn't answered here then check the issues on the github page
and, if there isn't an issue related to your problem/question then open a new one and
we will try to help.

----------------------------------------------------------------------

Usage Questions
---------------

No output to terminal
.....................

EQcorrscan uses `Logging <https://docs.python.org/3/howto/logging.html>`_
to handle output. If you are seeing no output to screen you
probably haven't set up your logging settings. A simple way to do
this is to run code like:

.. code-block:: python
    import logging

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

You can set different levels to get more or less output (DEBUG is the
most output, then INFO, then WARNNG, then ERROR, then CRITICAL). You will need to run this before any calls to EQcorrscan's functions to
get logging output.

---

No correlations computed
........................

Frequently the cause of no correlations being computed is that the
SEED ID (network.station.location.channel) for your template do not
match your continuous data. Check that they match, and try increasing
the logging output (above) to help you find the issue.

---

Everything is done multiple times!
..................................

EQcorrscan uses `multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_
under the hood, which will spawn multiple processes when called. To
ensure that programs do not run multiple times you should always
write your scripts with a form:

.. code-block:: python

    def main():
        # Do the mahi

    if __name__ == "__main__":
        main()

See the `multiprocessing note on "Safe importing of main module" <https://docs.python.org/3/library/multiprocessing.html>`_ for more info.

---

Making templates from SAC files
...............................

While there is support for making templates from SAC, it generally
isn't a good idea, unless you have SAC waveforms of the length that
you want to correlate with. This is for two reasons:
1. Because templates and continuous
data *must* be processed the same, including using the same length
FFT, for correlations to be accurate. If you make your template from
short data you will be forced to correlate with short chunks of data, which is not very efficient
2. The programming team are not SAC users, and do not understand the nuances of where picks can be saved in SAC headers.

Instead: you might find it better to convert your SAC picks into
another obspy readable pick/event format. You can try EQcorrscan's
basic sac_utils to do this, but not everything will be retained.

----------------------------------------------------------------------

Design Questions
----------------

Can I use different lengths for different channels in a template?
.................................................................

Not yet in EQcorrscan - we want this as well, but haven't had time to implement it.
If you want this then we would really appreciate the contribution! There are two
main issues with this that require some thought: 1) How correlations are
computed, and 2) how correlations are weighted in the correlation sum.

Why doesn't EQcorrscan have a GUI?
..................................

This would be cool, and it would be great if someone wants to contribute this,
however, the developers thus far have been focused on other things and don't have
unlimited time.

If you want this, and know how to program GUIs then please do contribute, it would
open EQcorrscan up to more users, which would be great!

Why do you have a functional and object-oriented API?
.....................................................

Simply legacy, when Calum started writing EQcorrscan he didn't know
anything about classes. The functional API is retained so that old
codes can still be run, but if you are starting from scratch please use the OO API where possible.



