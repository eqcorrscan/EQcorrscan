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

Logging

---

No correlations computed
........................

Check SEED IDs

---

Everything is done multiple times!
..................................

Multiprocessing issue - encapsulate in function

---

Making templates from SAC files
...............................

It probably isn't a good idea... Try converting the SAC file to an event and go from there.

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

Legacy - use the OO API where possible



