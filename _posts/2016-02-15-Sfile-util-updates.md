---
layout: post
title:  "Sfile_util changes"
tagline: "Important changes"
date:   2016-02-15 17:16:00
categories: updates
---
{% include JB/setup %}

Thanks to the hard work of Chet Hopp we have now integrated into the development
branch the [obspy.core.event](http://docs.obspy.org/master/packages/autogen/obspy.core.event.html) classes.
These allow us to read and write multiple formats, including quakeML, and NonLinLoc.
We have written wrappers to read and write important event and pick information
from and to nordic S-files (seisan format), however these are incomplete, and
do not include focal mechanism solutions, or various other things.

By migrating to obspy.core.event classes we hope this will simplify integration
with other projects, such as seisHub, the growing [ASDF](http://seismic-data.org/),
[hypoDDpy](https://github.com/krischer/hypoDDpy), [StreamPick](https://github.com/calum-chamberlain/StreamPick),
and other Python packages reliant on obspy.  The ability to read and write NonLinLoc
files, as well as seisan files opens up more location software (external to
Python).

On the negative side, this has meant a lot of work for the developers to migrate
PICK and EVENTINFO dependent functions to obspy.core.event usage.  It also means
that scripts that use the PICK or EVENTINFO classes **will no longer work**.  To
simplify transitions we have kept these classes and have written wrappers to
convert between the classes.  We have deliberately not used these wrappers for
the read/write S-file functions as we plan to move away entirely from these
classes.
