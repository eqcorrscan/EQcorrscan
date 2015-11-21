---
layout: page
title: About
tagline: Why EQcorrscan?
description: Why EQcorrscan?
permalink: /about/
group: navigation
---
{% include JB/setup %}

EQcorrscan is an earthquake detection package written in Python for OS X and Linux,
and distributed under the LGPL.  The purpose of this package is to detect repeating and
near-repeating earthquakes from continuous seismic data.  

The codes are designed for large-scale problems with multiple (hundreds to
thousands of) templates, run over multiple days of data.  However, it works
equally well for simple problems, with few templates and not much data - you just
don't see the benefit of the parallel processing so much.

At the heart of the routine, multiple templates are correlated with the same
day of data in parallel. using the Python
[multiprocessing](https://docs.python.org/2/library/multiprocessing.html) module.  
As such the parallel processing at this level is CPU bound (you can run as many
templates in parallel as your machine has CPUs).  

In addition to this CPU bound parallel computation, our correlation routines take
advantage of multi-threading by using the matchTemplate routine of the
[openCV](http://opencv.org/) image processing library.

At the top level, multiple days can be run in parallel on multiple nodes, or
machines.  For this, I have taken advantage of the [NeSi](https://www.nesi.org.nz/)
PAN HPC cluster.

This package has been tested on machines ranging from small dual-core laptops
and desktops, to large (multi-thousand core) cluster computers.  The codes
scale well when increasing the processing power, resulting in the ability to
use (tested) **600 templates, through 6.5 years of data in less than 10 hours clock time**.



## Design

These codes have been developed by Calum Chamberlain, during his time as a graduate
student at Victoria University of Wellington, New Zealand.  The codes have been
developed to:

* Work;
* Be readable;
* Have a well-documented [API](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest);
* Be quick and easy to install, and therefore - improve.

To this end, our codes are freely available on github, with both master and
development branches. They also have an initial set of utilities that allow integration
with [Seisan](http://seisan.info/), which was chosen to be the first integrated
standard observatory software, simply due to familiarity (we have plans for SAC,
and QuakeML integration at some point, but anything you want to integrate, you should!).
These utilities also have methods for:

* [clustering](http://eqcorrscan.readthedocs.org/en/latest/submodules/utils.clustering.html);
* [stacking](http://eqcorrscan.readthedocs.org/en/latest/submodules/utils.stacking.html);
* [magnitude calculation](http://eqcorrscan.readthedocs.org/en/latest/submodules/utils.mag_calc.html), including SVD derived magnitudes;
* [peak-finding in noisy data](http://eqcorrscan.readthedocs.org/en/latest/submodules/utils.findpeaks.html).

We do not currently claim that these codes are the best available, but we hope to
ignite collaboration from observational seismologists (and ideally computer
scientists) to help grow this project.  Ideally we would like to see this being
a useful test-base for people to develop and share methods.
