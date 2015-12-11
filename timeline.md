---
layout: page
title: Timeline
tagline: What next?
description: What next?
permalink: /timeline/
group: navigation
---
{% include JB/setup %}

EQcorrscan is a young project, and we have yet to integrate a great many things,
however this page serves to outline what we plan to integrate, with a rough
release timeline.  All of these dates are very rough and depend on how much
everyone contributes.  Any comments are welcome
on the issue [tracker](https://github.com/calum-chamberlain/EQcorrscan/issues/3).

##Release 0.0.10 - Jan/Feb 2016
* Incorporate obspy's
[event](https://docs.obspy.org/packages/autogen/obspy.core.event.html) files;
* Integrate
[QuakeML](https://github.com/calum-chamberlain/EQcorrscan/issues/23) event saving;
* Use Seishub for database management;
* Generate a more thorough tutorial using GeoNet data (rather than data downloaded
  with the package) which exemplifies the multi-parallel capabilites.

##Release 0.0.11 - Jun/July 2016
* Development of a subspace detector;
* Incorporate automatic re-picker for detected events;
* Integrate with hypoDDpy and StreamPick for standard observatory
practices including picking of templates, and locations.

##Release 0.1.0 - December 2016
* Development of [tests for all functions](https://github.com/calum-chamberlain/EQcorrscan/issues/11);
* Port to Windows systems;
* Provide Python 3.0 modules.

<!-- {% include JB/comments %} -->
