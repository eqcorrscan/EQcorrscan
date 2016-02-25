---
layout: page
title: Installation
permalink: /installation/
group: navigation
---

Since October 2015 EQcorrscan has been listed on [pypi](https://pypi.python.org/pypi/EQcorrscan/0.0.8); as such installation
is as simple as:

{% highlight bash %}
pip install EQcorrscan
{% endhighlight %}

Updating should be done using the following command:

{% highlight bash %}
pip install -U --no-deps EQcorrscan
{% endhighlight %}

This will ignore dependency updates, which have not been needed so far, if this
changes you will be notified along with the release update.

If you haven't installed openCV yourself, which EQcorrscan depends on, but
cannot install itself, it will warn you when it builds.  You should install
openCV.  Currently EQcorrscan uses openCV-2, openCV-3 looks very promising, but
Python wrappers are not yet complete.  To install openCV-2 with python bindings
you should use the instructions on the [openCV docs](http://docs.opencv.org/2.4/doc/tutorials/introduction/table_of_content_introduction/table_of_content_introduction.html#table-of-content-introduction).
