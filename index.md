---
layout: page
title: Welcome to EQcorrscan
tagline: Matched-filter earthquake detection in Python
---
{% include JB/setup %}

## Introduction

Welcome to EQcorrscan.

These pages make up the general documentation of the EQcorrscan project,
the full API can be found
[here](http://eqcorrscan.readthedocs.org/en/latest/?badge=latest). This
project is an open source, Python project for the detection and analysis
of repeating and near-repeating earthquakes. The codes are hosted on
[github](https://github.com/calum-chamberlain/EQcorrscan), with the
latest release also available on
[pypi](https://pypi.python.org/pypi/EQcorrscan/0.0.8).

This project is in constant development, although the core routines are
relatively stable. Before we migrate to a 0.1.0 release we plan on
implementing more tests and adapting the code to more pep8 friendly
style (currently we mostly follow pep8 guidelines, but we are always
learning!).

You can find more details about the project on the [About](about) page,
and details on how to install and update on the
[Installation](installation) page.

## Posts

<ul class="posts">
  {% for post in site.posts %}
    <li><span>{{ post.date | date_to_string }}</span> &raquo; <a href="{{ BASE_PATH }}{{ post.url }}">{{ post.title }}</a></li>
  {% endfor %}
</ul>
