---
layout: page
title: Contributing
tagline: Help us out!
description: How to contribute
permalink: /contributing/
group: navigation
---
{% include JB/setup %}


## First contributions

If you want to contribute to the project there are a number of ways to do it:

* There will be bugs in the code, some systems will break, some will give you
an odd answer, **be wary**.  If and when you do find bugs, if they are not big
then you should create an [issue](https://github.com/calum-chamberlain/EQcorrscan/issues),
or check the current issues and add your comments to that.  When flagging a bug
please give useful things like terminal output, error messages, what input you
used and if possible, highlight where in the source code you think the bug is.
* If you just want to get involved and give something back, we are in the
 process of developing tests for all our functions, currently the major ones
 are implemented, but there are lots left to go.
* If you want to learn about the codes, how they work, check where they could be
improved etc., then fork the [repo](https://github.com/calum-chamberlain/EQcorrscan)
and scan through them - make them more readable, improve the docs, anything you
like, it's all appreciated.
* If there is something you want to see in the package, you could check to see if
we have it on our [list of things to develop](https://github.com/calum-chamberlain/EQcorrscan/issues/3),
and if not, add it!

If you want to have a go at adding a feature (which would be great)
then please fork the [repo](https://github.com/calum-chamberlain/EQcorrscan) on
github, and create your own feature branch.  We are trying to follow the
[github flow](http://nvie.com/posts/a-successful-git-branching-model/) branching
model, so you might want to install their command-line
[tools](https://github.com/nvie/gitflow/wiki/Installation).

Then when you are done, create a pull request and we can integrate it back so
everyone can use it :)

## Style-guide

We are trying to migrate to a truly
[pep8-friendly](https://www.python.org/dev/peps/pep-0008/https://www.python.org/dev/peps/pep-0008/)
coding style, to that end I would recommend installing
[flake8](https://flake8.readthedocs.org/en/latest/).  Our main two developers
use the [atom](https://atom.io/https://atom.io/)
text editor along with [linter-flake8](https://atom.io/packages/linter-flake8),
although I often use vim as well.  Up to you.

We use [obspy](http://docs.obspy.org/http://docs.obspy.org/)
heavily, so try to follow their
[style guide](http://docs.obspy.org/coding_style.html), including
the standard import conventions for numpy and matplotlib:

{% highlight python %}
#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

{% endhighlight %}



As usual:

* Don't start and end function names with double underscores, these
are reserved for python;
* Never use dashes (I hate them);
* Use descriptive names, in contrast to the pep257 complaints about overly
complicated variable names, if it helps the readability then do it;
* If you want an internal or hidden function, start the function name with a
single underscore;
* Don't use CamelCase for module names, separate words with underscores.

## Documentation

As you will have seen from the [API](http://eqcorrscan.readthedocs.org/en/latest/),
docs are good.  We generate the docs using [sphinx](http://sphinx-doc.org/)
which can then be compiled either locally, or by
[readthedocs](https://readthedocs.org/).  In this way, docstrings are including
in the source code.

* We work on the principle that **If it's not documented, it doesn't exist**.
* We also *try to* work on the principle of **just enough docs**: don't over document the
functions within the script, if you think they need more explaining then the
chances are you are doing it wrong.  We work on this basis because the more
complicated the docs are, the less likely we are to maintain them.  
* If you find
yourself explaining a method in the docstrings: **don't**, put a link to a paper.
