EQcorrscan tutorial
===================
THIS TUTOTIAL IS NOT REALLY WRITTEN YET!


You must first set-up your parameter files in the *par* directory.  You can
leave these as the default settings for now, but study the parameters so that
you understand what each one is doing.

Then run the *LFE_search.py* routine in this top directory, running this will
run all the *core* routines to search for templates, generate templates and
compute the cross-channel cross-correlation values for the templates.  Finally
it will output detections from these templates.

The following is verbatim the LFEserach.py routine which outlines usage of this
package:

.. literalinclude:: ../LFEsearch.py
