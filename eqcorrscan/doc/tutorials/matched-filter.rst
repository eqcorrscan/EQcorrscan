Matched-filter detection
========================

Simple example
--------------

This example does not work out of the box, you will have to have your own templates
and data, and set things up for this.  However, in principle matched-filtering
can be as simple as:

.. code-block:: python

     from eqcorrscan.core import match_filter
     from eqcorrscan.utils import pre_processing
     from obspy import read

     # Read in and process the daylong data
     st = read('continuous_data')
     # Use the same filtering and sampling parameters as your template!
     st = pre_processing.dayproc(st, lowcut=2, highcut=10, filt_order=4,
                                 samp_rate=50,
                                 starttime=st[0].stats.starttime.date)
     # Read in the templates
     templates = []
     template_names = ['template_1', 'template_2']
     for template_file in template_names:
          templates.append(read(template_file))
     detections = match_filter.match_filter(template_names=template_names,
                                            template_list=templates, st=st,
                                            threshold=8, threshold_type='MAD',
                                            plotvar=False, cores=4)

This will create a list of detections, which are of class detection.  You can
write out the detections to a csv (colon separated) using the detection.write
method, set `append=True` to write all the detections to one file.  Beware though,
if this is set and the file already exists, it will just add on to the old file.

.. code-block:: python

     for detection in detections:
          detection.write('my_first_detections.csv', append=True)


Advanced example
----------------

In this section we will outline using the templates generated in the first tutorial
to scan for similar earthquakes within a day of data.  This small example does not truly exploit the parallel
operations within this package however, so you would be encouraged to think
about where parallel operations occur (*hint, the code can run one template
per CPU*), and why there are --instance and--splits flags in the other
scripts in the github repository (*hint, if you have heaps of memory
and CPUs you can do some brute force day parallelisation!*).

The main processing flow is outlined in the figure below, note the main speedups
in this process are achieved by running multiple templates at once, however this
increases memory usage.  If memory is a problem there are flags (mem_issue) in the
match_filter.py source that can be turned on - the codes will then write temporary
files, which is slower, but can allow for more data crunching at once, your trade-off,
your call.


.. image:: processing_flow.png
     :width: 600px
     :align: center
     :alt: processing_flow.png

.. literalinclude:: ../../tutorials/match_filter.py


SLURM example
-------------

When the authors of EQcorrscan work on large projects, we use grid computers with
the SLURM (Simple Linux Utility for Resource Management) job scheduler installed.
To facilitate ease of setup, what follows is an example of how we run this.

.. code-block:: bash

     #!/bin/bash
     #SBATCH -J MatchTest
     #SBATCH -A ##########
     #SBATCH --time=12:00:00
     #SBATCH --mem=7G
     #SBATCH --nodes=1
     #SBATCH --output=matchout_%a.txt
     #SBATCH --error=matcherr_%a.txt
     #SBATCH --cpus-per-task=16
     #SBATCH --array=0-49

     module load OpenCV/2.4.9-intel-2015a
     module load ObsPy/0.10.3rc1-intel-2015a-Python-2.7.9
     module load joblib/0.8.4-intel-2015a-Python-2.7.9

     srun python2.7 LFEsearch.py --splits 50 --instance $SLURM_ARRAY_TASK_ID


Where we use a script (LFEsearch.py) that accepts splits and instance flags,
this section of the script is as follows:

.. code-block:: python

     Split=False
     instance=False
     if len(sys.argv) == 2:
         flag=str(sys.argv[1])
         if flag == '--debug':
             Test=True
             Prep=False
         elif flag == '--debug-prep':
             Test=False
             Prep=True
         else:
             raise ValueError("I don't recognise the argument, I only know --debug and --debug-prep")
     elif len(sys.argv) == 5:
         # Arguments to allow the code to be run in multiple instances
         Split=True
         Test=False
         Prep=False
         args=sys.argv[1:len(sys.argv)]
         for i in xrange(len(args)):
             if args[i] == '--instance':
                 instance=int(args[i+1])
                 print 'I will run this for instance '+str(instance)
             elif args[i] == '--splits':
                 splits=int(args[i+1])
                 print 'I will divide the days into '+str(splits)+' chunks'

     elif not len(sys.argv) == 1:
         raise ValueError("I only take one argument, no arguments, or two flags with arguments")
     else:
         Test=False
         Prep=False
         Split=False

The full script is not included in EQcorrscan, but is available on request.


