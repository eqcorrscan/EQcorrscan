.. image:: EQcorrscan_logo.png
    :width: 600px
    :align: left
    :alt: EQcorrscan_logo.png
    :target: https://github.com/calum-chamberlain/EQcorrscan/releases

EQcorrscan
==========

A Python package for the detection and analysis of repeating and near-repeating seismicity.
EQcorrscan contains an efficient, multi-parallel,
:doc:`matched-filter </submodules/core.match_filter>` detection routine (template-matching), as well as
routines to implement :doc:`subspace </submodules/core.subspace>` detection.

Code is stored on |github|, the development branches are |github_dev|, or the
latest stable release can be found |releases_link|.

.. |github| raw:: html

    <a href="https://github.com/eqcorrscan/EQcorrscan" target="_blank">github</a>

.. |releases_link| raw:: html

  <a href="https://github.com/eqcorrscan/EQcorrscan/releases" target="_blank">here</a>

.. |github_dev| raw:: html

  <a href="https://github.com/eqcorrscan/EQcorrscan/tree/develop" target="_blank">here</a>

EQcorrscan uses |Obspy_link| bindings when reading and writing seismic data, and for handling most
of the event metadata, which ensures that detections can be easily migrated between
software.

.. |Obspy_link| raw:: html

  <a href="https://docs.obspy.org/" target="_blank">Obspy</a>

Also within this package are:

* :doc:`Correlation re-picking </submodules/core.lag_calc>`;
* :doc:`Clustering routines for seismic data </submodules/utils.clustering>`;
* :doc:`Peak finding algorithm (basic) </submodules/utils.findpeaks>`;
* :doc:`Stacking routines </submodules/utils.stacking>` including phase-weighted stacking based on Thurber at al. (2014);

This package is written by the EQcorrscan developers, and is 
distributed under the LGPL GNU Licence, Copyright EQcorrscan
developers 2017.

Citation
--------
If you use this package in your work, please cite the following paper:
Chamberlain, C. J., Hopp, C. J., Boese, C. M., Warren-Smith, E., Chambers, D., Chu, S. X., Michailos, K., Townend, J., EQcorrscan: Repeating and near-repeating earthquake detection and analysis in Python. Seismological Research Letters, `2017`

Contents:
---------

.. toctree::
   :numbered:
   :maxdepth: 2

   intro
   installation
   faq
   updates
   tutorial
   api
