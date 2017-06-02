"""
Default parameters for eqcorrscan, only used if the quickstart options are run.

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings
import os
import getpass

from obspy import UTCDateTime

import eqcorrscan


class EQcorrscanParameters:
    """
    Standard class for defining parameters in EQcorrscan.

    These parameters are set using the quickstart script of EQcorrscan for
    matched-filter analysis.  The quickstart functionality is currently in
    development.

    :type template_names: list
    :param template_names: List of str of template names.
    :type lowcut: float
    :param lowcut: Low-cut iun Hz
    :type highcut: float
    :param highcut: High-cut in Hz
    :type filt_order: int
    :param filt_order: Number of corners in filter
    :type samp_rate: float
    :param samp_rate: Desired sampling rate in Hz
    :type debug: int
    :param debug: 0-5 debug level, higher number is more output
    :type startdate: obspy.core.utcdatetime.UTCDateTime
    :param startdate: Start date for match-filter analysis.
    :type enddate: obspy.core.utcdatetime.UTCDateTime
    :param enddate: End date for matched-filter analysis
    :type archive: str
    :param archive: Path to archive, or FDSN name
    :type arc_type: str
    :param arc_type:
        Type of archive (see :func:`eqcorrscan.utils.archive_read` for
        available options.
    :type cores: int
    :param cores: Number of cores to parallel over
    :type plotvar: bool
    :param plotvar: Whether to show plots or not
    :type plotdir: str
    :param plotdir: Output plot location
    :type plot_format: str
    :param plot_format: Output plot type
    :type tempdir: str
    :param tempdir: Location of temporary files (if needed)
    :type threshold: float
    :param threshold: Threshold for matched-filter detection
    :type threshold_type: str
    :param threshold_type:
        Threshold type for matched-filter detection, see
        :func:`eqcorrscan.core.match_filter.match_filter` for available types.
    :type trigger_interval: float
    :param trigger_interval: Minimum trigger interval in seconds.
    """
    def __init__(self, template_names, lowcut, highcut, filt_order, samp_rate,
                 debug, startdate, enddate, archive, arc_type, cores, plotvar,
                 plotdir, plot_format, tempdir, threshold, threshold_type,
                 trigger_interval):
        """
        Standard parameter options.
        """
        if isinstance(template_names, list):
            self.template_names = [str(template_name)
                                   for template_name in template_names]
        else:
            self.template_names = list(str(template_names))
        self.lowcut = float(lowcut)
        self.highcut = float(highcut)
        self.filt_order = int(filt_order)
        self.samp_rate = float(samp_rate)
        if self.samp_rate <= 2 * self.highcut:
            msg = ('Highcut must be less than the Nyquist, setting to ' +
                   str((self.samp_rate / 2.0) - 1))
            warnings.warn(msg)
            self.highcut = (self.samp_rate / 2.0) - 1
        self.debug = int(debug)
        self.startdate = UTCDateTime(startdate)
        self.enddate = UTCDateTime(enddate)
        self.archive = str(archive)
        self.arc_type = str(arc_type)
        # Check that arc_type is in the list of allowed types
        if arc_type.lower() not in ['seishub', 'fdsn', 'day_vols']:
            msg = ' '.join(['arc_type:', arc_type, 'is not of coded type'])
            raise NotImplementedError(msg)
        self.cores = int(cores)
        self.plotvar = bool(plotvar)
        self.plotdir = str(plotdir)
        self.plot_format = str(plot_format)
        self.tempdir = str(tempdir)
        self.threshold = float(threshold)
        self.threshold_type = str(threshold_type)
        # Check that threshold_type is in the allowed types
        if threshold_type not in ['MAD', 'absolute', 'av_chan_corr']:
            msg = ' '.join(['threshold_type:', threshold_type,
                            'is not implemented'])
            raise NotImplementedError(msg)
        self.trigger_interval = float(trigger_interval)

    def __repr__(self):
        return "EQcorrscanParameters()"

    def __str__(self):
        """
        String output for EQcorrscanParameters.

        Note that we are not using the __dict__ attribute because it is \
        not recommened
        """
        print_str = ' '.join(["EQcorrscan parameters:",
                              "\n   template_names:", str(self.template_names),
                              "\n   lowcut:", str(self.lowcut),
                              "\n   highcut:", str(self.highcut),
                              "\n   filt_order:", str(self.filt_order),
                              "\n   samp_rate:", str(self.samp_rate),
                              "\n   debug:", str(self.debug),
                              "\n   startdate:", str(self.startdate),
                              "\n   enddate:", str(self.enddate),
                              "\n   archive:", str(self.archive),
                              "\n   arc_type:", str(self.arc_type),
                              "\n   cores:", str(self.cores),
                              "\n   plotvar:", str(self.plotvar),
                              "\n   plotdir:", str(self.plotdir),
                              "\n   plot_format:", str(self.plot_format),
                              "\n   tempdir:", str(self.tempdir),
                              "\n   threshold:", str(self.threshold),
                              "\n   threshold_type:", str(self.threshold_type),
                              "\n   trigger_interval:",
                              str(self.trigger_interval)
                              ])
        return print_str

    def write(self, outfile='../parameters/EQcorrscan_parameters.txt',
              overwrite=False):
        """
        Function to write the parameters to a file - user readable.

        :type outfile: str
        :param outfile: Full path to filename to store parameters in.
        :type overwrite: bool
        :param overwrite: Whether to overwrite the old file or not.
        """
        outpath = os.sep.join(outfile.split(os.sep)[0:-1])
        if len(outpath) > 0 and not os.path.isdir(outpath):
            msg = ' '.join([os.path.join(outfile.split(os.sep)[0:-1]),
                            'does not exist, check path.'])
            raise IOError(msg)
        # Make sure that the user wants to overwrite the old parameters
        if os.path.isfile(outfile) and not overwrite:
            responding = True
            while responding:
                print(' '.join([outfile, 'exists.  Overwrite? [y/N]']))
                option = raw_input()
                if option.upper() == 'N':
                    raise IOError('File exists, will not overwrite')
                elif option.upper() == 'Y':
                    responding = False
                else:
                    print('Must respond with y or n')
        f = open(outfile, 'w')
        # Write creation info.
        header = ' '.join(['# User:', getpass.getuser(),
                           '\n# Creation date:', str(UTCDateTime()),
                           '\n# EQcorrscan version:',
                           str(eqcorrscan.__version__),
                           '\n\n\n'])
        f.write(header)
        # Write parameter info in a user readable, and parsable format.
        parameters = self.__str__().split('\n')[1:]
        f.write('[eqcorrscan_pars]\n')
        for parameter in parameters:
            f.write(parameter.lstrip() + '\n')
        f.close()
        print('Written parameter file: ' + outfile)


def read_parameters(infile='../parameters/EQcorrscan_parameters.txt'):
    """
    Read the default parameters from file.

    :type infile: str
    :param infile: Full path to parameter file.

    :returns: parameters read from file.
    :rtype: :class:`eqcorrscan.utils.parameters.EQcorrscanParameters`
    """
    try:
        import ConfigParser
    except ImportError:
        import configparser as ConfigParser
    import ast
    f = open(infile, 'r')
    print('Reading parameters with the following header:')
    for line in f:
        if line[0] == '#':
            print(line.rstrip('\n').lstrip('\n'))
    f.close()
    config = ConfigParser.ConfigParser()
    config.read(infile)
    # Slightly tricky list reading
    template_names = list(ast.literal_eval(config.get("eqcorrscan_pars",
                                                      "template_names")))
    parameters = \
        EQcorrscanParameters(template_names=template_names,
                             lowcut=config.get("eqcorrscan_pars", "lowcut"),
                             highcut=config.get("eqcorrscan_pars", "highcut"),
                             filt_order=config.get("eqcorrscan_pars",
                                                   "filt_order"),
                             samp_rate=config.get("eqcorrscan_pars",
                                                  "samp_rate"),
                             debug=config.get("eqcorrscan_pars", "debug"),
                             startdate=config.get("eqcorrscan_pars",
                                                  "startdate"),
                             enddate=config.get("eqcorrscan_pars", "enddate"),
                             archive=config.get("eqcorrscan_pars", "archive"),
                             arc_type=config.get("eqcorrscan_pars",
                                                 "arc_type"),
                             cores=config.get("eqcorrscan_pars", "cores"),
                             plotvar=config.getboolean("eqcorrscan_pars",
                                                       "plotvar"),
                             plotdir=config.get("eqcorrscan_pars", "plotdir"),
                             plot_format=config.get("eqcorrscan_pars",
                                                    "plot_format"),
                             tempdir=ast.literal_eval(config.
                                                      get("eqcorrscan_pars",
                                                          "tempdir")),
                             threshold=config.get("eqcorrscan_pars",
                                                  "threshold"),
                             threshold_type=config.get("eqcorrscan_pars",
                                                       "threshold_type"),
                             trigger_interval=config.get("eqcorrscan_pars",
                                                         "trigger_interval")
                             )

    return parameters


if __name__ == '__main__':
    import doctest
    doctest.testmod()
