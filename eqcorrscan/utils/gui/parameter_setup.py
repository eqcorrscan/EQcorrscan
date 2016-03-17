"""
Functions to make a simple gui for setting default parameter values.  These \
will be written to a parameter file of the users choosing and will use the \
utils.parameters functions for writing and reading.

This file is part of EQcorrscan, copyright the authors 2016.

Licence: LGPL
"""

class ParameterSetup:
    def __init__(self, master, par=False):
        """
        GUI for selecting default parameters - will write parameters to file \
        of users choosing.

        :type master: Tk
        :param master: Tkinter window
        :type par: EQcorrscanParameters
        :param par: Default parameters to start-up with.
        """
        from tkinter import Label, Button, Entry, DoubleVar, StringVar, IntVar
        from tkinter import BooleanVar, OptionMenu
        from eqcorrscan.utils import parameters
        from obspy import UTCDateTime

        # Set the default par, only if they don't already exist.
        if not par:
            par = parameters.EQcorrscanParameters([''], 2, 10, 4, 20, 2,
                                                  '1900-01-01', '2300-01-01',
                                                  '', 'seishub', 4, False, '',
                                                  '', False, 8, 'MAD', 6)
        # Callback functions for all variables (ugly)

        def update_template_names(*args):
            par.template_names = [name.strip() for name in
                                  template_names.get().split(',')]
            template_names.set(', '.join(par.template_names))

        def update_lowcut(*args):
            par.lowcut = lowcut.get()
            lowcut.set(par.lowcut)

        def update_highcut(*args):
            par.highcut = highcut.get()
            highcut.set(par.highcut)

        def update_filt_order(*args):
            par.filt_order = filt_order.get()
            filt_order.set(par.filt_order)

        def update_samp_rate(*args):
            par.samp_rate = samp_rate.get()
            samp_rate.set(par.samp_rate)

        def update_debug(*args):
            par.debug = debug.get()
            debug.set(par.debug)

        def update_startdate(*args):
            par.startdate = UTCDateTime(startdate.get())
            startdate.set(str(par.startdate))

        def update_enddate(*args):
            par.enddate = UTCDateTime(enddate.get())
            enddate.set(str(par.enddate))

        def update_archive(*args):
            par.archive = archive.get()
            archive.set(par.archive)

        def update_arc_type(*args):
            par.arc_type = arc_type.get()
            arc_type.set(par.arc_type)

        def update_cores(*args):
            par.cores = cores.get()
            cores.set(par.cores)

        def update_plotvar(*args):
            par.plotvar = plotvar.get()
            plotvar.set(par.plotvar)

        def update_plot_format(*args):
            par.plot_format = plot_format.get()
            plot_format.set(par.plot_format)

        def update_tempdir(*args):
            par.tempdir = tempdir.get()
            tempdir.set(par.tempdir)

        def update_threshold(*args):
            par.threshold = threshold.get()
            threshold.set(par.threshold)

        def update_threshold_type(*args):
            par.threshold_type = threshold_type.get()
            threshold_type.set(par.threshold_type)

        def update_plotdir(*args):
            par.plotdir = plotdir.get()
            plotdir.set(par.plotdir)

        def update_trigger_interval(*args):
            par.trigger_interval = trigger_interval.get()
            trigger_interval.set(par.trigger_interval)
        # Set some grid parameters
        nrows = 25
        ncolumns = 3
        self.master = master
        master.title("EQcorrscan parameter setup")
        self.label = Label(master, text="EQcorrscan parameter setup")
        self.label.grid(column=0, columnspan=ncolumns, row=0)

        # Set up parameter input
        self.t_names_label = Label(master, text="Template names", anchor='e')
        self.t_names_label.grid(column=0, row=1)
        template_names = StringVar()
        template_names.set(', '.join(par.template_names))
        self.t_names_box = Entry(master, bd=2, textvariable=template_names)
        self.t_names_box.grid(column=1, row=1)
        template_names.trace("w", update_template_names)
        self.t_names_lookup = Button(master, text="Lookup",
                                     command=lambda: self.get_template_names(par))
        self.t_names_lookup.grid(column=2, row=1)

        self.lowcut_label = Label(master, text="Lowcut (Hz)", anchor='e')
        self.lowcut_label.grid(column=0, row=2)
        lowcut = DoubleVar()
        lowcut.set(par.lowcut)
        self.lowcut_box = Entry(master, bd=2, textvariable=lowcut)
        self.lowcut_box.grid(column=1, row=2)
        lowcut.trace("w", update_lowcut)

        self.highcut_label = Label(master, text="Highcut (Hz)", anchor='e')
        self.highcut_label.grid(column=0, row=3)
        highcut = DoubleVar()
        highcut.set(par.highcut)
        self.highcut_box = Entry(master, bd=2, textvariable=highcut)
        self.highcut_box.grid(column=1, row=3)
        highcut.trace("w", update_highcut)

        self.filt_order_label = Label(master, text="Filter order")
        self.filt_order_label.grid(column=0, row=4)
        filt_order = DoubleVar()
        filt_order.set(par.filt_order)
        self.filt_order_box = Entry(master, bd=2, textvariable=filt_order)
        self.filt_order_box.grid(column=1, row=4)
        filt_order.trace("w", update_filt_order)

        self.samp_rate_label = Label(master, text="Sample rate (Hz)")
        self.samp_rate_label.grid(column=0, row=5)
        samp_rate = DoubleVar()
        samp_rate.set(par.samp_rate)
        self.samp_rate_box = Entry(master, bd=2, textvariable=samp_rate)
        self.samp_rate_box.grid(column=1, row=5)
        samp_rate.trace("w", update_samp_rate)

        self.debug_label = Label(master, text="Debug")
        self.debug_label.grid(column=0, row=6)
        debug = IntVar()
        debug.set(par.debug)
        self.debug_box = Entry(master, bd=2, textvariable=debug)
        self.debug_box.grid(column=1, row=6)
        debug.trace("w", update_debug)

        self.startdate_label = Label(master, text="Start date (yyyy-mm-dd)")
        self.startdate_label.grid(column=0, row=6)
        startdate = StringVar()
        startdate.set(par.startdate)
        self.startdate_box = Entry(master, bd=2, textvariable=startdate)
        self.startdate_box.grid(column=1, row=6)
        startdate.trace("w", update_startdate)

        self.enddate_label = Label(master, text="End date (yyyy-mm-dd)")
        self.enddate_label.grid(column=0, row=8)
        enddate = StringVar()
        enddate.set(par.enddate)
        self.enddate_box = Entry(master, bd=2, textvariable=enddate)
        self.enddate_box.grid(column=1, row=8)
        enddate.trace("w", update_enddate)

        self.archive_label = Label(master, text="Archive")
        self.archive_label.grid(column=0, row=9)
        archive = StringVar()
        archive.set(par.archive)
        self.archive_box = Entry(master, bd=2, textvariable=archive)
        self.archive_box.grid(column=1, row=9)
        archive.trace("w", update_archive)
        self.archive_lookup = Button(master, text="Lookup",
                                     command=lambda: self.get_archive(par))
        self.archive_lookup.grid(column=2, row=9)

        self.arc_type_label = Label(master, text="Archive type")
        self.arc_type_label.grid(column=0, row=10)
        arc_type = StringVar()
        arc_type.set(par.arc_type)
        self.arc_type_box = OptionMenu(master, arc_type,
                                       "seishub", "fdsn", "day_vols")
        self.arc_type_box.grid(column=1, row=10)
        arc_type.trace("w", update_arc_type)

        self.cores_label = Label(master, text="Number of cores")
        self.cores_label.grid(column=0, row=11)
        cores = IntVar()
        cores.set(par.cores)
        self.cores_box = Entry(master, bd=2, textvariable=cores)
        self.cores_box.grid(column=1, row=11)
        cores.trace("w", update_cores)

        self.plotvar_label = Label(master, text="Plotting on/off")
        self.plotvar_label.grid(column=0, row=12)
        plotvar = BooleanVar()
        plotvar.set(par.plotvar)
        self.plotvar_box = Entry(master, bd=2, textvariable=plotvar)
        self.plotvar_box.grid(column=1, row=12)
        plotvar.trace("w", update_plotvar)

        self.plotdir_label = Label(master, text="Plot directory")
        self.plotdir_label.grid(column=0, row=13)
        plotdir = StringVar()
        plotdir.set(par.plotdir)
        self.plotdir_box = Entry(master, bd=2, textvariable=plotdir)
        self.plotdir_box.grid(column=1, row=13)
        plotdir.trace("w", update_plotdir)

        self.plot_format_label = Label(master, text="Plot format")
        self.plot_format_label.grid(column=0, row=14)
        plot_format = StringVar()
        plot_format.set(par.plot_format)
        self.plot_format_box = Entry(master, bd=2, textvariable=plot_format)
        self.plot_format_box.grid(column=1, row=14)
        plot_format.trace("w", update_plot_format)

        self.tempdir_label = Label(master, text="Temporary directory")
        self.tempdir_label.grid(column=0, row=15)
        tempdir = StringVar()
        tempdir.set(par.tempdir)
        self.tempdir_box = Entry(master, bd=2, textvariable=tempdir)
        self.tempdir_box.grid(column=1, row=15)
        tempdir.trace("w", update_tempdir)

        self.threshold_label = Label(master, text="Threshold")
        self.threshold_label.grid(column=0, row=16)
        threshold = DoubleVar()
        threshold.set(par.threshold)
        self.threshold_box = Entry(master, bd=2, textvariable=threshold)
        self.threshold_box.grid(column=1, row=16)
        threshold.trace("w", update_threshold)

        self.threshold_type_label = Label(master, text="Threshold type")
        self.threshold_type_label.grid(column=0, row=17)
        threshold_type = StringVar()
        threshold_type.set(par.threshold_type)
        self.threshold_type_box = Entry(master, bd=2,
                                        textvariable=threshold_type)
        self.threshold_type_box.grid(column=1, row=17)
        threshold_type.trace("w", update_threshold_type)

        self.trigger_interval_label = Label(master,
                                            text="Minimum trigger " +
                                            "interval (s)")
        self.trigger_interval_label.grid(column=0, row=18)
        trigger_interval = DoubleVar()
        trigger_interval.set(par.trigger_interval)
        self.trigger_interval_box = Entry(master, bd=2,
                                          textvariable=trigger_interval)
        self.trigger_interval_box.grid(column=1, row=18)
        trigger_interval.trace("w", update_trigger_interval)

        # End of user editable section, now we have read/write buttons
        self.read_button = Button(master, text="Read parameters",
                                  command=self.read_par)
        self.read_button.grid(column=0, row=nrows-2)

        self.write_button = Button(master, text="Write parameters",
                                   command=lambda: self.write_par(par))
        self.write_button.grid(column=1, row=nrows-2)

        self.close_button = Button(master, text="Quit",
                                   command=master.destroy)
        self.close_button.grid(column=0, columnspan=ncolumns, row=nrows-1)

    def read_par(self):
        """
        Function to open a file-browser and to select a parameter file.
        """
        from eqcorrscan.utils import parameters
        from tkFileDialog import askopenfilename
        parameter_filename = askopenfilename()
        try:
            par = parameters.read_parameters(parameter_filename)
            # Start a new instance
            self.master.destroy()
            run(par=par)
        except IOError:
            print 'No such file'
            return
        except TypeError:
            print 'Invalid parameter file'
            return

    def write_par(self, par):
        from eqcorrscan.utils import parameters
        from tkFileDialog import asksaveasfilename
        parameter_filename = asksaveasfilename()
        print(parameter_filename)
        # Set overwrite to true because asksavefilename already checks this.
        par.write(parameter_filename, overwrite=True)

    def get_template_names(self, par):
        from tkFileDialog import askopenfilenames
        par.template_names = askopenfilenames()
        self.master.destroy()
        run(par=par)

    def get_archive(self, par):
        from tkFileDialog import askdirectory
        par.archive = askdirectory()
        par.arc_type = 'day_vols'
        self.master.destroy()
        run(par=par)


def run(par=False):
    from tkinter import Tk
    root = Tk()
    parameter_gui = ParameterSetup(root, par)
    root.mainloop()
