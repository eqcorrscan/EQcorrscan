#!/usr/bin/python
r"""Functions to locate newly detected events or re-locate events using hyp2000

The functions herein are borrowed from Tobias Megies' GUI picking application,
Obspyck, the repository for which can be found at
https://github.com/megies/obspyck

The functions have been lightly modified to allow them to run in the absence
of an Obspyck instance.
"""


def doHyp2000(event):
    """
    Writes input files for hyp2000 and starts the hyp2000 program via a
    system call.
    """
    prog_dict = PROGRAMS['hyp_2000']
    files = prog_dict['files']
    precall = prog_dict['PreCall']
    precall(prog_dict)

    f = open(files['phases'], 'wt')
    phases_hypo71 = dicts2hypo71Phases()
    f.write(phases_hypo71)
    f.close()

    f2 = open(files['stations'], 'wt')
    stations_hypo71 = self.dicts2hypo71Stations()
    f2.write(stations_hypo71)
    f2.close()

    self.critical('Phases for Hypo2000:')
    self.catFile(files['phases'], self.critical)
    self.critical('Stations for Hypo2000:')
    self.catFile(files['stations'], self.critical)

    call = prog_dict['Call']
    (msg, err, returncode) = call(prog_dict)
    self.info(msg)
    self.error(err)
    self.critical('--> hyp2000 finished')
    self.catFile(files['summary'], self.critical)

def dicts2hypo71Stations(self):
    """
    Returns the station location information in hypo71
    stations file format as a string. This string can then be written to
    a file.
    """
    fmt = "%6s%02i%05.2f%1s%03i%05.2f%1s%4i\n"
    hypo71_string = ""

    for st in self.streams:
        stats = st[0].stats
        sta = stats.station
        lon = stats.coordinates.longitude
        lon_deg = int(lon)
        lon_min = (lon - lon_deg) * 60.
        lat = stats.coordinates.latitude
        lat_deg = int(abs(lat))
        lat_min = (abs(lat) - abs(lat_deg)) * 60.
        hem_NS = 'N'
        hem_EW = 'E'
        if lat < 0:
            hem_NS = 'S'
        if lon < 0:
            hem_NS = 'W'
        # hypo 71 format uses elevation in meters not kilometers
        ele = stats.coordinates.elevation
        hypo71_string += fmt % (sta, lat_deg, lat_min, hem_NS, lon_deg,
                                lon_min, hem_EW, ele)

    return hypo71_string

def dicts2hypo71Phases(self):
    """
    Returns the pick information in hypo71 phase file format
    as a string. This string can then be written to a file.

    Information on the file formats can be found at:
    http://geopubs.wr.usgs.gov/open-file/of02-171/of02-171.pdf p.30

    Quote:
    The traditional USGS phase data input format (not Y2000 compatible)
    Some fields were added after the original HYPO71 phase format
    definition.

    Col. Len. Format Data
     1    4  A4       4-letter station site code. Also see col 78.
     5    2  A2       P remark such as "IP". If blank, any P time is
                      ignored.
     7    1  A1       P first motion such as U, D, +, -, C, D.
     8    1  I1       Assigned P weight code.
     9    1  A1       Optional 1-letter station component.
    10   10  5I2      Year, month, day, hour and minute.
    20    5  F5.2     Second of P arrival.
    25    1  1X       Presently unused.
    26    6  6X       Reserved remark field. This field is not copied to
                      output files.
    32    5  F5.2     Second of S arrival. The S time will be used if this
                      field is nonblank.
    37    2  A2, 1X   S remark such as "ES".
    40    1  I1       Assigned weight code for S.
    41    1  A1, 3X   Data source code. This is copied to the archive
                      output.
    45    3  F3.0     Peak-to-peak amplitude in mm on Develocorder viewer
                      screen or paper record.
    48    3  F3.2     Optional period in seconds of amplitude read on the
                      seismogram. If blank, use the standard period from
                      station file.
    51    1  I1       Amplitude magnitude weight code. Same codes as P & S.
    52    3  3X       Amplitude magnitude remark (presently unused).
    55    4  I4       Optional event sequence or ID number. This number may
                      be replaced by an ID number on the terminator line.
    59    4  F4.1     Optional calibration factor to use for amplitude
                      magnitudes. If blank, the standard cal factor from
                      the station file is used.
    63    3  A3       Optional event remark. Certain event remarks are
                      translated into 1-letter codes to save in output.
    66    5  F5.2     Clock correction to be added to both P and S times.
    71    1  A1       Station seismogram remark. Unused except as a label
                      on output.
    72    4  F4.0     Coda duration in seconds.
    76    1  I1       Duration magnitude weight code. Same codes as P & S.
    77    1  1X       Reserved.
    78    1  A1       Optional 5th letter of station site code.
    79    3  A3       Station component code.
    82    2  A2       Station network code.
    84-85 2  A2     2-letter station location code (component extension).
    """

    fmtP = "%4s%1sP%1s%1i %15s"
    fmtS = "%12s%1sS%1s%1i\n"
    hypo71_string = ""

    for st in self.streams:
        net = st[0].stats.network
        sta = st[0].stats.station
        pick_p = self.getPick(network=net, station=sta, phase_hint='P')
        pick_s = self.getPick(network=net, station=sta, phase_hint='S')
        if not pick_p and not pick_s:
            continue
        if not pick_p:
            msg = ("Hypo2000 phase file format does not support S pick "
                   "without P pick. Skipping station: %s") % sta
            self.error(msg)
            continue

        # P Pick
        pick = pick_p
        t = pick.time
        hundredth = int(round(t.microsecond / 1e4))
        if hundredth == 100:  # XXX check!!
            t_p = t + 1
            hundredth = 0
        else:
            t_p = t
        date = t_p.strftime("%y%m%d%H%M%S") + ".%02d" % hundredth
        if pick.onset == 'impulsive':
            onset = 'I'
        elif pick.onset == 'emergent':
            onset = 'E'
        else:
            onset = '?'
        if pick.polarity == "positive":
            polarity = "U"
        elif pick.polarity == "negative":
            polarity = "D"
        else:
            polarity = "?"
        try:
            weight = int(pick.extra.weight.value)
        except:
            weight = 0
        hypo71_string += fmtP % (sta, onset, polarity, weight, date)

        # S Pick
        if pick_s:
            if not pick_p:
                err = "Warning: Trying to print a Hypo2000 phase file " + \
                      "with an S phase without P phase.\n" + \
                      "This case might not be covered correctly and " + \
                      "could screw our file up!"
                self.error(err)
            pick = pick_s
            t2 = pick.time
            # if the S time's absolute minute is higher than that of the
            # P pick, we have to add 60 to the S second count for the
            # hypo 2000 output file
            # +60 %60 is necessary if t.min = 57, t2.min = 2 e.g.
            mindiff = (t2.minute - t.minute + 60) % 60
            abs_sec = t2.second + (mindiff * 60)
            if abs_sec > 99:
                err = "Warning: S phase seconds are greater than 99 " + \
                      "which is not covered by the hypo phase file " + \
                      "format! Omitting S phase of station %s!" % sta
                self.error(err)
                hypo71_string += "\n"
                continue
            hundredth = int(round(t2.microsecond / 1e4))
            if hundredth == 100:
                abs_sec += 1
                hundredth = 0
            date2 = "%s.%02d" % (abs_sec, hundredth)
            if pick.onset == 'impulsive':
                onset2 = 'I'
            elif pick.onset == 'emergent':
                onset2 = 'E'
            else:
                onset2 = '?'
            if pick.polarity == "positive":
                polarity2 = "U"
            elif pick.polarity == "negative":
                polarity2 = "D"
            else:
                polarity2 = "?"
            try:
                weight2 = int(pick.extra.weight.value)
            except:
                weight2 = 0
            hypo71_string += fmtS % (date2, onset2, polarity2, weight2)
        else:
            hypo71_string += "\n"

    return hypo71_string
