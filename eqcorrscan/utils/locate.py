#!/usr/bin/python
r"""Functions to locate newly detected events or re-locate events using hyp2000

The functions herein are borrowed from Tobias Megies and Lion Krischer's GUI
picking application, Obspyck, the repository for which can be found at
https://github.com/megies/obspyck, and which is licensed under GPLv2. See the
following information:
-------------------------------------------------------------------------------
ObsPyck:
Author: Tobias Megies, Lion Krischer
Email: megies@geophysik.uni-muenchen.de
License: GPLv2
Copyright (C) 2010 Tobias Megies, Lion Krischer
-------------------------------------------------------------------------------

These functions have been adapted to allow them to run in the absence of an
Obspyck instance.
"""
#Create PROGRAMS dictionary as per ObsPyck/utils.py (ln. 194)
PROGRAMS = {
    'nlloc': {'filenames': {'exe': "NLLoc", 'phases': "nlloc.obs",
                            'summary': "nlloc.hyp",
                            'scatter': "nlloc.scat"}},
    'hyp_2000': {'filenames': {'exe': "hyp2000", 'control': "bay2000.inp",
                               'phases': "hyp2000.pha",
                               'stations': "stations.dat",
                               'summary': "hypo.prt"}},
    'focmec': {'filenames': {'exe': "rfocmec", 'phases': "focmec.dat",
                             'stdout': "focmec.stdout",
                             'summary': "focmec.out"}}}


def doHyp2000(event, plugindir):
    """
    Writes input files for hyp2000 and starts the hyp2000 program via a
    system call.
    """
    setup_external_programs(plugindir)

    # Removes the station/phase files in directory if they already exist
    prog_dict = PROGRAMS['hyp_2000']
    files = prog_dict['files']
    precall = prog_dict['PreCall']
    precall(prog_dict)
    # Writes phase file for hyp2000
    print 'Phases for Hypo2000:'
    f = open(files['phases'], 'wt')
    phases_hypo71 = dicts2hypo71Phases()
    f.write(phases_hypo71)
    f.close()
    # Writes station file for hyp2000
    print 'Stations for Hypo2000:'
    f2 = open(files['stations'], 'wt')
    stations_hypo71 = dicts2hypo71Stations()
    f2.write(stations_hypo71)
    f2.close()
    # Call hyp2000
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
    print hypo71_string
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
    print hypo71_string
    return hypo71_string


def loadHyp2000Data(self):
    files = PROGRAMS['hyp_2000']['files']
    lines = open(files['summary'], "rt").readlines()
    if lines == []:
        err = "Error: Hypo2000 output file (%s) does not exist!" % \
                files['summary']
        self.error(err)
        return
    # goto origin info line
    while True:
        try:
            line = lines.pop(0)
        except:
            break
        if line.startswith(" YEAR MO DA  --ORIGIN--"):
            break
    try:
        line = lines.pop(0)
    except:
        err = "Error: No location info found in Hypo2000 outputfile " + \
              "(%s)!" % files['summary']
        self.error(err)
        return

    year = int(line[1:5])
    month = int(line[6:8])
    day = int(line[9:11])
    hour = int(line[13:15])
    minute = int(line[15:17])
    seconds = float(line[18:23])
    time = UTCDateTime(year, month, day, hour, minute, seconds)
    lat_deg = int(line[25:27])
    lat_min = float(line[28:33])
    lat = lat_deg + (lat_min / 60.)
    if line[27] == "S":
        lat = -lat
    lon_deg = int(line[35:38])
    lon_min = float(line[39:44])
    lon = lon_deg + (lon_min / 60.)
    if line[38] == " ":
        lon = -lon
    depth = -float(line[46:51]) # depth: negative down!
    rms = float(line[52:57])
    errXY = float(line[58:63])
    errZ = float(line[64:69])

    # goto next origin info line
    while True:
        try:
            line = lines.pop(0)
        except:
            break
        if line.startswith(" NSTA NPHS  DMIN MODEL"):
            break
    line = lines.pop(0)

    #model = line[17:22].strip()
    gap = int(line[23:26])

    line = lines.pop(0)
    model = line[49:].strip()
    # this is to prevent characters that are invalid in QuakeML URIs
    # hopefully handled in the future by obspy/obspy#1018
    model = re.sub(r"[^\w\d\-\.\*\(\)\+\?_~'=,;#/&amp;]", '_', model)

    # assign origin info
    o = Origin()
    self.catalog[0].origins = [o]
    o.clear()
    o.method_id = "/".join([ID_ROOT, "location_method", "hyp2000", "2"])
    o.origin_uncertainty = OriginUncertainty()
    o.quality = OriginQuality()
    ou = o.origin_uncertainty
    oq = o.quality
    o.longitude = lon
    o.latitude = lat
    o.depth = depth * (-1e3)  # meters positive down!
    ou.horizontal_uncertainty = errXY
    ou.preferred_description = "horizontal uncertainty"
    o.depth_errors.uncertainty = errZ * 1e3
    oq.standard_error = rms #XXX stimmt diese Zuordnung!!!?!
    oq.azimuthal_gap = gap
    o.depth_type = "from location"
    o.earth_model_id = "%s/earth_model/%s" % (ID_ROOT, model)
    o.time = time

    # goto station and phases info lines
    while True:
        try:
            line = lines.pop(0)
        except:
            break
        if line.startswith(" STA NET COM L CR DIST AZM"):
            break

    oq.used_phase_count = 0
    oq.extra = AttribDict()
    oq.extra.usedPhaseCountP = {'value': 0, 'namespace': NAMESPACE}
    oq.extra.usedPhaseCountS = {'value': 0, 'namespace': NAMESPACE}
    used_stations = set()
    #XXX caution: we sometimes access the prior element!
    for i in range(len(lines)):
        # check which type of phase
        if lines[i][32] == "P":
            type = "P"
        elif lines[i][32] == "S":
            type = "S"
        else:
            continue
        # get values from line
        station = lines[i][0:6].strip()
        if station == "":
            station = lines[i-1][0:6].strip()
            distance = float(lines[i-1][18:23])
            azimuth = int(lines[i-1][23:26])
            #XXX TODO check, if incident is correct!!
            incident = int(lines[i-1][27:30])
        else:
            distance = float(lines[i][18:23])
            azimuth = int(lines[i][23:26])
            #XXX TODO check, if incident is correct!!
            incident = int(lines[i][27:30])
        used_stations.add(station)
        if lines[i][31] == "I":
            onset = "impulsive"
        elif lines[i][31] == "E":
            onset = "emergent"
        else:
            onset = None
        if lines[i][33] == "U":
            polarity = "positive"
        elif lines[i][33] == "D":
            polarity = "negative"
        else:
            polarity = None
        res = float(lines[i][61:66])
        weight = float(lines[i][68:72])

        # assign synthetic phase info
        pick = self.getPick(station=station, phase_hint=type)
        if pick is None:
            msg = "This should not happen! Location output was read and a corresponding pick is missing!"
            warnings.warn(msg)
        arrival = Arrival(origin=o, pick=pick)
        # residual is defined as P-Psynth by NLLOC!
        # XXX does this also hold for hyp2000???
        arrival.time_residual = res
        arrival.azimuth = azimuth
        arrival.distance = kilometer2degrees(distance)
        arrival.takeoff_angle = incident
        if onset and not pick.onset:
            pick.onset = onset
        if polarity and not pick.polarity:
            pick.polarity = polarity
        # we use weights 0,1,2,3 but hypo2000 outputs floats...
        arrival.time_weight = weight
        o.quality.used_phase_count += 1
        if type == "P":
            o.quality.extra.usedPhaseCountP['value'] += 1
        elif type == "S":
            o.quality.extra.usedPhaseCountS['value'] += 1
        else:
            self.error("Phase '%s' not recognized as P or S. " % type +
                       "Not incrementing P nor S phase count.")
    o.used_station_count = len(used_stations)


def setup_external_programs(pluginpath):
    """
    Sets up temdir, copies program files, fills in PROGRAMS dict, sets up
    system calls for programs.
    Depends on location of plugins, returns temporary directory.

    :param pluginpath: Path to associated plugins (namely hypo2000)
    :type pluginpath: string
    :returns: String representation of temporary directory with program files.
    """
    import os
    import platform
    import tempfile
    import shutil
    import subprocess
    import glob

    if pluginpath is None:
        pluginpath = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isdir(pluginpath):
        msg = "No such directory: '%s'" % pluginpath
        raise IOError(msg)
    tmp_dir = tempfile.mkdtemp(prefix="obspyck-")
    # set binary names to use depending on architecture and platform...
    env = os.environ
    architecture = platform.architecture()[0]
    system = platform.system()
    global SHELL
    if system == "Windows":
        SHELL = True
    else:
        SHELL = False
    # link velocity models ################################################
    if os.path.exists(os.path.join(pluginpath, "VELOCITY_MODELS")):
        os.symlink(os.path.join(pluginpath, "VELOCITY_MODELS"),
                   os.path.join(tmp_dir, "VELOCITY_MODELS"))
    # Setup external programs #############################################
    for prog_basename, prog_dict in PROGRAMS.iteritems():
        prog_srcpath = os.path.join(pluginpath, prog_basename)
        prog_tmpdir = os.path.join(tmp_dir, prog_basename)
        prog_dict['dir'] = prog_tmpdir
        shutil.copytree(prog_srcpath, prog_tmpdir, symlinks=True)
        prog_dict['files'] = {}
        for key, filename in prog_dict['filenames'].iteritems():
            prog_dict['files'][key] = os.path.join(prog_tmpdir, filename)
        prog_dict['files']['exe'] = "__".join(
                [prog_dict['filenames']['exe'], system, architecture])
        # setup clean environment
        prog_dict['env'] = {}
        prog_dict['env']['PATH'] = prog_dict['dir'] + os.pathsep + env['PATH']
        if 'SystemRoot' in env:
            prog_dict['env']['SystemRoot'] = env['SystemRoot']
    # Hyp2000 #############################################################
    prog_dict = PROGRAMS['hyp_2000']
    prog_dict['env']['HYP2000_DATA'] = prog_dict['dir'] + os.sep

    def tmp(prog_dict):
        files = prog_dict['files']
        for file in [files['phases'], files['stations'], files['summary']]:
            if os.path.isfile(file):
                os.remove(file)
        return
    prog_dict['PreCall'] = tmp

    def tmp(prog_dict):
        sub = subprocess.Popen(prog_dict['files']['exe'], shell=SHELL,
                               cwd=prog_dict['dir'], env=prog_dict['env'],
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        input = open(prog_dict['files']['control'], "rt").read()
        (msg, err) = sub.communicate(input)
        if system == "Darwin":
            returncode = sub.returncode
        else:
            returncode = sub.wait()
        return (msg, err, returncode)
    prog_dict['Call'] = tmp
    # NLLoc ###############################################################
    prog_dict = PROGRAMS['nlloc']

    def tmp(prog_dict):
        filepattern = os.path.join(prog_dict['dir'], "nlloc*")
        print filepattern
        for file in glob.glob(filepattern):
            os.remove(file)
        return
    prog_dict['PreCall'] = tmp

    def tmp(prog_dict, controlfilename):
        sub = subprocess.Popen([prog_dict['files']['exe'], controlfilename],
                               cwd=prog_dict['dir'], env=prog_dict['env'],
                               shell=SHELL, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        if system == "Darwin":
            returncode = sub.returncode
        else:
            returncode = sub.wait()
        msg = "".join(sub.stdout.readlines())
        err = "".join(sub.stderr.readlines())
        for pattern, key in [("nlloc.*.*.*.loc.scat", 'scatter'),
                             ("nlloc.*.*.*.loc.hyp", 'summary')]:
            pattern = os.path.join(prog_dict['dir'], pattern)
            newname = os.path.join(prog_dict['dir'], prog_dict['files'][key])
            for file in glob.glob(pattern):
                os.rename(file, newname)
        return (msg, err, returncode)
    prog_dict['Call'] = tmp
    # focmec ##############################################################
    prog_dict = PROGRAMS['focmec']
    def tmp(prog_dict):
        sub = subprocess.Popen(prog_dict['files']['exe'], shell=SHELL,
                cwd=prog_dict['dir'], env=prog_dict['env'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if system == "Darwin":
            returncode = sub.returncode
        else:
            returncode = sub.wait()
        msg = "".join(sub.stdout.readlines())
        err = "".join(sub.stderr.readlines())
        return (msg, err, returncode)
    prog_dict['Call'] = tmp
    #######################################################################
    return tmp_dir
