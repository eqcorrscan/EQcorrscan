from obspy import UTCDateTime, Stream, read, Catalog
from copy import copy
from datetime import time
from obspy.core.event import Event, ResourceIdentifier
from obspy.core.event.origin import Pick, Origin
from obspy.core.event.magnitude import Magnitude
from obspy.core.event.base import WaveformStreamID
from obspy.core.util.attribdict import AttribDict
def catalog_times(filename):
	cat_times = []
	cat_mags = []
	cat_depths = []
#filestr=filename[12:16]	
	event_ind = []
	with open(filename) as f:
		content = f.readlines()
		for i, line in enumerate(content):
			if 'J2012' in line:
				event_ind.append(i)
	for i, ind in enumerate(event_ind):
	#	new_ev = Event()
		line = content[ind]
		if i < len(event_ind)-1:
			nextind = event_ind[i+1]
		else:
			nextind = len(event_ind)
	#	year = int(line[1:5])
	#	month = int(line[5:7])
		day = int(line[7:9])
		hour = int(line[9:11])
		minute = int(line[11:13])
		second = int(line[13:15])
		micros = int(line[15:17])*10000
	#	print(micros)
		datestr = UTCDateTime(year,month,day,hour,minute,second,micros)
		lat = line[21:24]+' '+line[24:26]+'.'+line[26:28]
		longi = line[32:36]+' '+line[36:38]+'.'+line[38:40]
		depth = line[44:49]
		if not depth.isspace():
			depth = float((line[44:47]+'.'+line[47:49]).strip())	
		else:
			depth = 0
		mag = line[52:54]
		if not mag.isspace():
		#	mag = float(mag)
	#		print('mag 1 is ' + line[52:53])
	#		print('mag 2 is ' + line[53:54])
			mag = mag[0]+'.'+mag[1]
			mag = mag.replace("-", "-0")
			mag = mag.replace("A", "-1")
		 	mag = mag.replace("B", "-2")
			mag = mag.replace("C", "-3")
			mag = float(mag)
		else:
			mag = 0
		
		cat_times.append(datestr)
		cat_mags.append(mag)
		cat_depths.append(depth)
	return cat_times, cat_mags, cat_depths
		
def readpick_from_jma(filename,station):
	filestr=filename[12:16]	
	event_ind = []
	event_intermed = []
	with open(filename) as f:
		content = f.readlines()
		for i, line in enumerate(content):
			if 'J2012' in line:
				#this is an event index
				event_ind.append(i)
				if not line[44:49].isspace():
			#	if i < 10:
				#	print(line[44:49])
				#	print(int(line[44:49]))
					if int(line[44:49]) > 6000 and int(line[44:49]) < 30000:
						if line[60] == '1' or line[60] == '5':#natural EQ or LFE
							event_intermed.append(i)
	picks_catalog = Catalog()
	stat_mags = []
	print('intermediate depth events are ' + str(len(event_intermed)))
	for i, ind in enumerate(event_intermed):
		rid = "event-"+filestr+"-"+str(i)
		new_ev = Event(resource_id=rid)
		line = content[ind]
		whichevent = event_ind.index(ind) #which event is the intermediate depth one out of all
		if whichevent < len(event_ind)-1:
			nextind = event_ind[whichevent+1]
		else:
			nextind = len(content) #if it is the last index
		year = int(line[1:5])
		month = int(line[5:7])
		day = int(line[7:9])
		hour = int(line[9:11])
		minute = int(line[11:13])
		second = float(line[13:15]+'.'+line[15:17])
		odate = UTCDateTime(year,month,day,hour,minute,second)
		sec_errs = float((line[17:19]+'.'+line[19:21]).strip())
		lat = float((line[21:24]).strip())+float((line[24:26]+'.'+line[26:28]).strip())/60
		longitude = float(line[32:36].strip())+float((line[36:38]+'.'+line[38:40]).strip())/60
		lat_err = AttribDict()
		long_err = AttribDict()
		depth_err = AttribDict()
		lat_err.err = float((line[28:30]+'.'+line[30:32]).strip())
		long_err.err = float((line[40:42]+'.'+line[42:44]).strip())
		depth_err.err = float((line[49:50]+'.'+line[50:52]).strip())
	#lat = line[21:24]+' '+line[24:26]+'.'+line[26:28]
	#	longi = line[32:36]+' '+line[36:38]+'.'+line[38:40]
		#depth in kilometres
		depth = float((line[44:47]+'.'+line[47:49]).strip())
		mag = line[52:54]
		magtype = line[54]
		if not mag.isspace():
		#	mag = float(mag)
	#		print('mag 1 is ' + line[52:53])
	#		print('mag 2 is ' + line[53:54])
			mag = mag[0]+'.'+mag[1]
			mag = float(mag)
		else:
			mag = 0
	#	print(len(mag))
	#	mag1 = line[52:53]
	#	mag2 = line[53:54]
		
#		mag = mag.strip()
#		mag=mag.lstrip("0")
#		print(mag)
#		print('len mag is ' + str(len(mag)))
#		mag = float(mag)*0.1
		station_picks = content[ind+1:nextind]
		for j, stat in enumerate(station_picks):
#check to see that both picks are there
			if stat[12] == 'h' and 'P' in stat[15:18] and 'S' in stat[28:31] and stat[3:7] == station:
				stat_mags.append(mag)
				stat_code = stat[1:6]
				ppickh = stat[19:21]
				ppickm = stat[21:23]
				ppicks = stat[23:25]+'.'+stat[25:27]
				ppicktime = UTCDateTime(int(year), int(month), int(day), int(ppickh), int(ppickm), float(ppicks))
				spickm = stat[31:33]
				spicks = stat[33:35]+'.'+stat[35:27]
				spicktime = UTCDateTime(int(year), int(month), int(day), int(ppickh), int(spickm), float(spicks))
				wavid1 = WaveformStreamID(station_code=stat_code,channel_code='UU')
				wavid2 = WaveformStreamID(station_code=stat_code,channel_code='NN')
				wavid3 = WaveformStreamID(station_code=stat_code,channel_code='EE')
				pick1 = Pick(time=ppicktime, waveform_id=wavid1, phase_hint='P')
				pick2 = Pick(time=spicktime, waveform_id=wavid2, phase_hint='S')
				pick3 = Pick(time=spicktime, waveform_id=wavid3, phase_hint='S')	
				new_ev.picks.append(pick1)
				new_ev.picks.append(pick2)
				new_ev.picks.append(pick3)
				mmag = Magnitude(mag=mag,magnitude_type=magtype)
				new_ev.magnitudes.append(mmag)
				origin = Origin()
				origin.time = odate
				origin.latitude = lat
				origin.time_errors = lat_err
				origin.longitude = longitude
				origin.latitude_errors = long_err
				origin.longitude_errors
				origin.depth = depth
				origin.depth_errors = depth_err
				new_ev.origins.append(origin)
				picks_catalog.append(new_ev)
		#		print('appended event to the pick catalog')
#	xmlname = 'picks_'+station+'_'+str(month)+'_'+str(day)+'.xml'
	xmlname=filestr+station+'.xml'
#	xmlname = 'test.xml'
	picks_catalog.write(xmlname, format="QUAKEML")
#	return stat_mags
