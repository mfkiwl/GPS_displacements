import numpy
import datetime
import calendar
import math
import os
import sys
import obspy
import scipy
from scipy.signal import filtfilt
from scipy.signal import butter
from obspy.io.sac import SACTrace
sys.path.insert(1, '/home/jdegran/jd_gipsy_working')
from SNIVEL_tools import ecef2lla
from GPS_tools import covrot
###########################################################################
# for sac file structure for velocities processed with SNIVEL see, geodesy.ess.washington.edu/SNIVEL/output
# this script is taken an edited from SNIVEL_tools.py written by Brendan Crowell

# Constants ###############################################################

c = 299792458.0 # speed of light
fL1 = 1575.42e6 # L1 frequency
fL2 = 1227.60e6 # L2 frequency

# Functions ##########################################################################################################
#This takes displacements in x, y, z and converts them to north, east up
def dxyz2dneu(dx,dy,dz,lat,lon):
    lat = lat*math.pi/180
    lon = lon*math.pi/180
    dn = -numpy.sin(lat)*numpy.cos(lon)*dx-numpy.sin(lat)*numpy.sin(lon)*dy+numpy.cos(lat)*dz
    de = -numpy.sin(lon)*dx+numpy.cos(lon)*dy
    du = numpy.cos(lat)*numpy.cos(lon)*dx+numpy.cos(lat)*numpy.sin(lon)*dy+numpy.sin(lat)*dz
    return (dn, de, du)

def gpsleapsec(gpssec):
    leaptimes = numpy.array([46828800, 78364801, 109900802, 173059203, 252028804, 315187205, 346723206, 393984007, 425520008, 457056009, 504489610, 551750411, 599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017])
    leapseconds = numpy.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])
    a1 = numpy.where(gpssec > leaptimes)[0]
    leapsec = len(a1)
    return(leapsec)

#Time Series to CSV writer
def ts2csv(site,tdir,csvdir,startyear):
    csvfile = csvdir + '/' + site + '.csv'
    fno = open(csvfile,'w')
    isleap = calendar.isleap(int(startyear))
    ystr = str(startyear)
    x0 = 0
    y0 = 0
    z0 = 0
    dncorr=0
    decorr=0
    ducorr=0
    ind=0
    if str(isleap) == 'True':
        days = 366
    else:
        days = 365
    for i in range(0,days):
        dy = str(i).zfill(3)
        fnameo = tdir+'/'+site+'_'+dy+'_'+ystr[-2:]+'.txt'
        if (os.path.isfile(fnameo) == True):
            tl = list()
            dn = list()
            de = list()
            du = list()
            ne = list()
            ee = list()
            ue = list()
            with open(fnameo, 'rt') as g:
                rows = (line.split() for line in g)
                for grow in rows:
                    if (x0 == 0):
                        x0 = float(grow[2])
                        y0 = float(grow[3])
                        z0 = float(grow[4])
                        [lat,lon,alt]=ecef2lla(x0,y0,z0)
                        lat = lat*180.0/math.pi
                        lon = lon*180.0/math.pi
                    [ner,eer,uer]=covrot(float(grow[5]),float(grow[6]),float(grow[7]),lat,lon)
                    [dnorth,deast,dup]=dxyz2dneu(float(grow[2])-x0,float(grow[3])-y0,float(grow[4])-z0,lat,lon)
                    ne.append(ner)
                    ee.append(eer)
                    ue.append(uer)
                    dn.append(dnorth)
                    de.append(deast)
                    du.append(dup)
                    tl.append(float(grow[0]))
            dnarray = numpy.asarray(dn)
            dearray = numpy.asarray(de)
            duarray = numpy.asarray(du)
            nearray = numpy.asarray(ne)
            eearray = numpy.asarray(ee)
            uearray = numpy.asarray(ue)
            t = numpy.asarray(tl)
            wn = 1/numpy.power(nearray,2)
            we = 1/numpy.power(eearray,2)
            wu = 1/numpy.power(uearray,2)
            dnweight = numpy.sum(dnarray*wn)/numpy.sum(wn)
            deweight = numpy.sum(dearray*we)/numpy.sum(we)
            duweight = numpy.sum(duarray*wu)/numpy.sum(wu)
            #dnweight = numpy.median(dnarray)
            #deweight = numpy.median(dearray)
            #duweight = numpy.median(duarray)
            neweight = (numpy.sum(numpy.power(dnarray,2)*wn)/numpy.sum(wn) - numpy.power(dnweight,2))*math.pow(numpy.sum(wn),2)/(math.pow(numpy.sum(wn),2) - numpy.sum(numpy.power(wn,2)))
            eeweight = (numpy.sum(numpy.power(dearray,2)*we)/numpy.sum(we) - numpy.power(deweight,2))*math.pow(numpy.sum(we),2)/(math.pow(numpy.sum(we),2) - numpy.sum(numpy.power(we,2)))
            ueweight = (numpy.sum(numpy.power(duarray,2)*wu)/numpy.sum(wu) - numpy.power(duweight,2))*math.pow(numpy.sum(wu),2)/(math.pow(numpy.sum(wu),2) - numpy.sum(numpy.power(wu,2)))
            if (ind == 0):
                n1 = dnweight
                e1 = deweight
                u1 = duweight
                ind=ind+1
            os.system('/home/crowellb/GipsyX-1.3/bin/sec2date '+ str(numpy.mean(t)) + ' > /home/crowellb/tmp.txt')
            avgtime = numpy.loadtxt('/home/crowellb/tmp.txt',dtype=str)
            #print(avgtime[0], avgtime[1])
            if (dnweight < 100):
                nstr = "{0:.2f}".format((dnweight-n1)*1000+50)
                estr = "{0:.2f}".format((deweight-e1)*1000)
                ustr = "{0:.2f}".format((duweight-u1)*1000-50)
                nestr = "{0:.2f}".format(math.sqrt(neweight)*1000/2)
                eestr = "{0:.2f}".format(math.sqrt(eeweight)*1000/2)
                uestr = "{0:.2f}".format(math.sqrt(ueweight)*1000/2)
                fno.write(avgtime[0]+' 12:00:00'+','+nstr+','+nestr+','+estr+','+eestr+','+ustr+','+uestr+'\n')
    fno.close()
    return
#ts2csv('sc04','/gnssarchive/daily','sscsv',2021)

def writesac(dispfile, site, stalat, stalon, doy, year, samprate, event):
    a = numpy.loadtxt(dispfile)
    #tind = a[:, 0]
    gtime = a[:, 0]
    leapsec = gpsleapsec(gtime[0]) + 630763200

    # Get the start time of the file in UTC
    date = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(doy) - 1)
    gpstime = (numpy.datetime64(date) - numpy.datetime64('1980-01-06T00:00:00')) / numpy.timedelta64(1, 's')
    stime = (gtime[0] - leapsec) * numpy.timedelta64(1, 's') + numpy.datetime64('2000-01-01T12:00:00')
    sitem = stime.item()
    print(sitem)
    styr = sitem.year
    stdy = sitem.day
    stmon = sitem.month
    sthr = sitem.hour
    stmin = sitem.minute
    stsec = sitem.second

    # subtracts the average
    nunf = a[:, 1] - numpy.nanmean(a[:, 1]) #n unfiltered
    eunf = a[:, 2] - numpy.nanmean(a[:, 2])
    uunf = a[:, 3] - numpy.nanmean(a[:, 3])
    t = gtime - gtime[0] # subtracts the first epoch
    print(samprate) # check sample rate
    samplerate = float(samprate) # convert to float
    # butterworth digital and analog filter design
    # numerator, denominator = butter(order of filter, critical frequencies, low pass)
    bf, af = butter(4, 1.25 / 0.5 * samplerate, btype='low')
    nv = filtfilt(bf, af, nunf) # applies the filter forward and backward to the signal
    ev = filtfilt(bf, af, eunf)
    uv = filtfilt(bf, af, uunf)
    #plotvelocities(event, site, t, nv * 100, ev * 100, uv * 100)
    sr = "{0:.2f}".format(float(samprate))

    #output LXN
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXN.sac')
    headN = {'kstnm': site, 'kcmpnm': 'LXN', 'stla': float(stalat), 'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sacn = SACTrace(data=nv, **headN)
    if not os.path.exists(mainDir + 'output/sacfilt/' + event):
        os.makedirs(mainDir + 'output/sacfilt/' + event)
    sacn.write(mainDir + 'output/sacfilt/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.filt.LXN.sac')
    #sacn.write(event + '.' + site.upper() + '.' + sr + '.filt.LXN.sac')

    sacn = SACTrace(data=nunf, **headN)
    if not os.path.exists(mainDir + 'output/sac/' + event):
        os.makedirs(mainDir + 'output/sac/' + event)
    sacn.write(mainDir + 'output/sac/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.LXN.sac')
    #sacn.write(event + '.' + site.upper() + '.' + sr + '.LXN.sac')
    #### command on the above line works and creates the correct file in the cwd -- need to figure out how to write to NOT current directory

    #output LXE
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXE.sac')
    headE = {'kstnm': site, 'kcmpnm': 'LXE', 'stla': float(stalat), 'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sace = SACTrace(data=ev, **headE)
    sace.write(mainDir + 'output/sacfilt/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.filt.LXE.sac')
    sace = SACTrace(data=eunf, **headE)
    sace.write(mainDir + 'output/sac/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.LXE.sac')

    #output LXZ
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXZ.sac')
    headZ = {'kstnm': site, 'kcmpnm': 'LXZ', 'stla': float(stalat), 'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sacu = SACTrace(data=uv, **headZ)
    sacu.write(mainDir + 'output/sacfilt/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.filt.LXZ.sac')
    sacu = SACTrace(data=uunf, **headZ)
    sacu.write(mainDir + 'output/sac/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.LXZ.sac')

    # writesac('output/velocities_p156_354_2021.txt','p156','35','-120','354','2021',0.2,'ferndale')

# Actual Script #####################################################################################################
# first thing we need to do is take the input file (output file from GipsyX)
#    example: ac12_210_full.txt
#    J time                 POS              X                 Y                Z           dx        dy        dz
#    680810400.0000000    GPSPOS       -3450877.1679    -1284085.3131     5190636.3242    0.0168    0.0110    0.0187
#    680810401.0000000    GPSPOS       -3450877.1603    -1284085.3096     5190636.3184    0.0168    0.0110    0.0187
#    680810402.0000000    GPSPOS       -3450877.1644    -1284085.3076     5190636.3203    0.0168    0.0110    0.0187
#    680810403.0000000    GPSPOS       -3450877.1638    -1284085.3092     5190636.3101    0.0168    0.0110    0.0167
#    680810404.0000000    GPSPOS       -3450877.1626    -1284085.3073     5190636.3038    0.0168    0.0110    0.0167
#    680810405.0000000    GPSPOS       -3450877.1645    -1284085.3062     5190636.3124    0.0168    0.0110    0.0167
#    680810406.0000000    GPSPOS       -3450877.1626    -1284085.3052     5190636.3147    0.0168    0.0110    0.0167
#    680810407.0000000    GPSPOS       -3450877.1665    -1284085.3057     5190636.3219    0.0168    0.0110    0.0167
#    680810408.0000000    GPSPOS       -3450877.1638    -1284085.3046     5190636.3278    0.0176    0.0112    0.0168
#    680810409.0000000    GPSPOS       -3450877.1620    -1284085.3030     5190636.3305    0.0176    0.0112    0.0168

######################################################################################################################
# get user entered inputs
print('running sacOutput...')
print ("Enter the event name")
event = input()
print("enter the year of the event (yyyy)")
year = input()
print("enter the doy (ddd)")
doy = input()

# hard-coded variables
samprate = 0.2

# set the datapaths
mainDir = '/home/jdegran/jd_gipsy_working/'
eventDir = mainDir + event + '/'

#Create a variable for the file name containing all the stations that were downloaded and processed through GipsyX
filename = eventDir + 'downloadedSites_' + event + '_latlon.txt'

site = list()
siteLat = numpy.array([])
siteLon = numpy.array([])
#staAlt = np.array([])

#Open the file, read into each vector
infile = open(filename, 'r')
lines = infile.readlines()
for line in lines:
    sline = line.split(',')
    site.append(sline[0])
    siteLat = numpy.append(siteLat, float(sline[1]))
    siteLon = numpy.append(siteLon, float(sline[2]))
    #staAlt = numpy.append(staAlt, float(sline[3]))

infile.close()  #Always close the file!







# this will become a loop at some point, to loop through a file with all this information?
# write now let's do it manually
# site = 'ac12'
# doy = '210'
# year ='2021'
# siteLat = 54.83097 # station latitude --- could even pull these from numSites with the little script written earlier
# siteLon = -159.58954 # station longitude



# eventually this will be a directory to my ess.geodesy home, but testing in local python env
#filename = '/Users/jensen/PycharmProjects/GPS_displacements/' + site + '_' + doy + '_full.txt'
#filename = '/Users/jensendegrande/PycharmProjects/GPS_displacements/' + site + '_' + doy + '_full.txt'
#filename = '/home/jdegran/jd_gipsy_working/Alaska75/' + site + '_' + doy + '_full.txt'


# okay now we can make it a looooooop
numSites = len(site)
counter = 0
for i in range(0, numSites):
    filename = eventDir + site[i] + '_' + doy + '_full.txt' # this is the output file from gipsyX
    # load said file
    if os.stat(filename).st_size != 0:
        data = numpy.loadtxt(filename, dtype=float, usecols=(0,2,3,4,5,6,7))

        time = data[:,0]
        x = data[:,1]
        y = data[:,2]
        z = data[:,3]
        dx = data[:,4]
        dy = data[:,5]
        dz = data[:,6]

        # convert from x,y,z to n,e,u, use lat and lon of station (? double check this with brendan)
        [n, e, u] = dxyz2dneu(x,y,z,siteLat[i],siteLon[i])
        [dn, de, du] = dxyz2dneu(dx,dy,dz,siteLat[i],siteLon[i])

        # then create displacements file (like the velcoities files in SNIVEL) - a file with the time (JTime, dn, de, and du)
        dispData = numpy.column_stack((time,dn,de,du)) # write this to an actual file

        outFile = eventDir + 'displacements_' + site[i] + '_' + doy + '.txt' #outdirdata + '/' + SID + '_' + str(eqtime) + '.txt'
        print("Printing file " + outFile)
        output = open(outFile,'w+')

        length = len(dispData)
        # write the results of the conversion from x,y,z to e,n,u to the displacements file
        for j in range(0, length):
                timer = time
                north = dn
                east = de
                up = du
                output.write(str(j) + ' ' + str(timer[j]) + ' ' + str(north[j]) + ' ' + str(east[j]) + ' ' + str(up[j]) + '\n')

        output.close() # cool now we can write to sac format

        writesac(outFile, site[i], siteLat[i], siteLon[i], doy, year, samprate, event)
    else:
        print('file ' + filename + ' is empty.. moving on')
        counter+=1

print(str(counter) + ' event was not processed because file was empty')

# run writeSAC
# example call to function
# writesac('output/velocities_p156_354_2021.txt','p156','35','-120','354','2021',0.2,'ferndale')

