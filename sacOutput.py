import numpy
import datetime
import calendar
import math
import os
import obspy
import scipy
from scipy.signal import filtfilt
from scipy.signal import butter
from obspy.io.sac import SACTrace
###########################################################################
#constants
c = 299792458.0 # speed of light
fL1 = 1575.42e6 # L1 frequency
fL2 = 1227.60e6 # L2 frequency

# converts from month string to month in int
def month_converter(month):
    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    return months.index(month) + 1

# converts calendar date to doy
def doy_calc(year,month,day):
    isleap = calendar.isleap(year)
    if str(isleap) == 'True':
        dom = [31,29,31,30,31,30,31,31,30,31,30,31]
    else:
        dom = [31,28,31,30,31,30,31,31,30,31,30,31]
    doy = int(numpy.sum(dom[:(month-1)])+day)
    return(doy)

# converts from year and day to gps week and gps day of week
def gpsweekdow(year,doy):
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)
    gpstime = (numpy.datetime64(date) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's')
    gpsweek = int(gpstime/604800)
    gpsdow = math.floor((gpstime-gpsweek*604800)/86400)
    return(gpsweek, gpsdow)

# converts all of this to gps time
def gpstimeconvert(gpstime):
    gpsweek = int(gpstime/604800)
    gpsdow = math.floor((gpstime-gpsweek*604800)/86400)
    gpssow = gpstime-gpsweek*604800
    return(gpsweek, gpsdow, gpssow)



#This takes displacements in x, y, z and converts them to north, east up
def dxyz2dneu(dx,dy,dz,lat,lon):
    lat = lat*math.pi/180
    lon = lon*math.pi/180
    dn = -numpy.sin(lat)*numpy.cos(lon)*dx-numpy.sin(lat)*numpy.sin(lon)*dy+numpy.cos(lat)*dz
    de = -numpy.sin(lon)*dx+numpy.cos(lon)*dy
    du = numpy.cos(lat)*numpy.cos(lon)*dx+numpy.cos(lat)*numpy.sin(lon)*dy+numpy.sin(lat)*dz
    return (dn, de, du)

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

def writesac(velfile, site, stalat, stalon, doy, year, samprate, event):
    a = numpy.loadtxt(velfile)
    tind = a[:, 0]
    gtime = a[:, 1]
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
    nunf = a[:, 2] - numpy.nanmean(a[:, 2])
    eunf = a[:, 3] - numpy.nanmean(a[:, 3])
    uunf = a[:, 4] - numpy.nanmean(a[:, 4])
    t = gtime - gtime[0] # subtracts the first epoch
    print(samprate) # check sample rate
    samplerate = float(samprate) # convert to float
    # butterworth digital and analog filter design
    # numerator, denominator = butter(order of filter, critical frequencies, low pass)
    bf, af = butter(4, 1.25 / 0.5 * samplerate, btype='low')
    nv = filtfilt(bf, af, nunf) # applies the filter forward and backward to the signal
    ev = filtfilt(bf, af, eunf)
    uv = filtfilt(bf, af, uunf)
    plotvelocities(event, site, t, nv * 100, ev * 100, uv * 100)
    sr = "{0:.2f}".format(float(samprate))

    #output LXN
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXN.sac')
    headN = {'kstnm': site, 'kcmpnm': 'LXN', 'stla': float(stalat), 'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sacn = SACTrace(data=nv, **headN)
    sacn.write('output/sacfilt/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.filt.LXN.sac')
    sacn = SACTrace(data=nunf, **headN)
    sacn.write('output/sac/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.LXN.sac')

    #output LXE
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXE.sac')
    headE = {'kstnm': site, 'kcmpnm': 'LXE', 'stla': float(stalat), 'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sace = SACTrace(data=ev, **headE)
    sace.write('output/sacfilt/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.filt.LXE.sac')
    sace = SACTrace(data=eunf, **headE)
    sace.write('output/sac/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.LXE.sac')

    #output LXZ
    print('Writing SAC file ' + 'output/' + event + '.' + site + '.' + sr + '.LXZ.sac')
    headZ = {'kstnm': site, 'kcmpnm': 'LXZ', 'stla': float(stalat), 'stlo': float(stalon),
             'nzyear': int(year), 'nzjday': int(doy), 'nzhour': int(sthr), 'nzmin': int(stmin),
             'nzsec': int(stsec), 'nzmsec': int(0), 'delta': float(samprate)}
    sacu = SACTrace(data=uv, **headZ)
    sacu.write('output/sacfilt/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.filt.LXZ.sac')
    sacu = SACTrace(data=uunf, **headZ)
    sacu.write('output/sac/' + event + '/' + event + '.' + site.upper() + '.' + sr + '.LXZ.sac')

    # writesac('output/velocities_p156_354_2021.txt','p156','35','-120','354','2021',0.2,'ferndale')