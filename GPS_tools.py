#!/home/crowellb/anaconda3/bin/python3
import numpy
import datetime
import calendar
import math
import os

#####################################################################################
# Constants
c = 299792458.0  # speed of light
fL1 = 1575.42e6  # L1 frequency
fL2 = 1227.60e6  # L2 frequency


def month_converter(month):
    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    return months.index(month) + 1


def doy_calc(year, month, day):
    isleap = calendar.isleap(year)
    if str(isleap) == 'True':
        dom = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy = int(numpy.sum(dom[:(month - 1)]) + day)
    return (doy)


def gpsweekdow(year, doy):
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)
    gpstime = (numpy.datetime64(date) - numpy.datetime64('1980-01-06T00:00:00')) / numpy.timedelta64(1, 's')
    gpsweek = int(gpstime / 604800)
    gpsdow = math.floor((gpstime - gpsweek * 604800) / 86400)
    return (gpsweek, gpsdow)


def gpstimeconvert(gpstime):
    gpsweek = int(gpstime / 604800)
    gpsdow = math.floor((gpstime - gpsweek * 604800) / 86400)
    gpssow = gpstime - gpsweek * 604800
    return (gpsweek, gpsdow, gpssow)


def lla2ecef(lat, lon, alt):
    lat = lat * math.pi / 180
    lon = lon * math.pi / 180
    a = 6378137
    e = 8.1819190842622e-2

    N = a / numpy.sqrt(1 - numpy.power(e, 2) * numpy.power(numpy.sin(lat), 2))

    x = (N + alt) * numpy.cos(lat) * numpy.cos(lon)
    y = (N + alt) * numpy.cos(lat) * numpy.sin(lon)
    z = ((1 - numpy.power(e, 2)) * N + alt) * numpy.sin(lat)

    return (x, y, z)


def ecef2lla(x, y, z):
    a = 6378137
    e = 8.1819190842622e-2
    b = math.sqrt(math.pow(a, 2) * (1 - math.pow(e, 2)))
    ep = math.sqrt((math.pow(a, 2) - math.pow(b, 2)) / math.pow(b, 2))
    p = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    th = math.atan2(a * z, b * p)
    lon = math.atan2(y, x)
    lat = math.atan2((z + math.pow(ep, 2) * b * math.pow(math.sin(th), 3)),
                     (p - math.pow(e, 2) * a * math.pow(math.cos(th), 3)))
    N = a / math.sqrt(1 - math.pow(e, 2) * math.pow(math.sin(lat), 2))
    alt = p / math.cos(lat) - N
    return (lat, lon, alt)


def azi_elev(xsta, ysta, zsta, xsat, ysat, zsat):
    [latsta, lonsta, altsta] = ecef2lla(xsta, ysta, zsta)
    [latsat, lonsat, altsat] = ecef2lla(xsat, ysat, zsat)
    re = math.sqrt(math.pow(xsta, 2) + math.pow(ysta, 2) + math.pow(zsta, 2))
    rs = math.sqrt(math.pow(xsat, 2) + math.pow(ysat, 2) + math.pow(zsat, 2))
    gamma = math.acos(
        math.cos(latsta) * math.cos(latsat) * math.cos(lonsat - lonsta) + math.sin(latsta) * math.sin(latsat))
    elev = math.acos(math.sin(gamma) / math.sqrt(1 + math.pow(re / rs, 2) - 2 * re / rs * math.cos(gamma)))

    deltalon = lonsat - lonsta

    azi = math.atan2(math.sin(deltalon) * math.cos(latsat),
                     math.cos(latsta) * math.sin(latsat) - math.sin(latsta) * math.cos(latsat) * math.cos(deltalon))

    azi = azi * 180 / math.pi

    if (azi < 0):
        azi = azi + 360
    elev = elev * 180 / math.pi
    return (azi, elev)


# This takes displacements in x, y, z and converts them to north, east up

def dxyz2dneu(dx, dy, dz, lat, lon):
    lat = lat * math.pi / 180
    lon = lon * math.pi / 180
    dn = -numpy.sin(lat) * numpy.cos(lon) * dx - numpy.sin(lat) * numpy.sin(lon) * dy + numpy.cos(lat) * dz
    de = -numpy.sin(lon) * dx + numpy.cos(lon) * dy
    du = numpy.cos(lat) * numpy.cos(lon) * dx + numpy.cos(lat) * numpy.sin(lon) * dy + numpy.sin(lat) * dz
    return (dn, de, du)


# This takes covariances in x, y, z and converts them to north, east up

def covrot(cx, cy, cz, lat, lon):
    lat = lat * math.pi / 180
    lon = lon * math.pi / 180
    slat = numpy.sin(lat)
    slon = numpy.sin(lon)
    clat = numpy.cos(lat)
    clon = numpy.cos(lon)
    R = numpy.array([[-slon, -slat * clon, clat * clon], [clon, -slat * slon, clat * slon], [0, clat, slat]])
    Pxyz = numpy.array([[cx * cx, 0, 0], [0, cy * cy, 0], [0, 0, cz * cz]])
    Pneu = numpy.transpose(R) * Pxyz * R
    cn = numpy.sqrt(Pneu[1, 1])
    ce = numpy.sqrt(Pneu[0, 0])
    cu = numpy.sqrt(Pneu[2, 2])
    return (cn, ce, cu)


# Time Series to CSV writer
def ts2csv(site, tdir, csvdir, startyear, endyear):
    csvfile = csvdir + '/' + site + '.csv'
    fno = open(csvfile, 'w')
    x0 = 0
    y0 = 0
    z0 = 0
    dncorr = 0
    decorr = 0
    ducorr = 0
    ind = 0
    for y in range(int(startyear), int(endyear) + 1):
        isleap = calendar.isleap(int(y))
        ystr = str(y)
        if str(isleap) == 'True':
            days = 366
        else:
            days = 365
        for i in range(0, days + 1):
            dy = str(i).zfill(3)
            fnameo = tdir + '/' + site + '_' + dy + '_' + ystr[-2:] + '.txt'
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
                            [lat, lon, alt] = ecef2lla(x0, y0, z0)
                            lat = lat * 180.0 / math.pi
                            lon = lon * 180.0 / math.pi
                        [ner, eer, uer] = covrot(float(grow[5]), float(grow[6]), float(grow[7]), lat, lon)
                        [dnorth, deast, dup] = dxyz2dneu(float(grow[2]) - x0, float(grow[3]) - y0, float(grow[4]) - z0,
                                                         lat, lon)
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
                wn = 1 / numpy.power(nearray, 2)
                we = 1 / numpy.power(eearray, 2)
                wu = 1 / numpy.power(uearray, 2)
                dnweight = numpy.sum(dnarray * wn) / numpy.sum(wn)
                deweight = numpy.sum(dearray * we) / numpy.sum(we)
                duweight = numpy.sum(duarray * wu) / numpy.sum(wu)
                # dnweight = numpy.median(dnarray)
                # deweight = numpy.median(dearray)
                # duweight = numpy.median(duarray)
                neweight = (numpy.sum(numpy.power(dnarray, 2) * wn) / numpy.sum(wn) - numpy.power(dnweight,
                                                                                                  2)) * math.pow(
                    numpy.sum(wn), 2) / (math.pow(numpy.sum(wn), 2) - numpy.sum(numpy.power(wn, 2)))
                eeweight = (numpy.sum(numpy.power(dearray, 2) * we) / numpy.sum(we) - numpy.power(deweight,
                                                                                                  2)) * math.pow(
                    numpy.sum(we), 2) / (math.pow(numpy.sum(we), 2) - numpy.sum(numpy.power(we, 2)))
                ueweight = (numpy.sum(numpy.power(duarray, 2) * wu) / numpy.sum(wu) - numpy.power(duweight,
                                                                                                  2)) * math.pow(
                    numpy.sum(wu), 2) / (math.pow(numpy.sum(wu), 2) - numpy.sum(numpy.power(wu, 2)))
                if (ind == 0):
                    n1 = dnweight
                    e1 = deweight
                    u1 = duweight
                    ind = ind + 1
                os.system('/home/crowellb/GipsyX-1.3/bin/sec2date ' + str(numpy.mean(t)) + ' > /home/crowellb/tmp.txt')
                avgtime = numpy.loadtxt('/home/crowellb/tmp.txt', dtype=str)
                # print(avgtime[0], avgtime[1])
                if (dnweight < 100):
                    nstr = "{0:.2f}".format((dnweight - n1) * 1000 + 50)
                    estr = "{0:.2f}".format((deweight - e1) * 1000)
                    ustr = "{0:.2f}".format((duweight - u1) * 1000 - 50)
                    nestr = "{0:.2f}".format(math.sqrt(neweight) * 1000 / 2)
                    eestr = "{0:.2f}".format(math.sqrt(eeweight) * 1000 / 2)
                    uestr = "{0:.2f}".format(math.sqrt(ueweight) * 1000 / 2)
                    fno.write(avgtime[
                                  0] + ' 12:00:00' + ',' + nstr + ',' + nestr + ',' + estr + ',' + eestr + ',' + ustr + ',' + uestr + '\n')
    fno.close()
    return

# ts2csv('sc04','/gnssarchive3/daily','sscsv',2021,2022)




