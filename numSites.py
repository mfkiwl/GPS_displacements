# Quick little script to figure out what stations to include in GipsyX
# 1. Input earthquake event location in lat, lon, alt
# 2. Convert to x,y,z with lla2ecef
# 3. Input station file of all unavco stations
# 4. Convert to x,y,z
# 5. Compute distance between event location and station
# 6. If distance is <100 km == output the station name to file name

import numpy as np
import datetime
import calendar
import math
# import functions from SNIVEL_tools
from SNIVEL_tools import *

#Create a variable for the file name containing all the stations
filename = "stations.txt"

staName = list()
staLat = np.array([])
staLon = np.array([])
staAlt = np.array([])

#Open the file, read into each vector
infile = open(filename, 'r')
lines = infile.readlines()
for line in lines:
    sline = line.split(' ')
    staName.append(sline[0])
    staLat = np.append(staLat, float(sline[1]))
    staLon = np.append(staLon, float(sline[2]))
    staAlt = np.append(staAlt, float(sline[3]))

infile.close()  #Always close the file!


# earthquake event- eventually make the user enter this
# sand point, Alaska
eqLat = 54.602
eqLon = -159.626
alt = -28.4

# lla2ecef(lat,lon,alt):
[eqX, eqY, eqZ] = lla2ecef(eqLat,eqLon,alt) #outputs in meters!!!
eqX = eqX/1000
eqY = eqY/1000
eqZ = eqZ/1000 # convert to km

# unavco stations
# staName = ['ac01', 'test','abcd'] # vector of strings of unavco station names
# staLat = np.array([40.001, 100.000, -3.00])# some vector of latitudes for all stations
# staLon = np.array([-159.200, -138, 48])# some vector of longitudes for all stations
# alt = np.array([5, 24, 78]) # some vector for all latitudes

staLat = np.transpose(staLat)
staLon = np.transpose(staLon)
staAlt = np.transpose(staAlt)


# fun output file for the loop
numSites = open('sites.txt','w')

# a loooooop
counter = 0
numSta = len(staLat)
for i in range(0, numSta, 1):
    # lla2ecef(lat,lon,alt):
    [staX, staY, staZ] = lla2ecef(staLat, staLon, alt)

    staX = staX/1000
    staY = staY/1000
    staZ = staZ/1000

    dx = eqX - staX[i]
    dy = eqY - staY[i]
    dz = eqZ - staZ[i]

    # distance from eq to station
    dist = np.sqrt(dx**2 + dy**2 + dz**2)

    if (dist <= 300): # 100 km - this will change dependent on magnitude of eq
        #print(staName[i])
        numSites.write(staName[i]+ '\n')
        counter+=1

print(counter, 'many stations have been found within the radius')
print('no more stations within 300 km of event')
numSites.close()
# include something for if rinex exisyt = use station, if rinex does not exist, discard