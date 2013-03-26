#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 14:34:54 2012

@author: mag
"""

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 7, 5)
__modified__ = datetime.datetime(2012, 7, 5)
__version__  = "1.0"
__status__   = "Development"

# Adapted from code & formulas by David Z. Creemer and others
# http://www.zachary.com/blog/2005/01/12/python_zipcode_geo-programming
# http://williams.best.vwh.net/avform.htm
#

from math import sin,cos,atan,acos,asin,atan2,sqrt,pi, modf

# At the equator / on another great circle???
nauticalMilePerLat = 60.00721
nauticalMilePerLongitude = 60.10793

rad = pi / 180.0

milesPerNauticalMile = 1.15078
kmsPerNauticalMile = 1.85200

degreeInMiles = milesPerNauticalMile * 60
degreeInKms = kmsPerNauticalMile * 60

# Semi-axes of WGS-84 geoidal reference
WGS84_a = 6378.1370  # Major semiaxis [km]
WGS84_b = 6356.7523  # Minor semiaxis [km]


# Earth radius at a given latitude, according to the WGS-84 ellipsoid [m]
def WGS84EarthRadius(lat):
    # http://en.wikipedia.org/wiki/Earth_radius
    lat = lat * rad
    An = WGS84_a*WGS84_a * cos(lat)
    Bn = WGS84_b*WGS84_b * sin(lat)
    Ad = WGS84_a * cos(lat)
    Bd = WGS84_b * sin(lat)
    return sqrt( (An*An + Bn*Bn)/(Ad*Ad + Bd*Bd) )


def getDistanceLL(loc1, loc2):
   "Haversine formula - give coordinates as (lat_decimal,lon_decimal) tuples"
   lat1, lon1 = loc1
   lat2, lon2 = loc2
   # earth's mean radius = 6,371km
   earthradius = WGS84EarthRadius((lat1+lat2)/2)
   #if type(loc1[0]) == type(()):
   #   # convert from DMS to decimal
   #   lat1,lon1 = DMSToDecimal(loc1[0]),DMSToDecimal(loc1[1])
   #if type(loc2[0]) == type(()):
   #   lat2,lon2 = DMSToDecimal(loc2[0]),DMSToDecimal(loc2[1])
   # convert to radians
   lon1 = lon1 * rad
   lon2 = lon2 * rad
   lat1 = lat1 * rad
   lat2 = lat2 * rad
   # haversine formula
   dlon = lon2 - lon1
   dlat = lat2 - lat1
   a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2.0))**2
   c = 2.0 * atan2(sqrt(a), sqrt(1.0-a))
   km = earthradius * c
   return km

def getPixelResolution(lat, lon, shape, units="km"):
    """
    get the resolution of the pixel in Row/Col directions
    """
 
    if str(lat.size) == '4':
       col = ([lat[0], lon[0]], [lat[2], lon[2]])
       row = ([lat[0], lon[0]], [lat[3], lon[3]])
    else:
       col = ([lat[0,0], lon[0,0]], [lat[0,-1], lon[0,-1]])
       row = ([lat[0,0], lon[0,0]], [lat[-1,0], lon[-1,0]])

    PxResRow = getDistanceLL(col[0], col[1])/shape[1]
    PxResCol = getDistanceLL(row[0], row[1])/shape[0]

    if units == "deg":
        PxResRow = min(getCoordinateDiffForDistance(lat[0], lon[0], PxResRow, units="km"))
        PxResCol = min(getCoordinateDiffForDistance(lat[0], lon[0], PxResCol, units="km"))

    return PxResCol, PxResRow
    
def getDistancePx(loc1, loc2, lat, lon, shape):
    """Get the distance between two points given coordinates in loc[row,col]"""
    PxResCol, PxResRow = getPixelResolution(lat, lon, shape)
    km = sqrt( (PxResCol*(loc1[0]-loc2[0]))**2 + (PxResRow*(loc1[1]-loc2[1]))**2 )
    return km
    
def DecimalToDMS(decimalvalue):
    "convert a decimal value to degree,minute,second tuple"
    m,s = divmod(dd*3600,60)
    d,m = divmod(mnt,60)
#   d = modf(decimalvalue)[0]
#   m=0
#   s=0
    return (d,m,s)


def DMSToDecimal((degrees,minutes,seconds)):
   "Convert a value from decimal (float) to degree,minute,second tuple"
   d = abs(degrees) + (minutes/60.0) + (seconds/3600.0)
   if degrees < 0:
      return -d
   else:
      return d


def getCoordinateDiffForDistance(originlat, originlon, distance, units="km"):
   """return longitude & latitude values that, when added to & subtraced from
   origin longitude & latitude, form a cross / 'plus sign' whose ends are
   a given distance from the origin"""

   degreelength = 0

   if units == "km":
      degreelength = degreeInKms
   elif units == "miles":
      degreelength = degreeInMiles
   else:
      raise Exception("Units must be either 'km' or 'miles'!")

   lat = distance / degreelength
   lon = distance / (cos(originlat * rad) * degreelength)

   return (lat, lon)


def isWithinDistance(origin, loc, distance):
   "boolean for checking whether a location is within a distance"
   if getDistance(origin, loc) <= distance:
      return True
   else:
      return False
