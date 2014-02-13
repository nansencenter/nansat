#!/usr/bin/env python
# Name:    domain_test.py
# Purpose: Tutorial for Domain class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import os

from nansat import Domain

# input and output file names
from testio import testio
iPath, oPath = testio()
iFileName = os.path.join(iPath, 'gcps.tif')
print 'Input file: ', iFileName, '\n'
oFileName = os.path.join(oPath, 'output_domain_')
print 'Output file:', oFileName, '\n'

''' Domain class is a container for geographical reference of a raster

A Domain object describes all attributes of geographical
reference of a raster:
    width and height, pixel size,
    relation between pixel/line coordinates and geographical coordinates,
    type of data projection (e.g. geographical or stereographic)

'''

# Create Domain object. It describes the desired grid of reprojected image:
# projection, resolution, size, etc. In this case it is geographic projection;
# -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
           "-te 25 70 35 72 -ts 2000 2000")
d = Domain(4326, "-te 25 70 35 72 -ts 2000 2000")
d.write_map(oFileName + '01_latlong_map.png')
print 'Latlong Domain:', d, '\n'

# write to KML
#d.write_kml(kmlFileName=oFileName + '01_latlong_map.kml')

# read from KML
#d = Domain(oFileName + '01_latlong_map.kml')
#print 'Domain from KML', d

# get shape
print 'shape:', d.shape(), '\n'

# Generate two vectors with values of lat/lon for the border of domain
lonVec, latVec = d.get_border()
print 'lonVec :', lonVec, '\n'
print 'latVec :', latVec, '\n'

# Get upwards azimuth direction of domain.
bearing_center = d.upwards_azimuth_direction()
print 'bearing_center :', bearing_center, '\n'

# Create domain with stereographic projection
# -- Get corners of the image
lons, lats = d.get_corners()
meanLon = sum(lons, 0.0) / 4.
meanLat = sum(lats, 0.0) / 4.

# get string with WKT representation of the border polygon
print 'BorderPolygon:', d.get_border_wkt(), '\n'

# longitude and latitude grids representing the full data grid
longitude, latitude = d.get_geolocation_grids()
print 'longitude shape: :', longitude.shape, '\n'

# make KML file with image borders (to be opened in Googe Earth)
d.write_kml(kmlFileName=oFileName + '03_preview.kml')

# -- Create Domain with stereographic projection and resolution 1000m
srsString = "+proj=stere +lon_0=5 +lat_0=60"
extentString = '-te -1000000 -1000000 1000000 1000000 -tr 1000 1000'
d = Domain(srsString, extentString)
d.write_map(oFileName + '02_stereo_map.png')
print 'Stereo Domain:', d, '\n'


print '\n*** domain_test completed successfully. Output files are found here:' + oFileName

