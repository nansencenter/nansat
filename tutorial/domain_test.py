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

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import inspect, os

from nansat import Nansat, Domain, Mosaic
#from mosaic import Mosaic

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
def main():
    # Open an input file
    # Create a Nansat object <n> for futher high-level operations
    n = Nansat(iFileName)

    # Open an input file, specify which Mapper to use, set logging level
    n = Nansat(iFileName, mapperName='generic', logLevel=20)

    # get shape
    print 'shape:', n.shape(), '\n'

    # get lists with coordinates of the object borders
    print 'Border:', n.get_border(), '\n'

    # get string with WKT representation of the border polygon
    print 'BorderPolygon:', n.get_border_polygon(), '\n'

    # get list with coordinates of the object corners
    print 'Corners:', n.get_corners(), '\n'

    # longitude and latitude grids representing the full data grid
    longitude, latitude = n.get_geolocation_grids()
    print 'longitude :', longitude, '\n'
    print 'latitude :', latitude, '\n'

    # Create Domain object. It describes the desired grid of reprojected image:
    #    projection, resolution, size, etc. In this case it is geographic projection;
    #    -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
    dLatlong = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", "-te 25 70 35 72 -ts 2000 2000")
    dLatlong.write_map(oFileName + '01_latlong_map.png')
    print 'Latlong Domain:', dLatlong, '\n'

    # Generate two vectors with values of lat/lon for the border of domain
    lonVec, latVec = n.get_border()
    print 'lonVec :', lonVec, '\n'
    print 'latVec :', latVec, '\n'

    # Get upwards azimuth direction of domain.
    bearing_center = n.upwards_azimuth_direction()
    print 'bearing_center :', bearing_center, '\n'

    # Create domain with stereographic projection
    # -- Get corners of the image
    lons, lats = n.get_corners()
    meanLon = sum(lons, 0.0) / 4.
    meanLat = sum(lats, 0.0) / 4.
    srsString = "+proj=stere +lon_0=%f +lat_0=%f +k=1 +ellps=WGS84 +datum=WGS84 +no_defs" % (meanLon, meanLat)
    extentString = '-lle %f %f %f %f -tr 100 100' % (min(lons), min(lats), max(lons), max(lats))
    # -- Create Domain with stereographic projection, corner coordinates and resolution 1000m
    dStereo = Domain(srsString, extentString)
    dStereo.write_map(oFileName + '02_stereo_map.png')
    print 'Stereo Domain:', dStereo, '\n'

    # Reproject all GCPs to a new spatial reference system
    n.reproject_GCPs(srsString)
    for iGCP in range(10):
        print n.vrt.dataset.GetGCPs()[iGCP]
    print ''

    # make KML file with image borders (to be opened in Googe Earth)
    n.write_kml(kmlFileName=oFileName + '03_preview.kml')

    print '\n*** domain_test completed successfully. Output files are found here:' + oFileName

main()
