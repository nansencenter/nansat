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

from nansat import Nansat, Domain#, Mosaic
from mosaic import Mosaic

''' Domain class is a container for geographical reference of a raster

    A Domain object describes all attributes of geographical
    reference of a raster:
        width and height, pixel size,
        relation between pixel/line coordinates and geographical coordinates,
        type of data projection (e.g. geographical or stereographic)

    The core of Domain is a GDAL Dataset with no bands.
    It has only georeference information: rasterXsize, rasterYsize,
    GeoTransform and Projection or GCPs, etc. which fully describe dimentions
    and spatial reference of the grid.

    There are three ways to store geo-reference in a GDAL dataset:
      * Using GeoTransfrom
      * Using GCPs
      * Using Geolocation Array
    The relation between X/Y coordinates of the raster and latitude/longitude
    coordinates is defined by projection type and projection parameters.
    These pieces of information are therefore stored in Domain:
      * Type and parameters of projection +
      * GeoTransform or GCPs or GeolocationArrays

    Domain has methods for basic operations with georeference information:
      * creating georeference from input options;
      * fetching corner, border or full grids of X/Y coordinates;
      * making map of the georeferenced grid in a PNG or KML file;
      * and some more...

'''
# Get input and output file names
from testio import testio
iPath, iFileName, oPath, oFileName, shpFileName = testio()
oFileName = oFileName+'domain_'

# Open an input file
# Create a Nansat object <n> for futher high-level operations
n = Nansat(iFileName)

# Open an input file, specify which Mapper to use, set logging level
n = Nansat(iFileName, mapperName='generic', logLevel=10)

# Make image reprojected onto map of Northern Europe
# 1. Create Domain object. It describes the desired grid of reprojected image:
# projection, resolution, size, etc. In this case it is geographic projection;
# -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
# 2. Reproject the Nansat object
# 3. Make simple image
dLatlong = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", "-te 25 70 35 72 -ts 2000 2000")
dLatlong.write_map(oFileName + '01_latlong_map.png')
print 'Latlong Domain:', dLatlong
n.reproject(dLatlong)
n.write_figure(oFileName + '02_pro_latlon.png')

# Reprojected image into stereographic projection
# 1. Cancel previous reprojection
# 2. Get corners of the image
# 3. Create Domain with stereographic projection, corner coordinates and resolution 1000m
# 4. Reproject with cubic interpolation
# 5. Write image
n.reproject() # 1.
lons, lats = n.get_corners() # 2.
meanLon = sum(lons, 0.0) / 4.
meanLat = sum(lats, 0.0) / 4.
srsString = "+proj=stere +lon_0=%f +lat_0=%f +k=1 +ellps=WGS84 +datum=WGS84 +no_defs" % (meanLon, meanLat)
extentString = '-lle %f %f %f %f -tr 100 100' % (min(lons), min(lats), max(lons), max(lats))
dStereo = Domain(srsString, extentString) # 3.
dStereo.write_map(oFileName + '03_stereo_map.png')
print 'Stereo Domain:', dStereo
n.reproject(dStereo, 2) # 4.
n.write_figure(oFileName + '04_pro_stereo.png') # 5.

# Reproject onto grid of another file (destination file is projected)
n2 = Nansat('stere.tif')
n.reproject(n2)
n.write_figure(fileName=oFileName + '05_proj_1onto2.png', bands=[1,2,3], clim='hist')

# Reproject onto grid of another file (destination file is in swath projection)
n.reproject()
n2.reproject(n)
n2.write_figure(fileName=oFileName + '06_proj_2onto1.png', bands=[1,2,3], clim='hist')

# Reproject onto grids of lat/lon
lonGrid, latGrid = n.get_geolocation_grids()
dFromGrids = Domain(lon=lonGrid, lat=latGrid)
n2.reproject(dFromGrids)
n2.write_figure(fileName=oFileName + '07_proj_on_grid.png', bands=[1,2,3], clim='hist')

# reproject onto automatically generated domain
dstDomainAuto = Domain(srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", ds=n.raw.dataset)
n.reproject(dstDomainAuto)
n.write_figure(fileName=oFileName + '08_proj_1auto.png', bands=[1,2,3], clim='hist')

# export all data into NetCDF format
n.export(oFileName + '08.nc')

# export one band to GeoTIFF
n.export_band(oFileName + '08.tif', bandID=2, driver='GTiff')

# create new object from given domain and array
# 1. Reproject the current object
# 2. Get array with data
# 2. Create new Nansat object from the given array and for given domain
n.reproject(dStereo)
array = n[1]
nStereo = Nansat(domain=dStereo, array=array, parameters={'name': 'band1'})
print 'Stereo Nansat:', nStereo

# add band from array to existing object
# 0. Cancel reprojection. Adding bands works only on non-reprojected data
# 1. Get the data from the object and modify
# 2. Add band with modified data to the object
# 3. Check that the band is in the object
n.reproject()
array = n[1] * 10
n.add_band(array=array, parameters={'name': 'new_band', 'about': 'test'})
print 'Nansat with new band:', n

# add band from another file to existing object
n.add_band(fileName='stere.tif', bandID='L_645', parameters={'name': 'L_645_stere'})
print 'Nansat with new band from another file:', n

# Prepare the image for Google Earth exporting
# reproject image into Lat/Lon WGS84 (Simple Cylindrical) projection
# 1. Cancel previous reprojection
# 2. Get corners of the image and the pixel resolution
# 3. Create Domain with stereographic projection, corner coordinates 1000m
# 4. Reproject
n.reproject()                                                       # 1.
lons, lats = n.get_corners()                                        # 2.
srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs"
extentString = '-lle %f %f %f %f -ts 800 600' % (min(lons), min(lats), max(lons), max(lats))
d = Domain(srs=srsString, ext=extentString)                         # 3.
n.reproject(d)                                                      # 4.
# get array with watermask (landmask)
# it must be done after reprojection!
# 6. Get Nansat object with watermask
# 7. Get array from Nansat object. 0 - land, 1 - water
wm = n.watermask()                                                  # 6.
wmArray = wm[1]                                                     # 7.

# 8. Write the projected image with transparent land mask and image background
# transparentMask: boolean, defult = False
# If True, the masked pixels will be transparent when saving to png
#transparency: int
#transparency of the image background, set for PIL in Figure.save()
#default transparent color is [0,0,0]
n.write_figure(fileName=oFileName + '09_proj.png', bands=[1,2,3],
               mask_array=wmArray, mask_lut={0: [128, 128, 128]},
               clim='hist', transparency=[128, 128, 128])           # 8.

# make KML file for the exported image
n.write_kml_image(kmlFileName=oFileName + '10.kml', kmlFigureName=oFileName + '_proj.png')

# Perform batch averaging of several files
# 1. Create destination Nansat object with desired projection
nMosaic = Mosaic(domain=dStereo)
# 2. Perfom averaging
nMosaic.average(['gcps.tif', 'stere.tif'], bands=['L_645', 'L_555', 'L_469'])
# 3. Get mask of non-valid pixels
mask = nMosaic['mask']
# 4. Output averaged data using the mask
nMosaic.write_figure(fileName=oFileName + '11_mosaic.png', bands=['L_645', 'L_555', 'L_469'], clim='hist',
                        mask_array=mask, mask_lut={0:[128,128,128]})

print 'domain_test completed successfully. Output files are found here:' + oFileName

