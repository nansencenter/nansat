#!/usr/bin/env python
# Name:    tutorial.py
# Purpose: Tutorial with list of examples
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

from nansat import Nansat, Domain

# input and output file names
iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
iFileName = os.path.join(iPath, 'gcps.tif')
print 'Input file: ', iFileName
oFileName = os.path.join(iPath, 'tmpdata', 'output_')
print 'Output file prefix: ', oFileName

# Open an input file
# Create a Nansat object <n> for futher high-level operations
n = Nansat(iFileName)

# Open an input file, specify which Mapper to use, set logging level
n = Nansat(iFileName, mapperName='generic', logLevel=10)

# list bands and georeference of the object
print 'Raw Nansat:', n

# get dictionary with metadata from all bands 
print 'Bands:', n.bands()

# get size of the object (Y and X dimensions, to follow Numpy style)
print 'Shape:', n.shape()

# get list with coordinates of the object corners
print 'Corners:', n.get_corners()

# get lists with coordinates of the object borders
print 'Border:', n.get_border()

# get time of the image aquisition
print 'Time:', n.get_time()[0]

# Get band data and do some operations
# 1. Get data from 1st band as numpy array
# 2. Plot the array (pyplot image is save to a PNG file)
# 3. Save as Matlab file
a = n[1]
plt.imshow(a);plt.colorbar();plt.savefig(oFileName + '_imshow.png');plt.close()
savemat(oFileName + '.mat', {'band_1': a})

# make simple indexed image from 1st band with default colormap
n.write_figure(oFileName + '.png')

# make RGB image from bands 1,2,3 with brightness correction
n.write_figure(oFileName + '_rgb.png', bands=[1,2,3], clim='hist', ratio=0.9)

# make indexed image with legend
n.write_figure(oFileName + '_legend.png', legend=True, titleString='NANSAT Tutorial', LEGEND_HEIGHT=0.3, fontSize=10)

# get array with watermask (landmask)
# 1. Get Nansat object with watermask
# 2. Get array from Nansat object. 0 - land, 1 - water
wm = n.watermask()
wmArray = wm[1]

# write figure with land overlay (gray color) and apply brightness gamma correction
n.write_figure(oFileName + '_land.png', clim='hist', mask_array=wmArray, mask_lut={0: [128, 128, 128]}, logarithm=True, gamma=3)

# add logo to image to the lower left corner
# (make sure file is in the current folder)
n.write_figure(oFileName + '_logo.png', logoFileName='nansat_logo_s.png', logoLocation=[10, -10], logoSize=[70, 70])

# write figure with lat/lon grid
# 1. Get lat/lon arrays from Nansat object (may take some time)
# 2. Make figure with lat/lon grids
lonGrid, latGrid = n.get_geolocation_grids()
n.write_figure(oFileName + '_latlon.png', latGrid=latGrid, lonGrid=lonGrid, latlonGridSpacing=10, latlonLabels=10)

# make small preview
# 1. Reduce size of the Nansat object ten times
# 2. Make simple grayscaled image with brightness correction
# 3. Resize back to original resolution
n.resize(0.1)
n.write_figure(oFileName + '_preview.png', clim='hist', cmapName='gray')
n.resize()

# enlarge the image with bicubic-spline interpolation
n.resize(2, eResampleAlg=3)
n.write_figure(oFileName + '_large.png', clim='hist', cmapName='gray')
n.resize()

# make KML file with image borders (to be opened in Googe Earth)
n.write_kml(kmlFileName=oFileName + '_preview.kml')

# make image with map of the file location
n.write_map(oFileName + '_map.png')


# Make image reprojected onto map of Northern Europe
# 1. Create Domain object. It describes the desired grid of reprojected image:
# projection, resolution, size, etc. In this case it is geographic projection;
# -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
# 2. Reproject the Nansat object
# 3. Make simple image
dLatlong = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", "-te 25 70 35 72 -ts 2000 2000")
dLatlong.write_map(oFileName + '_latlong_map.png')
print 'Latlong Domain:', dLatlong
n.reproject(dLatlong)
n.write_figure(oFileName + '_pro_latlon.png')

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
dStereo.write_map(oFileName + '_stereo_map.png')
print 'Stereo Domain:', dStereo
n.reproject(dStereo, 2) # 4.
n.write_figure(oFileName + '_pro_stereo.png') # 5.

# Reproject onto grid of another file (destination file is projected)
n2 = Nansat('stere.tif')
n.reproject(n2)
n.write_figure(fileName=oFileName + '_proj_1onto2.png', bands=[1,2,3], clim='hist')

# Reproject onto grid of another file (destination file is in swath projection)
n.reproject()
n2.reproject(n)
n2.write_figure(fileName=oFileName + '_proj_2onto1.png', bands=[1,2,3], clim='hist')

# Reproject onto grids of lat/lon
dFromGrids = Domain(lon=lonGrid, lat=latGrid)
n2.reproject(dFromGrids)
n2.write_figure(fileName=oFileName + '_proj_on_grid.png', bands=[1,2,3], clim='hist')

# reproject onto automatically generated domain
dstDomainAuto = Domain(srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", ds=n.raw.dataset)
n.reproject(dstDomainAuto)
n.write_figure(fileName=oFileName + '_proj_1auto.png', bands=[1,2,3], clim='hist')

# export all data into NetCDF format
n.export(oFileName + '_0.nc')

# export one band to GeoTIFF
n.export_band(oFileName + '_2.tif', bandID=2, driver='GTiff')

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
n.write_figure(fileName=oFileName + '_proj.png', bands=[1,2,3],
               mask_array=wmArray, mask_lut={0: [128, 128, 128]},
               clim='hist', transparency=[128, 128, 128])           # 8.

# make KML file for the exported image
n.write_kml_image(kmlFileName=oFileName + '.kml', kmlFigureName=oFileName + '_proj.png')

# Perform batch averaging of several files
# 1. Create destination Nansat object with desired projection
nMosaic = Nansat(domain=dStereo)
# 2. Perfom averaging
nMosaic.mosaic(['gcps.tif', 'stere.tif'], bands=['L_645', 'L_555', 'L_469'])
# 3. Get mask of non-valid pixels
mask = nMosaic['mask']
# 4. Output averaged data using the mask
nMosaic.write_figure(fileName=oFileName + '_mosaic.png', bands=['L_645', 'L_555', 'L_469'], clim='hist',
                        mask_array=mask, mask_lut={0:[128,128,128]})

print 'Tutorial completed successfully. Output files are found here:' + oFileName
