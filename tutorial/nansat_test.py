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

import matplotlib.pyplot as plt
from scipy.io import savemat
import inspect, os

from nansat import Nansat

# input and output file names
iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
iFileName = os.path.join(iPath, 'gcps.tif')
print 'Input file: ', iFileName
oPath = os.path.join(iPath, 'tmpdata')
print 'Output path:', oPath
if not os.path.exists(oPath):
    os.mkdir(oPath)
oFileName = os.path.join(oPath, 'output_')
print 'Output file:', oFileName
print ''

# Open an input file
# Create a Nansat object <n> for futher high-level operations
n = Nansat(iFileName)

# Open an input file, specify which Mapper to use, set logging level
n = Nansat(iFileName, mapperName='generic', logLevel=10)

# list bands and georeference of the object
print 'Raw Nansat:', n
print ''

# get dictionary with metadata from all bands
print 'Bands:', n.bands()
print ''

# get size of the object (Y and X dimensions, to follow Numpy style)
print 'Shape:', n.shape()
print ''

# get list with coordinates of the object corners
print 'Corners:', n.get_corners()
print ''

# get lists with coordinates of the object borders
print 'Border:', n.get_border()
print ''

# get time of the image aquisition
print 'Time:', n.get_time()[0]
print ''

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
n.write_figure(oFileName + '_land.png', clim='hist', mask_array=wmArray, mask_lut={2: [128, 128, 128]}, logarithm=True, gamma=3)

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


n.write_nansatmap(oFileName + '_nansatMap.png', )

print 'nansat_test completed successfully. Output files are found here:' + oFileName
