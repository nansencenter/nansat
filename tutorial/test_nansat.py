#!/usr/bin/env python
# Name:    nansat_test.py
# Purpose: Tutorial for nansat
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

from nansat import Nansat, Domain
from mosaic import Mosaic

# input and output file names
from testio import testio
iPath, oPath = testio()
iFileName = os.path.join(iPath, 'gcps.tif')
print 'Input file: ', iFileName
oFileName = os.path.join(oPath, 'output_nansat_')
print 'Output file:', oFileName

'''Nansat is a container for geospatial data, performs all high-level operations

n = Nansat(fileName) opens the file with satellite or model data for
reading, adds scientific metadata to bands, and prepares the data for
further handling.

The instance of Nansat class (the object <n>) contains information
about geographical reference of the data (e.g raster size, pixel
resolution, type of projection, etc) and about bands with values of
geophysical variables (e.g. water leaving radiance, normalized radar
cross section, chlrophyll concentraion, etc). The object <n> has methods
for high-level operations with data. E.g.:
* reading data from file (Nansat.__getitem__);
* visualization (Nansat.write_figure);
* changing geographical reference (Nansat.reproject);
* exporting (Nansat.export)
* and much more...

Nansat inherits from Domain (container of geo-reference information)

'''

# Open an input file
# Create a Nansat object <n> for futher high-level operations
n = Nansat(iFileName)

# Open an input file, specify which Mapper to use, set logging level
n = Nansat(iFileName, mapperName='generic', logLevel=10)

# list bands and georeference of the object
print 'Raw Nansat:', n, '\n'

# get dictionary with metadata from all bands
print 'Bands:', n.bands(), '\n'

# get time of the image aquisition
print 'Time:', n.get_time()[0], '\n'

# set GlobalMetadata
n.set_metadata(key='GlobalKey', value='GlobalVal')
# get Global Metadata
print 'Global Metadata:', n.get_metadata(), '\n'

# set BandMetadata to the 1st band
n.set_metadata(key='BandKey', value='BandVal', bandID=1)
# get 1st Band Metadata
print '1st Band Metadata:', n.get_metadata(bandID=1), '\n'

# add a band from file (copy the 2nd band to the end (4th band)
n.add_band(fileName=n.fileName, bandID=2)
# add a band from numpy array (copy the 1st band to the end (5th band))
n.add_band(array=n[1], parameters={'name':'Name1', 'wkt':'copy from the 1st band array'})
# print band list
n.list_bands()
# get GDAL raster band (2nd band)
band = n.get_GDALRasterBand(bandID=2)

# Get band data and do some operations
# -- Get data from 1st band as numpy array
a = n[1]
# -- Plot the array (pyplot image is save to a PNG file)
plt.imshow(a);plt.colorbar();plt.savefig(oFileName + '01_imshow.png');plt.close()
# -- Save as Matlab file
savemat(oFileName + '01.mat', {'band_1': a})

# Resize the data to 50%
n.resize(0.5)
# make simple indexed image from 1st band with default colormap
n.write_figure(oFileName + '02.png', clim='hist')
# undo resize
n.resize()

# Resize the data to 50% using CubicSpline
n.resize_lite(0.5, eResampleAlg=3)
# make simple indexed image from 1st band with default colormap
n.write_figure(oFileName + '02CubicSpline.png', clim='hist')
# undo resize
n.resize()


# make image with map of the file location
n.write_map(oFileName + '04_map.png')

# Writes an 8-bit GeoTiff image for a given band
n.write_geotiffimage(oFileName + '05_geotiff.tif', bandID=1)

# create a NetCDF file with all bands
n.export(oFileName + '06a.nc')
n.export(oFileName + '06b.nc', bottomup=True)

# create a GTiff file with one band (default driver is NetCDF)
n.export_band(oFileName + '07.tif', bandID=1, driver='GTiff')

# get array with watermask (landmask)
# -- Get Nansat object with watermask
wm = n.watermask()[1]

# -- Reproject with cubic interpolation
d = Domain(4326, "-te 27 70.3 31 71.5 -ts 300 300")
n.reproject(d, 2)
# -- Write image
n.write_figure(oFileName + '08_pro.png', clim='hist')

# Get transect of the 1st and 2nd bands corresponding to the given points
values, lonlat, pixlinCoord =n.get_transect(points=((29.287, 71.153), (29.275, 71.145), (29.210, 71.154)), transect=False, bandList=[1, 2])
# print the results
print '1stBandVal  2ndBandVal       pix/lin         lon/lat '
for i in range (len(values[0])):
    print '%6d %10d %13.2f /%6.2f  %7.2f /%6.2f' %(values[0][i], values[1][i], pixlinCoord[0][i], pixlinCoord[1][i], lonlat[0][i], lonlat[1][i])
print ''

print '\n***nansat_test completed successfully. Output files are found here:' + oFileName

