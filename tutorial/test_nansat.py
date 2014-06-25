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
import os

from nansat import *

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

# create nansat object from grids of lon, lat and values
# first get the grids from existing nansat objects
lonGrid, latGrid = n.get_geolocation_grids()
valGrid = n[1]
# next create a domain
d = Domain(lon=lonGrid, lat=latGrid)
# at last generate nansat from given domain and array
n2 = Nansat(domain=d, array=valGrid, parameters={'name': 'new_band'})
print n


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

# add a band from numpy array (copy the 1st band to the end (5th band))
n.add_band(n[1], parameters={'name': 'Name1',
                             'info':  'copy from the 1st band array'})
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

# Resize the data to 50% using averaging
n.resize(0.5)
# make simple indexed image from 1st band with default colormap
n.write_figure(oFileName + '02.png', clim='hist')
# undo resize
n.undo()

# Smooth the data by resizing to 50% and to 200% using CubicSpline
n.resize(0.5, eResampleAlg=3)
n.resize(2.0, eResampleAlg=3)
# make simple indexed image from 1st band with default colormap
n.write_figure(oFileName + '02CubicSpline.png', clim='hist')
# undo smoothing (undo two resizing)
n.undo(2)

# make image with map of the file location
n.write_map(oFileName + '04_map.png')

# Writes an 8-bit GeoTiff image for a given band
n.write_geotiffimage(oFileName + '05_geotiff.tif', bandID=1)

# create a NetCDF file with all bands
n.export(oFileName + '06a.nc')
n.export(oFileName + '06b.nc', bottomup=True)

# create a GTiff file with one band (default driver is NetCDF)
n.export_band(oFileName + '07.tif', bandID=1, driver='GTiff')

n.crop(lonlim=[28, 29], latlim=[70.7, 71])
# get array with watermask (landmask)
# -- Get Nansat object with watermask
wm = n.watermask()[1]

# -- Reproject with cubic interpolation
d = Domain(4326, "-te 27 70.3 31 71.5 -ts 300 300")
n.reproject(d, 2)
# -- Write image
n.write_figure(oFileName + '08_pro.png', clim='hist')
print n.vrt
n.undo()
print n.vrt

# crop image using xOff, yOff, xSize, ySize
n.crop(100, 110, 50, 60)
n.write_figure(oFileName + '09_crop.png', clim='hist')
n.undo()

# crop image using lonlim, latlim
n.crop(lonlim=[28, 29], latlim=[70.7, 71])

# automatically reproject cropped image
d = Domain(4326, ds=n.vrt.dataset)
n.reproject(d)
n.write_figure(oFileName + '09_crop_pro.png', clim='hist')
# undo croping and reproject
n.undo(100)

# Get transect of the 1st and 2nd bands corresponding to the given points
points=((29.287, 71.153),
        (29.275, 71.145),
        (29.210, 71.154))
#import pdb; pdb.set_trace()
values, lonlat, pixlinCoord = n.get_transect(points,
                                             transect=False,
                                             bandList=[1, 2])
# print the results
print '1stBandVal  2ndBandVal       pix/lin         lon/lat '
for i in range (len(values[0])):
    print '%6d %10d %13.2f /%6.2f  %7.2f /%6.2f' % (values[0][i],
                                                    values[1][i],
                                                    pixlinCoord[0][i],
                                                    pixlinCoord[1][i],
                                                    lonlat[0][i],
                                                    lonlat[1][i])
print ''

ogrObject = n.get_transect(points, returnOGR=True)
ogrObject.export(oFileName + '_10_transect.shp')

# export into THREDDS friendly file
iFileName = os.path.join(iPath, 'stere.tif')
n = Nansat(iFileName)
n.export2thredds(oFileName + 'thredds.nc', [1,2,3])


print '\n***nansat_test completed successfully. Output files are found here:' + oFileName

