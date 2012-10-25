#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat


from nansat import Nansat
from domain import Domain

iPath = '/data/gdal_test/'
oPath = '/data/'
fileName = 'MER_FRS_1CNPDK20110503_105735_000000733102_00109_47968_7898.N1'
oFileName = oPath + fileName

# create Nansat object
n = Nansat(iPath + fileName)

# list bands and georeference of the object
print n

# get dictionary with all bands metadata
print n.bands()

# get size of the object (Y and X dimensions, to follow Numpy style)
print n.shape()

# get list with coordinates of the object corners
print n.get_corners()

# get lists with coordinates of the object borders
print n.get_border()

# Get band data and do some operations
# 1. Get data from 1st band as numpy array
# 2. Plot the array (NB: close the shown image to continue processing)
# 3. Save as Matlab file
a = n[1]
plt.imshow(a);plt.colorbar();plt.show()
savemat(oFileName + '.mat', {'band_1': a})

# make simple indexed image from 1st band with default colormap
n.write_figure(oFileName + '.png')

# make RGB image from bands 6,5,2 with brightness correction
n.write_figure(oFileName + '_rgb.png', bands=[6,5,2], clim='hist', ratio=0.7, logarithm=True, gamma=3)

# make indexed image with legend
n.write_figure(oFileName + '_legend.png', bands='radiance_900', legend=True, titleString='NANSAT Tutorial', fontSize=40)

# get array with watermask (landmask)
# 1. Get Nansat object with watermask
# 2. Get array from Nansat object. 0 - land, 1 - water
wm = n.watermask()
wmArray = wm[1]

# write figure with land overlay (gray color)
n.write_figure(oFileName + '_land.png', mask_array=wmArray, mask_lut={0: [0.5, 0.5, 0.5]})

# add logo to image to the lower left corner
# (make sure file is in the current folder)
n.write_figure(oFileName + '_logo.png', logoFileName='nansat_logo_s.png', logoLocation=[10, -10], logoSize=[150, 100])

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

# make KML file with image borders (to be opened in Googe Earth)
n.write_kml(kmlFileName=oPath + fileName + '_preview.kml')

# make image with map of the file location
n.write_map(oFileName + '_map.png')

# Make image reprojected onto map of Northern Europe
# 1. Create Domain object. It describes the desired grid of reprojected image:
# projection, resolution, size, etc. In this case it is geographic projection;
# -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
# 2. Reproject the Nansat object
# 3. Make simple image
dLatlong = Domain(srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", ext="-lle -10 50 30 70 -ts 2000 2000")
dLatlong.write_map(oFileName + '_latlong_map.png')
print dLatlong
n.reproject(dLatlong)
n.write_figure(oFileName + '_pro_latlon.png')

# Reprojected image into stereographic projection
# 1. Cancel previous reprojection
# 2. Get corners of the image
# 3. Create Domain with stereographic projection, corner coordinates and resolution 1000m
# 4. Reproject
# 5. Write image
n.reproject() # 1.
lons, lats = n.get_corners() # 2.
meanLon = sum(lons, 0.0) / 4.
meanLat = sum(lats, 0.0) / 4.
srsString = "+proj=stere +lon_0=%f +lat_0=%f +k=1 +ellps=WGS84 +datum=WGS84 +no_defs" % (meanLon, meanLat)
extentString = '-lle %f %f %f %f -tr 1000 1000' % (min(lons), min(lats), max(lons), max(lats))
dStereo = Domain(srs=srsString, ext=extentString) # 3.
dStereo.write_map(oFileName + '_stereo_map.png')
print dStereo
n.reproject(dStereo) # 4.
n.write_figure(oFileName + '_pro_stereo.png') # 5.

# reproject object onto given lat/lon arrays
n.reproject(Domain(lon=lonGrid, lat=latGrid))

# export all data into NetCDF format (may take some time and lot's of space)
n.export(oFileName + '.nc')

# export only few bands into NetCDF format
n.export(oFileName + '.nc', bands=[1,2,3])

# export only few bands and remove some metadata
n.export(oFileName + '.nc', bands=[1,2,3], rmMetadata=['SourceFilename'])

# export one band to GeoTIFF
n.export(oFileName + '.nc', bands=[4], driver='GeoTIFF')




# add band from array to existing object
# 0. Cancel reprojection. Adding bands works only on non-reprojected data
# 1. Get the data from the object and modify
# 2. Add band with modified data to the object
# 3. Check that the band is in the ibject
n.reproject()
arrayMeta = n.get_metadata(bandID=1)
array = n[1] * 10
n.add_band(array=array, parameters=arrayMeta) #{'BandName': 'new_band_1', 'description': 'Modified values of first band'})
print n

# add band from another file to existing object
# (actually the same file in this example but it may be any other file)
n.add_band(fileName=iPath+fileName, bandID=1, parameters={'BandName': 'one_more_band_1'})

# create new object from given domain and array
# 1. Reproject the current object and get data
# 2. Create new Nansat object
n.reproject(dStereo)
array = n[1]
nNew = Nansat(domain=dStereo, array=array, parameters={'BandName': 'first_band'})
