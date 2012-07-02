#!/usr/bin/env python

from nansat import Nansat
from domain import Domain
import matplotlib.pyplot as plt
from scipy.io import savemat

iPath = '/data/gdal_test/'
oPath = '/data/'
fileName = 'MER_FRS_1CNPDK20110503_105735_000000733102_00109_47968_7898.N1'

# create Nansat object
n = Nansat(iPath + fileName)

# list bands and georeference of the object
print n

# Get band data and do some operations
# 1. Get data from 1st band as numpy array
# 2. Plot the array (NB: close the shown image to continue processing)
# 3. Save as Matlab file
a = n[1]
plt.imshow(a);plt.colorbar();plt.show()
savemat(oPath + fileName + '.mat', {'band_1': a})

# make simple indexed image from 1st band
n.write_figure(oPath + fileName + '.png')

# make RGB image from bands 6,5,2 with brightness correction
n.write_figure(oPath + fileName + '_rgb.png', bands=[6,5,2], clim='hist', ratio=0.7, logarithm=True, gamma=3)

# make indexed image with legend
n.write_figure(oPath + fileName + '_legend.png', bands='radiance_900', clim='hist', ratio=0.7, legend=True, titleString='NANSAT Tutorial', fontSize=40)

# make small preview in three steps:
# 1. Reduce size of the Nansat object ten times
# 2. Make simple grayscaled image with brightness correction
# 3. Resize back to original resolution 
n.resize(0.1)
n.write_figure(oPath + fileName + '_preview.png', clim='hist', ratio=0.7, cmapName='gray')
n.resize()

# make KML file with image borders (to be opened in Googe Earth)
n.write_kml(kmlFileName=oPath + fileName + '_preview.kml')

# make image with map of the file location
n.write_map(oPath + fileName + '_map.png')

# Make image, reprojected onto map of Northern Europe in three steps:
# 1. Create Domain object. It describes the desired grid of reprojected image:
# projection, resolution, size, etc. In this case it is geographic projection;
# -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
# 2. Reproject the Nansat object
# 3. Make simple image
d = Domain(srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", ext="-lle -10 50 30 70 -ts 2000 2000")
print d
n.reproject(d)
n.write_figure(oPath + fileName + '_pro_latlon.png')


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
d = Domain(srs=srsString, ext=extentString) # 3.
print d
n.reproject(d) # 4.
n.write_figure(oPath + fileName + '_pro_stereo.png') # 5.
