#!/usr/bin/env python
# Name:         tutorial.py
# Author:       Anton Korosov
# Licence:      GPL v.3
# This file is part of NANSAT. You can redistribute it or modify
# under the terms of GNU General Public License, v.3
# http://www.gnu.org/licenses/gpl-3.0.html

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import inspect, os

from nansat import Nansat, Domain, Mosaic
from testio import testio

# Set input and output file names
iPath, oPath = testio()
iFileName = os.path.join(iPath, 'gcps.tif')
oFileName = os.path.join(oPath, 'tutor_')

# Open an input satellite image with Nansat 
n = Nansat(iFileName)

# List bands and georeference of the object
print n

# Write picture with map of the file location
n.write_map(oFileName + 'map.png')

# Write indexed picture with data from the first band
n.write_figure(oFileName + '.png', clim='hist')

# Reproject input image onto map of Norwegian Coast
# 1. Create domain describing the desired map
# 2. Transform the original satellite image
# 3. Write the transfromed image into RGB picture
dLatlong = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                  "-te 27 70.2 31 71.5 -ts 500 500")
n.reproject(dLatlong)
n.write_figure(oFileName + 'pro.png', bands=[1,2,3], clim=[0, 100])

# Export projected satelite image into NetCDF format
n.export(oFileName + '.nc')

# Collect values from interactively drawn transect
# 1. draw transect interactively
# 2. plot the values
values, lonlat, pixlinCoord =n.get_transect()
plt.plot(lonlat[0], values[0], '.-');plt.show()

# run tests of other nansat components
import test_domain
import test_nansat
import test_figure
import test_nansatmap
import test_nansatshape
import test_mosaic
import test_pointbrowser
