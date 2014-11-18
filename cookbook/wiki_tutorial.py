#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:         wiki_tutorial.py
# Purpose:      Code for the tutorial on the github wiki
#
# Author:       Anton Korosov, Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	11.11.2014
# Last modified:11.11.2014 10:27
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import sys, os
home = os.path.expanduser("~")

import numpy as np
import matplotlib.pyplot as plt

from nansat.nansatmap import Nansatmap
from nansat.nansat import Nansat, Domain

iFileName = os.path.join(home,
                'python/nansat/nansat/tests/data/gcps.tif')

# Open an input satellite image with Nansat
n = Nansat(iFileName)

# List bands and georeference of the object
print n

# Write picture with map of the file location
n.write_map('map.png')

# Write indexed picture with data from the first band
n.write_figure('rgb.png', clim='hist')

# Reproject input image onto map of Norwegian Coast
# 1. Create domain describing the desired map
# 2. Transform the original satellite image
# 3. Write the transfromed image into RGB picture
dLatlong = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                          "-te 27 70.2 31 71.5 -ts 2000 2000")
n.reproject(dLatlong)
n.write_figure('pro.png', bands=[1,2,3], clim=[0, 100])

# Export projected satelite image into NetCDF format
n.export('gcps_projected.nc')

# Collect values from interactively drawn transect
# 1. draw transect interactively
# 2. plot the values
values, lonlat, pixlinCoord =n.get_transect()
plt.plot(lonlat['shape0']['longitude'], values['1:L_645']['shape0'], '.-');plt.show()


