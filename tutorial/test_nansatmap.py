#!/usr/bin/env python
# Name:    nansatmap_test.py
# Purpose: Tutorial for nansatmap class
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

from nansat import Nansat, Domain, Nansatmap

# input and output file names
from testio import testio
iPath, oPath = testio()
iFileName = os.path.join(iPath, 'map.tif')
print 'Input file: ', iFileName
oFileName = os.path.join(oPath, 'output_nansatmap_')
print 'Output file:', oFileName

''' Nansatmap class perform opeartions with graphical files:
create, add legend and geolocation_grids, save.

The core of Nansatmap is Basemap.
Nansatmap object is created based on Domain object.
This class has methods as: filled contour plot, line contour plot,
pseudo-color plot, quiver plot and save.

'''

# Create a Nansat object (n)
n = Nansat(iFileName)
# Get data from 1st, 2nd and 3rd bands as numpy array (u,v and w)
u = n[1]; v = n[2]; w = n[3]

# Create Nansatmap object from nansat (domain) object
nMap = Nansatmap(n)
# draw filled contour plot
nMap.contourf(w, v=range(4,22,2))
# draw black smooth contour plot with labels
nMap.contour(w, smooth=True, fontsize=8, colors='k')
# add colorbar
nMap.add_colorbar(fontsize=10)
# add geocoordinates
nMap.drawgrid()
# save to a file
nMap.save(oFileName+'01_contourf_contour.png', landmask=False)

# Create Nansatmap object from nansat (domain) object
nMap = Nansatmap(n, resolution='l')
# pseudo-color plot over the map
nMap.pcolormesh(w)
# quiver plot
nMap.quiver(u, v, quivectors=20)
# save to a file
nMap.save(oFileName+'02_pcolormesh_quiver.png')

# use Nansatmap for converting lon/lat into x/y
# 1. Create domain over the area of interest in stereographic projection
extentString = '-lle -10 50 20 70 -tr 1000 1000'
srsString = '+proj=stere +lon_0=10 +lat_0=60 +k=1 +ellps=WGS84 +datum=WGS84 +no_defs'
d = Domain(srsString, extentString)
# 2. Create nansatmap object from tha domain
nmap = Nansatmap(d)
# 3. Use the created object to convert from lon/lat into x/y:
x, y = nmap([0, 2, 4], [63, 64, 65])
# 4 or from x/y into lon/lat
lon, lat = nmap(x, y, inverse=True)

print '\n*** nansatmap_test completed successfully. Output files are found here:' + oFileName

