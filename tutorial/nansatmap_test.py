#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import inspect, os

from nansat import Nansat, Domain, Nansatmap
#from nansat_map import NansatMap

# input and output file names
iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
iFileName = os.path.join(iPath, 'map.tif')
print 'Input file: ', iFileName
oFileName = os.path.join(iPath, 'tmpdata', 'outputmap_')
print 'Output file prefix: ', oFileName

# Create a Nansat object (n)
n = Nansat(iFileName)
# Get data from 1st, 2nd and 3rd bands as numpy array (u,v and w)
u = n[1]; v = n[2]; w = n[3]

nMap = Nansatmap(n)
# draw filled contour plot
nMap.contourf(w)
# draw black smooth contour plot with labels
nMap.contour(w, smooth=True, fontsize=8, colors='k')
# add colorbar
nMap.add_colorbar(fontsize=10)
# add geocoordinates
nMap.drawgrid()
# save to a file
nMap.save(oFileName+'contourf_contour.png', landmask=False)


nMap = Nansatmap(n, resolution='h')
# pseudo-color plot over the map
nMap.pcolormesh(w)
# quiver plot
nMap.quiver(u, v, quivectors=20)
nMap.save(oFileName+'pcolormesh_quiver.png')

# use Nansatmap for converting lon/lat into x/y
# 1. Create domain over the area of interest in stereographic projection
extentString = '-lle -10 50 20 70 -tr 1000 1000'
srsString = '+proj=stere +lon_0=10 +lat_0=60 +k=1 +ellps=WGS84 +datum=WGS84 +no_defs'
d =Domain(srsString, extentString)
# 2. Create nansatmap object from tha domain
nmap = Nansatmap(d)
# 3. Use the created object to convert from lon/lat into x/y:
x, y = nmap([0, 2, 4], [63, 64, 65])
# 4 or from x/y into lon/lat
lon, lat = nmap(x, y, inverse=True)


#n.write_map(grid='w', contour='w', smooth=true, quiver=('u','v'), step=5)
