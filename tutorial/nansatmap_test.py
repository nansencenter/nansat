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




#n.write_map(grid='w', contour='w', smooth=true, quiver=('u','v'), step=5)
