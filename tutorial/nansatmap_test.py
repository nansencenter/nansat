#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import inspect, os

from nansat import Nansat, Domain, Nansatmap

#asumak
#iDir = 'c:/Users/asumak/Data/input/NCEP_GRIB/'
#oDir = 'c:/Users/asumak/Data/output/tutorial'

#antonk
iDir = '/Data/sat/GDAL_test/'
oDir = '/data/'

iFileName = os.path.join(iDir, 'gfs.t06z.master.grbf03')
oFileName = os.path.join(oDir, os.path.split(iFileName)[1])
print 'iFileName: ', iFileName
print 'oFileName: ', oFileName

# Create a Nansat object (n)
n = Nansat(iFileName)
# Get data from 1st, 2nd and 3rd bands as numpy array (u,v and w)
u = n[1]; v = n[2]; w = n[3]

# Create a Nansatmap object from the Nansat object (n)
nMap = Nansatmap(n)

# 1. Draw the contour_plots (line)
# 2. Draw continent
# 3. Save to a file
nMap.contour_plots(w)                                       # 1
nMap.process()                                              # 2
nMap.save(oFileName+'_contour1.png')                                 # 3
plt.close()
# or
nMap.process(contour_data=w)                                # 1 & 2
nMap.save(oFileName+'_contour2.png')                                 # 3
plt.close()

nMap.process(contour_data=w, contour_smoothing=True,
             contour_label=True, contour_linesfontsize=8,
             title='NCEP wind speed and direction',
             title_fontsize=15)                                     # 1 & 2
nMap.save(oFileName+'_contour3.png')                                 # 3
plt.close()

nMap.process(contour_data=w, contour_smoothing=True,
             contour_style='fill', colorbar_fontsize=8)             # 1 & 2
nMap.save(oFileName+'_contour4.png')                                 # 3
plt.close()

# 1. Make a pseudo-color plot over the map
# 2. Draw quiver_plots
# 3. Draw continent and geocoordinate grids and add the color bar
# 4. Show the map
nMap.put_color(w)                                                   # 1
nMap.quiver_plots([u, v])                                           # 2
nMap.process(geocoordinates=True, lat_fontsize=8, lon_fontsize=8)   # 3
nMap.save(oFileName+'_quiver1.png')                                 # 3
plt.close()

# or
nMap.process(color_data=w, quituiver_data=[u, v],
             geocoordinates=True, lat_fontsize=8, lon_fontsize=8)
nMap.save(oFileName+'_quiver2.png')
plt.close()
