#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import inspect, os

from nansat import Nansat, Domain
from nansat_map import NansatMap

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

# Create a NansatMap object from the Nansat object (n)
nMap = NansatMap(n)

# 1. Draw the contour_plots (line)
# 2. Draw continent
# 3. Save to a file
# 4. Clear the figure
nMap.contour_plots(w)                                       # 1
nMap.process()                                              # 2
nMap.save(oFileName+'contour1.png')                         # 3
nMap.fig.clear()                                            # 4
# or
nMap.process(contour_data=w)                                # 1 & 2
nMap.save(oFileName+'contour2.png')                         # 3
nMap.fig.clear()                                            # 4

# 1. Smooth the data and draw the contour plots (line) and the labels
# 2. Draw continent and write the title
# 3. Save to a file
# 4. Clear the figure
nMap.contour_plots(w, contour_smoothing=True,
                   contour_label=True, contour_linesfontsize=8)     # 1
nMap.process(title='NCEP wind speed and direction',
             title_fontsize=15)                                     # 2
nMap.save(oFileName+'contour3.png')                                 # 3
nMap.fig.clear()                                                    # 4
# or
nMap.process(contour_data=w, contour_smoothing=True,
             contour_label=True, contour_linesfontsize=8,
             title='NCEP wind speed and direction',
             title_fontsize=15)                                     # 1 & 2
nMap.save(oFileName+'contour4.png')                                 # 3
nMap.fig.clear()                                                    # 4

# 1. Draw contour_plots (fill)
# 2. Draw continent and color bar
# 3. Save to a file
# 4. Clear the figure
nMap.contour_plots(w, contour_smoothing=True, contour_style='fill')     # 1
nMap.process(colorbar_fontsize=8)                                       # 2
nMap.save(oFileName+'contour_fill1.png')                                # 3
nMap.fig.clear()                                                        # 4
# or
nMap.process(contour_data=w, contour_smoothing=True,
             contour_style='fill', colorbar_fontsize=8)             # 1 & 2
nMap.save(oFileName+'contour_fill2.png')                            # 3
nMap.fig.clear()                                                    # 4

# 1. Make a pseudo-color plot over the map
# 2. Draw quiver_plots
# 3. Draw continent and geocoordinate grids and add the color bar
# 4. Save to a file
# 5. Clear the figure
nMap.put_color(w)                                                   # 1
nMap.quiver_plots([u, v])                                           # 2
nMap.process(geocoordinates=True, lat_fontsize=8, lon_fontsize=8)   # 3
nMap.save(oFileName+'color_quiver1.png')                            # 4
nMap.fig.clear()                                                    # 5
# or
nMap.process(color_data=w, quituiver_data=[u, v],
             geocoordinates=True, lat_fontsize=8, lon_fontsize=8)  # 1, 2 & 3
nMap.save(oFileName+'color_quiver2.png')                           # 4
nMap.fig.clear()                                                   # 5
