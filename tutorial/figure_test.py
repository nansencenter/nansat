#!/usr/bin/env python
# Name:    tutorial.py
# Purpose: Tutorial with list of examples
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
import inspect, os

from nansat import Nansat, Figure

# input and output file names
iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#iFileName = os.path.join(iPath, 'map.tif')
iFileName = os.path.join(iPath, 'gcps.tif')
print 'Input file: ', iFileName
oFileName = os.path.join(iPath, 'tmpdata', 'outputfig_')
print 'Output file prefix: ', oFileName

# Create a Nansat object (n)
n = Nansat(iFileName)
# get numpy array from the Nansat object
array = n[1]
print array

# Create a Figure object (fig)
fig = Figure(array)
# Set minimum and maximum values and add color-bar and title
fig.process(cmin=10, cmax=60)
# Save a figure
fig.save(oFileName+'01.png')

# Create a Figure object (fig)
fig = Figure(array)
# compute min and max valuse from ratio
clim = fig.clim_from_histogram(ratio=0.9)
# Set cmin and cmax values
fig.process(cmin=clim[0], cmax=clim[1], legend=True, titleString='NANSAT Tutorial', LEGEND_HEIGHT=0.3, fontSize=10)
# Save a figure
fig.save(oFileName+'02.png')

# create 3D array
array = [array, n[2], n[3]]
# Create a Figure object (fig) from 3D array
fig = Figure(array)
#fig.process(clim='hist', ratio=0.9)
#fig.process()
# make RGB image from bands 1,2,3 with brightness correction
fig.save(oFileName + '_rgb.png', bands=[1])

