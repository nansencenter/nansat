#!/usr/bin/env python
# Name:    figure_test.py
# Purpose: Tutorial for figrue class
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
import os

from nansat import Nansat, Figure

# input and output file names
from testio import testio
iPath, oPath = testio()
iFileName = os.path.join(iPath, 'gcps.tif')
print 'Input file: ', iFileName
oFileName = os.path.join(oPath, 'output_figure_')
print 'Output file:', oFileName

'''Figrue class creates figure(png, jpg, tif, bmp) from numpy array

Figure object is created by numpy array.
This class consequently create a figrue from 2D numpy array (1band) or
an RGB figure from 3D numpy array (3bands)

This class has fllolowing methods:
    estimate min/max, apply logarithmic scaling, convert to uint8,
    append legend, save to a file

'''

# Create a Nansat object (n)
n = Nansat(iFileName)
# get numpy array from the Nansat object
array = n[2]

# Create a Figure object (fig)
fig = Figure(array)
# Set minimum and maximum values
fig.process(cmin=10, cmax=60)
# Save the figure
fig.save(oFileName + '01_clim.png')

# Create a Figure object (fig)
fig = Figure(array)
# Set minimum and maximum values and custom colormap
fig.process(cmin=10, cmax=60, cmapName='obpg')
# Save the figure
fig.save(oFileName + '01_clim_obpg.png')

# Create a Figure object (fig)
fig = Figure(array)
# Compute min and max valuse from histogram
clim = fig.clim_from_histogram(ratio=0.7)
# Set cmin and cmax values
fig.process(cmin=clim[0], cmax=clim[1])
# Save the figure
fig.save(oFileName + '02_clim.png')

# Create a Figure object (fig)
fig = Figure(array)
# Compute min and max valuse from histogram
clim = fig.clim_from_histogram(ratio=0.7)
# Set cmin and cmax values
fig.process(cmin=10, cmax=60, logarithm=True, gamma=2)
# Save the figure
fig.save(oFileName + '03_logscale.png')

# Create a Figure object (fig)
fig = Figure(array)
# Make indexed image with legend
fig.process(cmin=10, cmax=60, legend=True, titleString='NANSAT figure_test',
         LEGEND_HEIGHT=0.3, fontSize=10)
# Save the figure
fig.save(oFileName + '04_title.png')

# Create a Figure object (fig)
fig = Figure(array)
# add logo to image to the lower left corner
# (make sure file is in the current folder)
fig.process(cmin=10, cmax=60, logoFileName='nansat_logo_s.png',
            logoLocation=[10, -35], logoSize=[20, 20],
            legend=True, LEGEND_HEIGHT=0.3)
# Save the figure
fig.save(oFileName + '05_logo.png')

# Create a Figure object (fig)
fig = Figure(array)
# Get lat/lon arrays from Nansat object (may take some time)
lonGrid, latGrid = n.get_geolocation_grids()
# Make figure with lat/lon grids
fig.process(cmin=10, cmax=60, latGrid=latGrid, lonGrid=lonGrid,
            latlonGridSpacing=10, latlonLabels=10)
# save the fig
fig.save(oFileName + '06_latlon.png', )

# Create a Figure object (fig)
fig = Figure(array)
# Get Nansat object with watermask
wm = n.watermask()
# Get array from Nansat object. 0 - land, 1 - water
wmArray = wm[1]
# Compute min and max valuse from ratio using non masked pixels only
clim = fig.clim_from_histogram(ratio=0.8, mask_array=wmArray,
                                          mask_lut={2: [128, 128, 128]})
# Make figure with land overlay (gray) and apply brightness logarithm correction
fig.process(cmin=clim[0], cmax=clim[1])
# save the fig
fig.save(oFileName + '07_land.png')

# create 3D numpy array
array = np.array([n[1], n[2], n[3]])
# Create a Figure object (fig) from 3D array
fig = Figure(array)
# Compute min and max valuse from ratio
clim = fig.clim_from_histogram(ratio=0.9)
# Set cmin and cmax values
fig.process(cmin=clim[0], cmax=clim[1])
# make RGB image from bands 1,2,3 with brightness correction
fig.save(oFileName + '08_rgb.png', bands=[1, 2, 3])

print '\n***figure_test completed successfully. Output files are found here:' + oFileName
