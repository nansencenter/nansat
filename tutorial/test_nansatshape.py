#!/usr/bin/env python
# Name:    nansatshape_test.py
# Purpose: Tutorial for nansatshape class
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

import os
import numpy as np

from nansat.nansatshape import Nansatshape
from nansat.nsr import NSR

# input and output file names
from testio import testio
iPath, oPath = testio()
oFileName = os.path.join(oPath, 'output_nansatshape_')
print 'Output file:', oFileName

''' Nansatshape class reads and writes ESRI-shape files

The core of Nansatshape is a OGR. the main functions of the class are
1. Create empty object in memory and add data (fields and geometory).
2. Open shape file and read the data.

Nansatshape support points

'''
# Create a nansatShape object with point geometry, lonlat spatial reference
nansatOGR = Nansatshape(srs=NSR())

# Create numpy array for coordinates of points
coordinates = np.array([[5.3,  30.2, 116.4, 34.0],
                        [60.4, 60.0,  39.9, 18.4]])

# Create structured numpy array for fieldValues
fieldValues = np.zeros(4, dtype={'names':   ['id', 'name', 'random'],
                                 'formats': ['i4', 'a10',  'f8']})
fieldValues['id'] = [1, 2, 3, 4]
fieldValues['name'] = ['Bergen', 'St.Petersburg', 'Beijin', 'Cape Town']
fieldValues['random'] = [1.3, 5.7, 11.13, 17.19]

# add fields to the feature and geometry at once
nansatOGR.add_features(fieldValues, coordinates)

# save to a file
nansatOGR.export(oFileName + '.shp')

# Create a nansatShape from shape file
nansatOGR2 = Nansatshape(oFileName + '.shp')
# Get corner points (geometries of featuers) in the layer
points, latlon = nansatOGR2.get_points()
# print corner points
print 'Corner Points ---'
for iPoint in points:
    print iPoint

print '\n*** nansatshape_test completed successfully. Output files are found here:' + oFileName
