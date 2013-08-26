#!/usr/bin/env python
# Name:    mosaic_test.py
# Purpose: Tutorial for nansat.Mosaic
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      23.08.2013
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

from nansat import Nansat, Domain, Mosaic

from testio import testio
iPath, iFileName, oPath, oFileName = testio()

print iPath

# Create target domain
domain = Domain(4326, '-lle 27 70 31 72 -ts 1400 1300')

# A. Perform averaging of several files
# 1. Create destination Nansat object with desired projection
nMosaic = Mosaic(domain=domain)
# 2. Perfom averaging
nMosaic.average(['gcps.tif', 'stere.tif'], bands=['L_645', 'L_555', 'L_469'])
# 3. Get mask of valid pixels
mask = nMosaic['mask']
# 4. Output averaged data using the mask
nMosaic.write_figure(fileName=oFileName + '_mosaic.png', bands=['L_645', 'L_555', 'L_469'], clim='hist',
                        mask_array=mask, mask_lut={0:[128,128,128]})
# 5. Get values of standard deviation from averaging of input files
L_469_std = nMosaic['L_469_std']


# B. calculate median from the first band
nMosaic.median(['gcps.tif', 'stere.tif'])


# C. fill the result with the latest image
nMosaic.latest(['gcps.tif', 'stere.tif'])
