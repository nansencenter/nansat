#!/usr/bin/env python
# Name:         tutorial.py
# Author:       Anton Korosov
# Licence:      GPL v.3
# This file is part of NANSAT. You can redistribute it or modify
# under the terms of GNU General Public License, v.3
# http://www.gnu.org/licenses/gpl-3.0.html

import matplotlib.pyplot as plt
import os

from nansat import Nansat, Domain
from testio import testio

# Set input and output file names
iPath, oPath = testio()
iFileName1 = os.path.join(iPath, 'gcps.tif')
iFileName2 = os.path.join(iPath, 'stere.tif')
oFileName = os.path.join(oPath, 'tutor_')

# Open an input satellite image with Nansat
n1 = Nansat(iFileName1, logLevel=10)
n2 = Nansat(iFileName2, logLevel=10)

n1.vrt.export('n1.vrt.txt')

n1.write_figure(oFileName + 'n1.png', clim=[0, 60])
n1.resize(0.1, eResampleAlg=3)
n1.resize(10, eResampleAlg=3)
n1.write_figure(oFileName + 'n1x01x10.png', clim=[0, 60])
n1.undo(100)

n2.write_figure(oFileName + 'n2.png', clim=[0, 60])
n2.reproject(n1, eResampleAlg=2)
n2.write_figure(oFileName + 'n2_on_n1.png', clim=[0, 60])


n1.resize(2)
n1.write_figure(oFileName + 'n1x2.png', clim=[0, 60])
n2.reproject(n1, eResampleAlg=2)
n2.write_figure(oFileName + 'n2x10_on_n1x2.png', clim=[0, 60])

n2.undo(100)
n2.resize(0.2, eResampleAlg=-1)
n2.write_figure(oFileName + 'n2x05.png', clim=[0, 60])
n2.reproject(n1, eResampleAlg=2)
n2.write_figure(oFileName + 'n2x05_on_n1x20.png', clim=[0, 60])
