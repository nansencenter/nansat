#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Anton Korosov, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:18.11.2014 11:48
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
import unittest
import warnings
import os
import sys
import glob
from types import ModuleType, FloatType
import datetime

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from nansat import Nansat, Domain, Nansatmap, NSR
from nansat.tools import gdal, OptionError

import nansat_test_data as ntd

IS_CONDA = 'conda' in os.environ['PATH']


class NansatmapTest(unittest.TestCase):
    def setUp(self):
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        plt.switch_backend('Agg')

        if not os.path.exists(self.test_file_stere):
            raise ValueError('No test data available')

    def test_create_map(self):
        ''' should simply create a Nansatmap instance '''
        n = Nansat(self.test_file_stere, logLevel=40)
        nmap = Nansatmap(n)

        self.assertEqual(type(nmap), Nansatmap)

    def test_imshow(self):
        ''' Should use Nansatmap.imshow '''
        n = Nansat(self.test_file_stere, logLevel=40)
        b1 = n[1]
        nmap = Nansatmap(n)
        nmap.imshow(b1, cmap='ak01')
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_imshow.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_imshow_random(self):
        ''' Should use Nansatmap.imshow '''
        n = Nansat(self.test_file_stere, logLevel=40)
        b1 = n[1]
        nmap = Nansatmap(n)
        nmap.imshow(b1/5, cmap='random')
        nmap.add_colorbar()
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_imshow_random.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_pcolormesh(self):
        ''' Should use Nansatmap.pcolormesh '''
        n = Nansat(self.test_file_stere, logLevel=40)
        b1 = n[1]
        nmap = Nansatmap(n)
        nmap.pcolormesh(b1)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_pcolormesh.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_add_labels(self):
        size, npo = 100, 10
        xy = np.random.randint(0, size, npo*2).reshape(npo, 2)
        z = np.random.randint(0, size, npo)
        xg, yg = np.meshgrid(range(size), range(size))
        zg = griddata(xy, z, np.dstack([xg, yg]), method='nearest')
        dstDomain = Domain(NSR().wkt, '-te -10 -10 10 10 -ts 100 100')

        nmap = Nansatmap(dstDomain)
        nmap.imshow(zg, cmap='random')
        nmap.add_zone_labels(zg, fontsize=10)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_zonelables.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

if __name__ == "__main__":
    unittest.main()
