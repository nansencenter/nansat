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
from scipy.io.netcdf import netcdf_file

from nansat import Nansat, Domain, Nansatmap
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
        nmap.imshow(b1)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_imshow.png')
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

if __name__ == "__main__":
    unittest.main()
