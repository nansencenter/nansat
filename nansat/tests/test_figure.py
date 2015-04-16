#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:16.03.2015 13:19
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

from nansat import Figure, Nansat, Domain
from nansat.tools import gdal, OptionError

import nansat_test_data as ntd

IS_CONDA = 'conda' in os.environ['PATH']


class FigureTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        plt.switch_backend('Agg')

        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')

    def test_init_array(self):
        f = Figure(np.zeros((10,10)), logLevel=40)

        self.assertEqual(type(f), Figure)

    def test_add_latlon_grids(self):
        n = Nansat(self.test_file_gcps)
        b = n[1]
        lon, lat = n.get_geolocation_grids()
        f = Figure(np.zeros((10,10)), logLevel=40)

        self.assertEqual(type(f), Figure)

if __name__ == "__main__":
    unittest.main()
