#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:     Anton Korosov
#
# Created:      18.06.2014
# Last modified:29.09.2014
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Nansat, Domain
from nansat.nansat import _import_mappers

import nansat_test_data as ntd


class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')

        if not os.path.exists(self.test_file):
            raise ValueError('No test data available')

    def test_resize_by_pixelsize(self):
        n = Nansat(self.test_file, logLevel=40)
        n.resize(pixelsize=500, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_factor(self):
        n = Nansat(self.test_file, logLevel=40)
        n.resize(0.5, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_width(self):
        n = Nansat(self.test_file, logLevel=40)
        n.resize(width=100, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_height(self):
        n = Nansat(self.test_file, logLevel=40)
        n.resize(height=500, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_reproject(self):
        n = Nansat(self.test_file, logLevel=40)
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-lle 15 -35 35 -25 -ts 500 500")
        n.reproject(d)
        self.assertEqual(type(n[1]), np.ndarray)

    def tearDown(self):
        # if any test plots are created, they could be deleted here
        pass
