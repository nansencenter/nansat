#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:  18.06.2014
# Last modified:27.08.2014 10:58
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Nansat, Domain
from nansat.nansat import _import_mappers


class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_data = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'data/gcps.tif')

        if not os.path.exists(self.test_data):
            raise ValueError('No test data available')

    def test_resize_by_pixelsize(self):
        n = Nansat(self.test_data, logLevel=40)
        n.resize(pixelsize=500, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_factor(self):
        n = Nansat(self.test_data, logLevel=40)
        n.resize(0.5, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_width(self):
        n = Nansat(self.test_data, logLevel=40)
        n.resize(width=100, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_height(self):
        n = Nansat(self.test_data, logLevel=40)
        n.resize(height=500, eResampleAlg=1)
        self.assertEqual(type(n[1]), np.ndarray)

    def test_reproject(self):
        n = Nansat(self.test_data, logLevel=40)
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-lle 15 -35 35 -25 -ts 500 500")
        n.reproject(d)
        self.assertEqual(type(n[1]), np.ndarray)

    def tearDown(self):
        # if any test plots are created, they could be deleted here
        pass
