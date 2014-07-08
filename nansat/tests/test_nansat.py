#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:08.07.2014 10:11
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Nansat, Domain
from nansat.nansat import nansatMappers 
import nansat_test_archive as tna

class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noData:
            raise ValueError('No test data available')

    def test_mapper_imports(self):
        for mapper in nansatMappers:
            assert type(__import__('nansat.mappers.'+mapper)) is ModuleType

    def test_pixel_functions(self):
        n = Nansat(self.test_data.asar[0])
        if sys.version_info < (2, 7):
            type(n['sigma0_VV']) == np.ndarray
        else:
            self.assertIsInstance(n['sigma0_VV'], np.ndarray)

    def test_resize_by_pixelsize(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(pixelsize=500, eResampleAlg=1)
        self.assertIsInstance(n[1], np.ndarray)

    def test_resize_by_factor(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(0.5, eResampleAlg=1)
        self.assertIsInstance(n[1], np.ndarray)

    def test_resize_by_width(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(width=100, eResampleAlg=1)
        self.assertIsInstance(n[1], np.ndarray)

    def test_resize_by_height(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(height=500, eResampleAlg=1)
        self.assertIsInstance(n[1], np.ndarray)

    def test_reproject(self):
        n = Nansat(self.test_data.asar[0])
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", 
                "-te 15 -35 35 -25 -ts 500 500") 
        n.reproject(d)
        self.assertIsInstance(n[1], np.ndarray)

    def tearDown(self):
        # if any test plots are created, they could be deleted here
        pass

class DomainTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noData:
            raise ValueError('No test data available')

    def test_get_pixelsize_meters(self):
        n = Nansat(self.test_data.asar[0])
        dX, dY = n.get_pixelsize_meters()
        self.assertIsInstance(dX, FloatType)
        self.assertIsInstance(dY, FloatType)

    def tearDown(self):
        #self.test_data.delete_downloaded()
        pass



if __name__=='__main__':
    unittest.main()




