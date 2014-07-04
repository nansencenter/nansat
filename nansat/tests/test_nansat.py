#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:04.07.2014 16:04
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import *

import nansat_test_archive as tna

class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noData:
            raise ValueError('No test data available')

    def test_mapper_imports(self):
        for folder in sys.path:
            for mapper in glob.glob(folder + '/mapper_*.py'):
                mm = os.path.basename(mapper.replace('.py',''))
                assert type(__import__(mm)) is ModuleType

    def test_pixel_functions(self):
        n = Nansat(self.test_data.asar[0])
        if sys.version_info < (2, 7):
            type(n['sigma0_VV']) == np.ndarray
        else:
            self.assertIsInstance(n['sigma0_VV'], np.ndarray)

    def tearDown(self):
        # This will delete the downloaded data after every test - it
        # seems to be the only solution if we don't want to keep the test data
        # Later we could make a data folder under tests and add that to
        # .gitignore
        self.test_data.delete_downloaded()

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
        self.test_data.delete_downloaded()



if __name__=='__main__':
    unittest.main()




