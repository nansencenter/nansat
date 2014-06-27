#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:27.06.2014 11:00
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import *
import numpy as np

from nansat import *

import nansat_test_archive as tna
test_data = tna.TestData()

class NansatTest(unittest.TestCase):
    def setUp(self):
        if test_data.noData:
            raise ValueError('No test data available')

    def test_mapper_imports(self):
        for folder in sys.path:
            for mapper in glob.glob(folder + '/mapper_*.py'):
                mm = os.path.basename(mapper.replace('.py',''))
                assert type(__import__(mm)) is ModuleType

    def test_pixel_functions(self):
        n = Nansat(test_data.asar[0])
        if sys.version_info < (2, 7):
            type(n['sigma0_VV']) == np.ndarray
        else:
            self.assertIsInstance(n['sigma0_VV'], np.ndarray)

if __name__=='__main__':
    unittest.main()




