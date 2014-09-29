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


class DomainTest(unittest.TestCase):
    def setUp(self):
        self.test_data = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'data/gcps.tif')

        if not os.path.exists(self.test_data):
            raise ValueError('No test data available')

    def test_init_from_strings(self):
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-te 25 70 35 72 -ts 2000 2000")
        self.assertEqual(type(d), Domain)
