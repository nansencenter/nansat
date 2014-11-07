#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
import unittest
import warnings
import os
import sys
import glob
from types import ModuleType, FloatType
import numpy as np

from nansat import NSR
from nansat.tools import osr, ProjectionError


class NSRTest(unittest.TestCase):
    def test_init_empty(self):
        nsr = NSR()

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 0)

    def test_init_from_none(self):
        nsr = NSR(None)

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 5)

    def test_init_from_0(self):
        nsr = NSR(0)

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 0)

    def test_init_from_EPSG(self):
        nsr = NSR(4326)

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 0)
        self.assertTrue('4326' in nsr.ExportToWkt())

    def test_init_from_proj4(self):
        nsr = NSR('+proj=longlat')

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 0)
        self.assertTrue('longlat' in nsr.ExportToProj4())

    def test_init_from_wkt(self):
        nsr = NSR(osr.SRS_WKT_WGS84)

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 0)
        self.assertTrue('longlat' in nsr.ExportToProj4())

    def test_init_from_NSR(self):
        nsr = NSR(NSR(osr.SRS_WKT_WGS84))

        self.assertEqual(type(nsr), NSR)
        self.assertEqual(nsr.Validate(), 0)
        self.assertTrue('longlat' in nsr.ExportToProj4())

    def test_dont_init_from_invalid(self):
        self.assertRaises(ProjectionError, NSR, -10)
        self.assertRaises(ProjectionError, NSR, 'some crap')
