#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Anton Korosov, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:  18.06.2014
# Last modified:27.08.2014 10:58
# Copyright:    (c) NERSC
# License:
#------------------------------------------------------------------------------
import unittest
from nansat import NSR
from nansat.utils import osr

from nansat.exceptions import NansatProjectionError


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

    def test_init_from_proj4_unicode(self):
        nsr = NSR(u'+proj=longlat')

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
        self.assertRaises(NansatProjectionError, NSR, -10)
        self.assertRaises(NansatProjectionError, NSR, 'some crap')
        ss = osr.SpatialReference()
        self.assertRaises(NansatProjectionError, NSR, ss)
