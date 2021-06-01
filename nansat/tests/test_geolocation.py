#------------------------------------------------------------------------------
# Name:         test_vrt.py
# Purpose:      Test the VRT class
#
# Author:       Anton Korosov
#
# Created:      15.01.2018
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
import unittest
import os

import numpy as np

from nansat.vrt import VRT
from nansat.geolocation import Geolocation
from nansat.utils import gdal, osr
from nansat.tests import nansat_test_data as ntd

class GeolocationTest(unittest.TestCase):
    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')

    def test_init(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        x_vrt = VRT.from_array(lon)
        y_vrt = VRT.from_array(lat)

        ga = Geolocation(x_vrt, y_vrt)

        self.assertIsInstance(ga, Geolocation)
        self.assertEqual(ga.data['X_DATASET'], x_vrt.filename)
        self.assertEqual(ga.data['Y_DATASET'], y_vrt.filename)
        self.assertEqual(ga.data['LINE_OFFSET'], '0')
        self.assertEqual(ga.data['LINE_STEP'], '1')
        self.assertEqual(ga.data['PIXEL_OFFSET'], '0')
        self.assertEqual(ga.data['PIXEL_STEP'], '1')
        srs = osr.SpatialReference()
        status = srs.ImportFromWkt(ga.data['SRS'])
        self.assertEqual(status, 0)
        self.assertEqual(srs.ExportToProj4().strip(), '+proj=longlat +datum=WGS84 +no_defs')
        self.assertEqual(ga.data['X_BAND'], '1')
        self.assertEqual(ga.data['Y_BAND'], '1')
        self.assertEqual(ga.x_vrt, x_vrt)
        self.assertEqual(ga.y_vrt, y_vrt)

    def test_from_dataset(self):
        ds = gdal.Open(self.test_file)
        g = Geolocation.from_dataset(ds)
        self.assertIsInstance(g, Geolocation)

    def test_from_filenames(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        x_vrt = VRT.from_array(lon)
        y_vrt = VRT.from_array(lat)
        g = Geolocation.from_filenames(x_vrt.filename, y_vrt.filename)
        self.assertIsInstance(g, Geolocation)
        self.assertEqual(g.data['X_DATASET'], x_vrt.filename)
        self.assertEqual(g.data['Y_DATASET'], y_vrt.filename)
        self.assertEqual(g.data['LINE_OFFSET'], '0')
        self.assertEqual(g.data['LINE_STEP'], '1')
        self.assertEqual(g.data['PIXEL_OFFSET'], '0')
        self.assertEqual(g.data['PIXEL_STEP'], '1')
