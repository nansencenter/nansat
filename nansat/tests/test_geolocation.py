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

import gdal
import numpy as np

from nansat.vrt import VRT
from nansat.geolocation import Geolocation
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
        self.assertEqual(ga.data['SRS'], 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",'
                                         '6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUT'
                                         'HORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORI'
                                         'TY["EPSG","8901"]],UNIT["degree",0.0174532925199433'
                                         ',AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
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
