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

from nansat.geolocation import Geolocation
import nansat_test_data as ntd

class Mock(object):
    pass

class GeolocationTest(unittest.TestCase):
    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')

    def test_init(self):
        x_vrt = Mock()
        x_vrt.fileName = 'aaa'
        y_vrt = Mock()
        y_vrt.fileName = 'bbb'

        ga = Geolocation(x_vrt, y_vrt)

        self.assertIsInstance(ga, Geolocation)

    def test_from_dataset(self):
        ds = gdal.Open(self.test_file)
        ga = Geolocation.from_dataset(ds)
        
        self.assertIsInstance(ga, Geolocation)
