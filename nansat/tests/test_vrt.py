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
import logging
import os

import numpy as np
import gdal

from nansat.vrt import VRT
import nansat_test_data as ntd


class VRTTest(unittest.TestCase):
    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')

    def test_init(self):
        vrt = VRT()

        self.assertIsInstance(vrt, VRT)
        self.assertIsInstance(vrt.fileName, str)
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertIsInstance(vrt.logger, logging.Logger)
        self.assertIsInstance(vrt.driver, gdal.Driver)
        self.assertEqual(vrt.band_vrts, {})
        self.assertEqual(vrt.tps, False)
        self.assertTrue(vrt.vrt is None)
        self.assertTrue(vrt.xml.startswith('<VRTDataset rasterXSize="1" rasterYSize="1"'))

    def test_del(self):
        vrt = VRT()
        filename_vrt = vrt.fileName
        filename_raw = vrt.fileName.replace('.vrt', '.raw')
        vrt = None
        
        self.assertEqual(gdal.Unlink(filename_vrt), -1)
        self.assertEqual(gdal.Unlink(filename_raw), -1)

    def test_init_metadata(self):
        vrt1 = VRT(metadata={'aaa': 'bbb'})
        self.assertEqual(vrt1.dataset.GetMetadata()['aaa'], 'bbb')

    def test_init_nomem(self):
        vrt = VRT(nomem=True)

        self.assertIsInstance(vrt, VRT)
        self.assertIsInstance(vrt.fileName, str)
        self.assertTrue(os.path.exists(vrt.fileName))

    def test_from_gdal_dataset(self):
        ds = gdal.Open(self.test_file)
        vrt = VRT.from_gdal_dataset(ds)

        self.assertIsInstance(vrt, VRT)
        self.assertIsInstance(vrt.fileName, str)
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertEqual(vrt.dataset.RasterXSize, ds.RasterXSize)
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())

    def test_from_dataset_params(self):
        ds = gdal.Open(self.test_file)
        vrt = VRT.from_dataset_params(ds.RasterXSize,
                                      ds.RasterYSize,
                                      ds.GetGeoTransform(),
                                      ds.GetProjection(),
                                      ds.GetGCPs(),
                                      ds.GetGCPProjection())
                                      
        filename_vrt = vrt.fileName
        filename_raw = vrt.fileName.replace('.vrt', '.raw')
        self.assertIsInstance(vrt, VRT)
        self.assertIsInstance(filename_vrt, str)
        self.assertIsInstance(filename_raw, str)
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertEqual(vrt.dataset.RasterXSize, ds.RasterXSize)
        self.assertIsInstance(vrt.dataset.GetGCPs(), tuple)
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())
        
    def test_from_array(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)
                                      
        self.assertIsInstance(vrt, VRT)
        self.assertIsInstance(vrt.fileName, str)
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertEqual(vrt.dataset.RasterXSize, array.shape[1])
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())

    def test_from_lonlat(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)

        self.assertIsInstance(vrt, VRT)
        self.assertIsInstance(vrt.fileName, str)
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertEqual(vrt.dataset.RasterXSize, 10)
        self.assertEqual(vrt.dataset.RasterYSize, 30)
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())

    def test_copy_empty_vrt(self):
        vrt1 = VRT()
        vrt2 = vrt1.copy()

        self.assertIsInstance(vrt2, VRT)
        self.assertIsInstance(vrt2.fileName, str)

    def test_copy_vrt_with_band(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        vrt2 = vrt1.copy()

        self.assertIsInstance(vrt2, VRT)
        self.assertIsInstance(vrt2.fileName, str)
        self.assertEqual(vrt2.dataset.RasterCount, 1)

    def test_export(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)
        vrt.export('temp.vrt.xml')
        self.assertTrue(os.path.exists('temp.vrt.xml'))
