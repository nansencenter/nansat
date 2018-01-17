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
from mock import patch

import numpy as np
import gdal

from nansat.vrt import VRT
import nansat_test_data as ntd


class VRTTest(unittest.TestCase):
    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')

    @patch.object(VRT, '_make_filename', return_value='/vsimem/filename.vrt')
    def test_init(self, _make_filename_mock):
        vrt = VRT()

        self.assertIsInstance(vrt, VRT)
        self.assertEqual(vrt.fileName, '/vsimem/filename.vrt')
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertIsInstance(vrt.logger, logging.Logger)
        self.assertIsInstance(vrt.driver, gdal.Driver)
        self.assertEqual(vrt.band_vrts, {})
        self.assertEqual(vrt.tps, False)
        self.assertTrue(vrt.vrt is None)
        self.assertTrue(vrt.xml.startswith('<VRTDataset rasterXSize="1" rasterYSize="1"'))
        _make_filename_mock.called_once()

    def test_del(self):
        vrt = VRT()
        self.assertEqual(gdal.Unlink(vrt.fileName), 0)

        vrt = VRT()
        filename_vrt = vrt.fileName
        vrt = None
        self.assertEqual(gdal.Unlink(filename_vrt), -1)

    def test_init_metadata(self):
        vrt1 = VRT(metadata={'aaa': 'bbb'})
        self.assertEqual(vrt1.dataset.GetMetadata()['aaa'], 'bbb')

    def test_init_nomem(self):
        vrt = VRT(nomem=True)

        self.assertTrue(os.path.exists(vrt.fileName))

    def test_from_gdal_dataset(self):
        ds = gdal.Open(self.test_file)
        vrt = VRT.from_gdal_dataset(ds)

        self.assertEqual(vrt.dataset.RasterXSize, ds.RasterXSize)
        self.assertEqual(vrt.dataset.RasterYSize, ds.RasterYSize)
        self.assertEqual(vrt.dataset.GetProjection(), ds.GetProjection())
        self.assertEqual(vrt.dataset.GetGeoTransform(), ds.GetGeoTransform())
        self.assertEqual(vrt.dataset.GetGCPProjection(), ds.GetGCPProjection())
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())

    def test_from_dataset_params(self):
        ds = gdal.Open(self.test_file)
        vrt = VRT.from_dataset_params(ds.RasterXSize,
                                      ds.RasterYSize,
                                      ds.GetGeoTransform(),
                                      ds.GetProjection(),
                                      ds.GetGCPs(),
                                      ds.GetGCPProjection())

        self.assertEqual(vrt.dataset.RasterXSize, ds.RasterXSize)
        self.assertEqual(vrt.dataset.RasterYSize, ds.RasterYSize)
        self.assertEqual(vrt.dataset.GetProjection(), ds.GetProjection())
        self.assertEqual(vrt.dataset.GetGeoTransform(), ds.GetGeoTransform())
        self.assertEqual(vrt.dataset.GetGCPProjection(), ds.GetGCPProjection())
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())

    def test_from_array(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)

        self.assertEqual(vrt.dataset.RasterXSize, array.shape[1])
        self.assertEqual(vrt.dataset.RasterYSize, array.shape[0])
        self.assertEqual(vrt.dataset.RasterCount, 1)
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())
        self.assertEqual(gdal.Unlink(vrt.fileName.replace('.vrt', '.raw')), 0)

    def test_from_lonlat(self):
        geo_keys = ['LINE_OFFSET', 'LINE_STEP', 'PIXEL_OFFSET', 'PIXEL_STEP', 'SRS',
                    'X_BAND', 'X_DATASET', 'Y_BAND', 'Y_DATASET']
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)

        self.assertEqual(vrt.dataset.RasterXSize, 10)
        self.assertEqual(vrt.dataset.RasterYSize, 30)
        self.assertIn('fileName', vrt.dataset.GetMetadata().keys())
        geo_metadata = vrt.dataset.GetMetadata('GEOLOCATION')
        for geo_key in geo_keys:
            self.assertEqual(vrt.geolocation.data[geo_key], geo_metadata[geo_key])
        self.assertIsInstance(vrt.geolocation.x_vrt, VRT)
        self.assertIsInstance(vrt.geolocation.y_vrt, VRT)
        self.assertEqual(vrt.geolocation.x_vrt.fileName, geo_metadata['X_DATASET'])
        self.assertEqual(vrt.geolocation.y_vrt.fileName, geo_metadata['Y_DATASET'])

    def test_copy_empty_vrt(self):
        vrt1 = VRT()
        vrt2 = vrt1.copy()

        self.assertIsInstance(vrt2, VRT)
        self.assertIsInstance(vrt2.fileName, str)
        self.assertEqual(vrt2.dataset.RasterXSize, vrt1.dataset.RasterXSize)
        self.assertEqual(vrt2.dataset.RasterYSize, vrt1.dataset.RasterYSize)
        self.assertEqual(vrt2.dataset.GetProjection(), vrt1.dataset.GetProjection())
        self.assertEqual(vrt2.dataset.GetGeoTransform(), vrt1.dataset.GetGeoTransform())
        self.assertEqual(vrt2.dataset.GetGCPProjection(), vrt1.dataset.GetGCPProjection())
        self.assertIn('fileName', vrt2.dataset.GetMetadata().keys())

    def test_copy_vrt_with_band(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        vrt2 = vrt1.copy()

        self.assertEqual(vrt2.dataset.RasterCount, 1)

    def test_export(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)
        vrt.export('temp.vrt.xml')
        self.assertTrue(os.path.exists('temp.vrt.xml'))

    def test_create_band(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        vrt2 = VRT(x_size=array.shape[1], y_size=array.shape[0])
        self.assertEqual(vrt2.dataset.RasterCount, 0)
        vrt2._create_band({'SourceFilename': vrt1.fileName})
        self.assertEqual(vrt2.dataset.RasterCount, 1)

    def test_make_source_bands_xml(self):
        array = gdal.Open(self.test_file).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        src1 = {'SourceFilename': vrt1.fileName}
        src2 = VRT._make_source_bands_xml(src1)
        self.assertIn('XML', src2)
        self.assertEqual(src2['SourceFilename'], vrt1.fileName)
        self.assertEqual(src2['SourceBand'], 1)
        self.assertEqual(src2['LUT'], '')
        self.assertEqual(src2['NODATA'], '')
        self.assertEqual(src2['SourceType'], 'ComplexSource')
        self.assertEqual(src2['ScaleRatio'], 1.0)
        self.assertEqual(src2['ScaleOffset'], 0.0)
        self.assertEqual(src2['DataType'], 1)
        self.assertEqual(src2['xSize'], 200)
        self.assertEqual(src2['ySize'], 190)

        with self.assertRaises(AttributeError):
            src2 = VRT._make_source_bands_xml({})

    def test_set_add_band_options(self):
        # case 1
        srcs = [{'SourceFilename': 'filename', 'SourceBand': 1}]
        dst = []
        options = VRT._set_add_band_options(srcs, dst)
        self.assertEqual(options, [])
        # case 2
        srcs = [{'SourceFilename': 'filename',
                 'SourceBand': 0,
                 'ImageOffset': 0,
                 'PixelOffset': 0,
                 'LineOffset': 0,
                 'ByteOrder': 'i'}]
        options = VRT._set_add_band_options(srcs, dst)
        self.assertIn('subclass=VRTRawRasterBand', options)
        self.assertIn('SourceFilename=filename', options)
        self.assertIn('ImageOffset=0', options)
        self.assertIn('PixelOffset=0', options)
        self.assertIn('LineOffset=0', options)
        self.assertIn('ByteOrder=i', options)


if __name__ == "__main__":
    unittest.main()
