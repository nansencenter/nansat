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
import xml.etree.ElementTree as ET

import numpy as np
import gdal

from nansat.vrt import VRT
import nansat_test_data as ntd


class VRTTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')

    @patch.object(VRT, '_make_filename', return_value='/vsimem/filename.vrt')
    def test_init(self, _make_filename_mock):
        vrt = VRT()

        self.assertIsInstance(vrt, VRT)
        self.assertEqual(vrt.filename, '/vsimem/filename.vrt')
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
        self.assertEqual(gdal.Unlink(vrt.filename), 0)

        vrt = VRT()
        filename_vrt = vrt.filename
        vrt = None
        self.assertEqual(gdal.Unlink(filename_vrt), -1)

    def test_init_metadata(self):
        vrt1 = VRT(metadata={'aaa': 'bbb'})
        self.assertEqual(vrt1.dataset.GetMetadata()['aaa'], 'bbb')

    def test_init_nomem(self):
        vrt = VRT(nomem=True)

        self.assertTrue(os.path.exists(vrt.filename))

    def test_from_gdal_dataset(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt = VRT.from_gdal_dataset(ds)

        self.assertEqual(vrt.dataset.RasterXSize, ds.RasterXSize)
        self.assertEqual(vrt.dataset.RasterYSize, ds.RasterYSize)
        self.assertEqual(vrt.dataset.GetProjection(), ds.GetProjection())
        self.assertEqual(vrt.dataset.GetGeoTransform(), ds.GetGeoTransform())
        self.assertEqual(vrt.dataset.GetGCPProjection(), ds.GetGCPProjection())
        self.assertIn('filename', vrt.dataset.GetMetadata().keys())

    def test_from_dataset_params(self):
        ds = gdal.Open(self.test_file_gcps)
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
        self.assertIn('filename', vrt.dataset.GetMetadata().keys())

    def test_from_array(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)

        self.assertEqual(vrt.dataset.RasterXSize, array.shape[1])
        self.assertEqual(vrt.dataset.RasterYSize, array.shape[0])
        self.assertEqual(vrt.dataset.RasterCount, 1)
        self.assertIn('filename', vrt.dataset.GetMetadata().keys())
        self.assertEqual(gdal.Unlink(vrt.filename.replace('.vrt', '.raw')), 0)

    def test_from_lonlat(self):
        geo_keys = ['LINE_OFFSET', 'LINE_STEP', 'PIXEL_OFFSET', 'PIXEL_STEP', 'SRS',
                    'X_BAND', 'X_DATASET', 'Y_BAND', 'Y_DATASET']
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)

        self.assertEqual(vrt.dataset.RasterXSize, 10)
        self.assertEqual(vrt.dataset.RasterYSize, 30)
        self.assertIn('filename', vrt.dataset.GetMetadata().keys())
        geo_metadata = vrt.dataset.GetMetadata('GEOLOCATION')
        for geo_key in geo_keys:
            self.assertEqual(vrt.geolocation.data[geo_key], geo_metadata[geo_key])
        self.assertIsInstance(vrt.geolocation.x_vrt, VRT)
        self.assertIsInstance(vrt.geolocation.y_vrt, VRT)
        self.assertEqual(vrt.geolocation.x_vrt.filename, geo_metadata['X_DATASET'])
        self.assertEqual(vrt.geolocation.y_vrt.filename, geo_metadata['Y_DATASET'])

    def test_copy_empty_vrt(self):
        vrt1 = VRT()
        vrt2 = vrt1.copy()

        self.assertIsInstance(vrt2, VRT)
        self.assertIsInstance(vrt2.filename, str)
        self.assertEqual(vrt2.dataset.RasterXSize, vrt1.dataset.RasterXSize)
        self.assertEqual(vrt2.dataset.RasterYSize, vrt1.dataset.RasterYSize)
        self.assertEqual(vrt2.dataset.GetProjection(), vrt1.dataset.GetProjection())
        self.assertEqual(vrt2.dataset.GetGeoTransform(), vrt1.dataset.GetGeoTransform())
        self.assertEqual(vrt2.dataset.GetGCPProjection(), vrt1.dataset.GetGCPProjection())
        self.assertIn('filename', vrt2.dataset.GetMetadata().keys())

    def test_copy_vrt_with_band(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        vrt2 = vrt1.copy()

        self.assertEqual(vrt2.dataset.RasterCount, 1)

    def test_export(self):
        tmpfilename = os.path.join(ntd.tmp_data_path, 'temp.vrt.xml')
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)
        vrt.export(tmpfilename)
        self.assertTrue(tmpfilename)
        tree = ET.parse(tmpfilename)
        root = tree.getroot()

        self.assertEqual(root.tag, 'VRTDataset')
        self.assertEqual(root.keys(), ['rasterXSize', 'rasterYSize'])
        self.assertEqual([e.tag for e in root], ['Metadata', 'VRTRasterBand'])

    def test_create_band(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        vrt2 = VRT(x_size=array.shape[1], y_size=array.shape[0])
        self.assertEqual(vrt2.dataset.RasterCount, 0)
        vrt2._create_band({'SourceFilename': vrt1.filename})
        self.assertEqual(vrt2.dataset.RasterCount, 1)

    def test_make_source_bands_xml(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        src1 = {'SourceFilename': vrt1.filename}
        src2 = VRT._make_source_bands_xml(src1)
        self.assertIn('XML', src2)
        self.assertEqual(src2['SourceFilename'], vrt1.filename)
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

    def test_remove_geotransform(self):
        ds = gdal.Open('NETCDF:"%s":UMass_AES'%self.test_file_arctic)
        vrt = VRT.copy_dataset(ds)
        self.assertTrue('<GeoTransform>' in vrt.xml)
        vrt._remove_geotransform()
        self.assertFalse('<GeoTransform>' in vrt.xml)

    def test_set_geotransform_for_resize(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt._set_geotransform_for_resize()

        self.assertEqual(vrt.dataset.GetMetadata('GEOLOCATION'), {})
        self.assertEqual(vrt.dataset.GetGCPs(), ())
        self.assertEqual(vrt.dataset.GetGeoTransform(), (0.0, 1.0, 0.0, 30, 0.0, -1.0))

    def test_set_gcps_geolocation_geotransform_with_geolocation(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt._create_band({'SourceFilename': vrt.geolocation.x_vrt.filename})
        vrt._set_gcps_geolocation_geotransform()
        self.assertFalse('<GeoTransform>' in vrt.xml)
        self.assertEqual(vrt.dataset.GetGCPs(), ())

    def test_set_gcps_geolocation_geotransform_with_gcps(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt._create_band({'SourceFilename': vrt.geolocation.x_vrt.filename})
        vrt._remove_geolocation()
        vrt._set_gcps_geolocation_geotransform()
        self.assertFalse('<GeoTransform>' in vrt.xml)
        self.assertIsInstance(vrt.dataset.GetGCPs(), (list, tuple))
        self.assertTrue(len(vrt.dataset.GetGCPs()) > 0)
        self.assertEqual(vrt.dataset.GetMetadata('GEOLOCATION'), {})

    def test_set_gcps_geolocation_geotransform_with_geotransform(self):
        ds = gdal.Open('NETCDF:"%s":UMass_AES'%self.test_file_arctic)
        vrt = VRT.copy_dataset(ds)
        vrt._set_gcps_geolocation_geotransform()
        self.assertEqual(vrt.dataset.GetGeoTransform(),
                        (-1000000.0, 25000.0, 0.0, 5000000.0, 0.0, -25000.0))
        self.assertEqual(vrt.dataset.GetMetadata('GEOLOCATION'), {})
        self.assertEqual(vrt.dataset.GetGCPs(), ())

    def test_update_warped_vrt_xml(self):
        pass

    def test_set_fake_gcps(self):
        pass

if __name__ == "__main__":
    unittest.main()
