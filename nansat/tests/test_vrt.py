# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

from __future__ import absolute_import, unicode_literals
import unittest
import logging
import os
import sys
if sys.version_info.major == 2:
    from mock import patch
else:
    from unittest.mock import patch

import xml.etree.ElementTree as ET
import numpy as np
import gdal
import pythesint as pti

from nansat.node import Node
from nansat.vrt import VRT
from . import nansat_test_data as ntd


class VRTTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.test_file_complex = os.path.join(ntd.test_data_path, 'complex.nc')

        self.nsr_wkt = ('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",'
                        '6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUT'
                        'HORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORI'
                        'TY["EPSG","8901"]],UNIT["degree",0.0174532925199433'
                        ',AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')

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
        self.assertTrue(vrt.xml.startswith(b'<VRTDataset rasterXSize="1" rasterYSize="1"'))
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
        self.assertIn('filename', list(vrt.dataset.GetMetadata().keys()))

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
        self.assertIn('filename', list(vrt.dataset.GetMetadata().keys()))

    def test_from_array(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)

        self.assertEqual(vrt.dataset.RasterXSize, array.shape[1])
        self.assertEqual(vrt.dataset.RasterYSize, array.shape[0])
        self.assertEqual(vrt.dataset.RasterCount, 1)
        self.assertIn('filename', list(vrt.dataset.GetMetadata().keys()))
        self.assertEqual(gdal.Unlink(vrt.filename.replace('.vrt', '.raw')), 0)

    def test_from_lonlat(self):
        geo_keys = ['LINE_OFFSET', 'LINE_STEP', 'PIXEL_OFFSET', 'PIXEL_STEP', 'SRS',
                    'X_BAND', 'X_DATASET', 'Y_BAND', 'Y_DATASET']
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt = VRT.from_lonlat(lon, lat)

        self.assertEqual(vrt.dataset.RasterXSize, 10)
        self.assertEqual(vrt.dataset.RasterYSize, 30)
        self.assertIn('filename', list(vrt.dataset.GetMetadata().keys()))
        geo_metadata = vrt.dataset.GetMetadata(str('GEOLOCATION'))
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
        self.assertIn('filename', list(vrt2.dataset.GetMetadata().keys()))

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
        self.assertIn('rasterXSize', list(root.keys()))
        self.assertIn('rasterYSize', list(root.keys()))

        self.assertEqual([e.tag for e in root], ['Metadata', 'VRTRasterBand'])

    def test_create_band(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt1 = VRT.from_array(array)
        vrt2 = VRT(x_size=array.shape[1], y_size=array.shape[0])
        self.assertEqual(vrt2.dataset.RasterCount, 0)
        vrt2.create_band({'SourceFilename': vrt1.filename})
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
        ds = gdal.Open('NETCDF:"%s":UMass_AES' % self.test_file_arctic)
        vrt = VRT.copy_dataset(ds)
        self.assertTrue(b'<GeoTransform>' in vrt.xml)
        vrt._remove_geotransform()
        self.assertFalse(b'<GeoTransform>' in vrt.xml)

    def test_set_geotransform_for_resize(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt._set_geotransform_for_resize()

        self.assertEqual(vrt.dataset.GetMetadata(str('GEOLOCATION')), {})
        self.assertEqual(vrt.dataset.GetGCPs(), ())
        self.assertEqual(vrt.dataset.GetGeoTransform(), (0.0, 1.0, 0.0, 30, 0.0, -1.0))

    def test_set_gcps_geolocation_geotransform_with_geolocation(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt.create_band({str('SourceFilename'): vrt.geolocation.x_vrt.filename})
        vrt._set_gcps_geolocation_geotransform()
        self.assertFalse(b'<GeoTransform>' in vrt.xml)
        self.assertEqual(vrt.dataset.GetGCPs(), ())

    def test_set_gcps_geolocation_geotransform_with_gcps(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt.create_band({'SourceFilename': vrt.geolocation.x_vrt.filename})
        vrt._remove_geolocation()
        vrt._set_gcps_geolocation_geotransform()
        self.assertFalse(b'<GeoTransform>' in vrt.xml)
        self.assertIsInstance(vrt.dataset.GetGCPs(), (list, tuple))
        self.assertTrue(len(vrt.dataset.GetGCPs()) > 0)
        self.assertEqual(vrt.dataset.GetMetadata(str('GEOLOCATION')), {})

    def test_set_gcps_geolocation_geotransform_with_geotransform(self):
        ds = gdal.Open('NETCDF:"%s":UMass_AES' % self.test_file_arctic)
        vrt = VRT.copy_dataset(ds)
        vrt._set_gcps_geolocation_geotransform()
        self.assertEqual(vrt.dataset.GetGeoTransform(),
                         (-1000000.0, 25000.0, 0.0, 5000000.0, 0.0, -25000.0))
        self.assertEqual(vrt.dataset.GetMetadata(str('GEOLOCATION')), {})
        self.assertEqual(vrt.dataset.GetGCPs(), ())

    def test_update_warped_vrt_xml(self):
        dataset = gdal.Open('NETCDF:"%s":UMass_AES' % self.test_file_arctic)
        warped_dataset = gdal.AutoCreateWarpedVRT(dataset, None, str(self.nsr_wkt), 0)
        warped_vrt = VRT.copy_dataset(warped_dataset)
        x_size = 100
        y_size = 200
        geo_transform = (0.0, 1.0, 0.0, 200.0, 0.0, -1.0)
        block_size = 64
        working_data_type = 'Float32'
        warped_vrt._update_warped_vrt_xml(x_size, y_size, geo_transform, block_size,
                                          working_data_type)

        self.assertEqual(warped_vrt.dataset.RasterXSize, x_size)
        self.assertEqual(warped_vrt.dataset.RasterYSize, y_size)
        self.assertEqual(warped_vrt.dataset.GetGeoTransform(), geo_transform)
        self.assertEqual(warped_vrt.dataset.GetRasterBand(1).GetBlockSize(),
                         [block_size, block_size])
        self.assertIn(b'<WorkingDataType>Float32</WorkingDataType>', warped_vrt.xml)

    def test_set_fake_gcps_empty(self):
        ds = gdal.Open('NETCDF:"%s":UMass_AES' % self.test_file_arctic)
        vrt = VRT.copy_dataset(ds)

        dst_wkt = vrt._set_fake_gcps(self.nsr_wkt, [], 1)
        self.assertEqual(dst_wkt, self.nsr_wkt)
        self.assertEqual(len(vrt.dataset.GetGCPs()), 0)

    def test_set_fake_gcps(self):
        ds = gdal.Open('NETCDF:"%s":UMass_AES' % self.test_file_arctic)
        gcps = gdal.Open(self.test_file_gcps).GetGCPs()
        vrt = VRT.copy_dataset(ds)

        dst_wkt = vrt._set_fake_gcps(self.nsr_wkt, gcps, 1)
        self.assertEqual(dst_wkt, None)
        self.assertEqual(len(vrt.dataset.GetGCPs()), len(gcps))
        self.assertEqual([gcp.GCPPixel for gcp in gcps],
                         [gcp.GCPX for gcp in vrt.dataset.GetGCPs()])
        self.assertEqual([gcp.GCPLine for gcp in gcps],
                         [gcp.GCPY for gcp in vrt.dataset.GetGCPs()])

    def test_get_dst_band_data_type(self):
        self.assertEqual(VRT._get_dst_band_data_type([], {'dataType': 'Float32'}), 'Float32')
        self.assertEqual(VRT._get_dst_band_data_type([1, 2, 3], {}), gdal.GDT_Float32)
        self.assertEqual(VRT._get_dst_band_data_type([{'ScaleRatio': 2}], {}), gdal.GDT_Float32)
        self.assertEqual(VRT._get_dst_band_data_type([{'LUT': [1, 2, 3]}], {}), gdal.GDT_Float32)
        self.assertEqual(VRT._get_dst_band_data_type([{}], {}), gdal.GDT_Float32)
        self.assertEqual(VRT._get_dst_band_data_type([{'DataType': 'Float32'}], {}), 'Float32')

    def test_create_band_name(self):
        wkv = pti.get_wkv_variable('sigma0')
        ds = gdal.Open('NETCDF:"%s":UMass_AES' % self.test_file_arctic)
        vrt = VRT.copy_dataset(ds)
        self.assertEqual(vrt._create_band_name({'name': 'name1'}), ('name1', {}))
        self.assertEqual(vrt._create_band_name({'wkv': 'sigma0'}), ('sigma0', wkv))
        self.assertEqual(vrt._create_band_name({'wkv': 'sigma0', 'suffix': 'HH'}),
                         ('sigma0_HH', wkv))
        self.assertEqual(vrt._create_band_name({'name': 'UMass_AES'}),
                         ('UMass_AES_000', {}))

    def test_leave_few_bands(self):
        ds = gdal.Open(os.path.join(ntd.test_data_path, 'gcps.tif'))
        vrt = VRT.copy_dataset(ds)
        vrt.leave_few_bands([1, 'L_469'])
        self.assertEqual(vrt.dataset.RasterCount,2)
        self.assertEqual(vrt.dataset.GetRasterBand(1).GetMetadataItem('name'), 'L_645')
        self.assertEqual(vrt.dataset.GetRasterBand(2).GetMetadataItem('name'), 'L_469')

    def test_find_complex_band(self):
        a = np.random.randn(100,100)
        vrt1 = VRT.from_array(a)
        vrt2 = VRT.from_array(a.astype(np.complex64))

        vrt3 = VRT.from_gdal_dataset(vrt1.dataset)
        vrt3.create_bands([{'src':{'SourceFilename': vrt1.filename}},
                           {'src':{'SourceFilename': vrt2.filename}}])

        self.assertEqual(vrt1._find_complex_band(), None)
        self.assertEqual(vrt2._find_complex_band(), 1)
        self.assertEqual(vrt3._find_complex_band(), 2)

    def test_split_complex_bands(self):
        a = np.random.randn(100,100)
        vrt1 = VRT.from_array(a.astype(np.complex64))
        vrt2 = VRT.from_array(a)
        vrt3 = VRT.from_array(a.astype(np.complex64))

        vrt4 = VRT.from_gdal_dataset(vrt1.dataset)
        vrt4.create_bands([{'src':{'SourceFilename': vrt1.filename}, 'dst':{'name': 'vrt1'}},
                           {'src':{'SourceFilename': vrt2.filename}, 'dst':{'name': 'vrt2'}},
                           {'src':{'SourceFilename': vrt3.filename}, 'dst':{'name': 'vrt3'}},])

        vrt4.split_complex_bands()

        self.assertEqual(vrt4.dataset.RasterCount,5)
        self.assertEqual(vrt4.dataset.GetRasterBand(1).GetMetadataItem('name'), 'vrt2')
        self.assertEqual(vrt4.dataset.GetRasterBand(2).GetMetadataItem('name'), 'vrt1_real')
        self.assertEqual(vrt4.dataset.GetRasterBand(3).GetMetadataItem('name'), 'vrt1_imag')
        self.assertEqual(vrt4.dataset.GetRasterBand(4).GetMetadataItem('name'), 'vrt3_real')
        self.assertEqual(vrt4.dataset.GetRasterBand(5).GetMetadataItem('name'), 'vrt3_imag')

    def test_create_geolocation_bands(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt.create_geolocation_bands()

        self.assertEqual(vrt.dataset.RasterCount, 2)
        self.assertEqual(vrt.dataset.GetRasterBand(1).GetMetadataItem('name'), 'longitude')
        self.assertEqual(vrt.dataset.GetRasterBand(2).GetMetadataItem('name'), 'latitude')
        self.assertTrue(np.allclose(vrt.dataset.GetRasterBand(1).ReadAsArray(), lon))
        self.assertTrue(np.allclose(vrt.dataset.GetRasterBand(2).ReadAsArray(), lat))

    def test_fix_band_metadata(self):
        ds = gdal.Open(os.path.join(ntd.test_data_path, 'gcps.tif'))
        vrt = VRT.copy_dataset(ds)
        self.assertIn('standard_name', vrt.dataset.GetRasterBand(1).GetMetadata())
        self.assertIn('time', vrt.dataset.GetRasterBand(1).GetMetadata())
        vrt.fix_band_metadata(['standard_name', 'time'])
        self.assertNotIn('standard_name', vrt.dataset.GetRasterBand(1).GetMetadata())
        self.assertNotIn('time', vrt.dataset.GetRasterBand(1).GetMetadata())

    def test_fix_global_metadata(self):
        ds = gdal.Open(os.path.join(ntd.test_data_path, 'gcps.tif'))
        vrt = VRT.copy_dataset(ds)
        vrt.dataset.SetMetadataItem('test', '"test"')
        vrt.fix_global_metadata(['AREA_OR_POINT'])
        self.assertNotIn('AREA_OR_POINT', vrt.dataset.GetMetadata())
        self.assertEqual('&quot;test&quot;', vrt.dataset.GetMetadataItem('test'))

    def test_hardcopy_bands(self):
        ds = gdal.Open(os.path.join(ntd.test_data_path, 'gcps.tif'))
        vrt = VRT.copy_dataset(ds)
        vrt.hardcopy_bands()

        #import ipdb; ipdb.set_trace()
        self.assertTrue(np.allclose(vrt.dataset.ReadAsArray(), ds.ReadAsArray()))
        band_nodes = Node.create(vrt.xml).nodeList('VRTRasterBand')
        self.assertEqual(band_nodes[0].node('SourceFilename').value, vrt.band_vrts[1].filename)
        self.assertEqual(band_nodes[1].node('SourceFilename').value, vrt.band_vrts[2].filename)
        self.assertEqual(band_nodes[2].node('SourceFilename').value, vrt.band_vrts[3].filename)



if __name__ == "__main__":
    unittest.main()
