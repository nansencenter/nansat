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
from mock import patch, PropertyMock, Mock, MagicMock, DEFAULT

import xml.etree.ElementTree as ET
import warnings

import numpy as np
import gdal
import pythesint as pti

from nansat.node import Node
from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.tests.nansat_test_base import NansatTestBase

from nansat.exceptions import NansatProjectionError


class VRTTest(NansatTestBase):
    nsr_wkt = ('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",'
                '6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUT'
                'HORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORI'
                'TY["EPSG","8901"]],UNIT["degree",0.0174532925199433'
                ',AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')

    @patch.object(VRT, '_make_filename', return_value='/vsimem/filename.vrt')
    def test_init(self, mock_make_filename):
        metadata={'key': 'value'}
        vrt = VRT(metadata=metadata)

        self.assertIsInstance(vrt, VRT)
        self.assertEqual(vrt.filename, '/vsimem/filename.vrt')
        self.assertIsInstance(vrt.dataset, gdal.Dataset)
        self.assertIsInstance(vrt.logger, logging.Logger) # just for testing mocking
        self.assertIsInstance(vrt.driver, gdal.Driver)
        self.assertEqual(vrt.band_vrts, {})
        self.assertEqual(vrt.tps, False)
        self.assertTrue(vrt.vrt is None)
        self.assertTrue(vrt.xml.startswith('<VRTDataset rasterXSize="1" rasterYSize="1"'))
        self.assertTrue(mock_make_filename.called_once())
        self.assertEqual(vrt.dataset.GetMetadata(), metadata)

    @patch.object(VRT, '_make_filename', return_value='filename.vrt')
    def test_del(self, mock_make_filename):
        vrt = VRT()
        vrt = None
        self.assertFalse(os.path.exists('filename.vrt'))

    @patch.object(VRT, '_init_from_gdal_dataset')
    def test_from_gdal_dataset(self, _init_from_gdal_dataset):
        ds = gdal.Open(self.test_file_gcps)
        vrt = VRT.from_gdal_dataset(ds)
        self.assertIsInstance(vrt, VRT)
        self.assertTrue(_init_from_gdal_dataset.called_once())

    @patch.object(VRT, '_add_geolocation')
    def test_init_from_gdal_dataset(self, _add_geolocation):
        vrt = VRT()
        ds = gdal.Open(self.test_file_gcps)
        vrt._init_from_gdal_dataset(ds)

        self.assertEqual(vrt.dataset.RasterXSize, ds.RasterXSize)
        self.assertEqual(vrt.dataset.RasterYSize, ds.RasterYSize)
        self.assertEqual(vrt.dataset.GetProjection(), ds.GetProjection())
        self.assertEqual(vrt.dataset.GetGeoTransform(), ds.GetGeoTransform())
        self.assertEqual(vrt.dataset.GetGCPProjection(), ds.GetGCPProjection())
        self.assertIn('filename', list(vrt.dataset.GetMetadata().keys()))
        self.assertTrue(_add_geolocation.called_once())

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
        vrt = VRT.from_lonlat(lon, lat, n_gcps=25)

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
        self.assertEqual(len(vrt.dataset.GetGCPs()), 25)

    def test_from_lonlat_no_gcps(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt = VRT.from_lonlat(lon, lat, add_gcps=False)
        self.assertEqual(len(vrt.dataset.GetGCPs()), 0)

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

    def test_copy_vrt_pixel_func(self):
        vrt1 = VRT()
        vrt1_xml = '''
        <VRTDataset rasterXSize="200" rasterYSize="200">
            <VRTRasterBand dataType="Byte" band="1">
                <ComplexSource>
                    <SourceFilename relativeToVRT="0">%s</SourceFilename>
                    <SourceBand>1</SourceBand>
                    <SourceProperties RasterXSize="200" RasterYSize="200" DataType="Byte" BlockXSize="200" BlockYSize="13" />
                    <SrcRect xOff="0" yOff="0" xSize="200" ySize="200" />
                    <DstRect xOff="0" yOff="0" xSize="200" ySize="200" />
                </ComplexSource>
            </VRTRasterBand>
            <VRTRasterBand dataType="Float32" band="2" subClass="VRTDerivedRasterBand">
                <ComplexSource>
                    <SourceFilename relativeToVRT="0">%s</SourceFilename>
                    <SourceBand>1</SourceBand>
                    <SourceProperties RasterXSize="200" RasterYSize="200" DataType="Byte" BlockXSize="128" BlockYSize="128" />
                    <SrcRect xOff="0" yOff="0" xSize="200" ySize="200" />
                    <DstRect xOff="0" yOff="0" xSize="200" ySize="200" />
                </ComplexSource>
                <PixelFunctionType>sqrt</PixelFunctionType>
            </VRTRasterBand>
        </VRTDataset>
       ''' % (self.test_file_gcps, vrt1.filename)
        vrt1.write_xml(vrt1_xml)
        vrt2 = vrt1.copy()

        self.assertFalse(os.path.basename(vrt1.filename) in vrt2.xml)

    def test_copy_geolocation(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt1 = VRT.from_lonlat(lon, lat)

        vrt2 = vrt1.copy()
        vrt1 = None
        self.assertTrue(vrt2.geolocation.x_vrt is not None)
        self.assertTrue(vrt2.geolocation.y_vrt is not None)

    def test_export(self):
        array = gdal.Open(self.test_file_gcps).ReadAsArray()[1, 10:, :]
        vrt = VRT.from_array(array)
        vrt.export(self.tmp_filename)
        self.assertTrue(self.tmp_filename)
        tree = ET.parse(self.tmp_filename)
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

        with self.assertRaises(KeyError):
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
        self.assertTrue('<GeoTransform>' in vrt.xml)
        vrt._remove_geotransform()
        self.assertFalse('<GeoTransform>' in vrt.xml)

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
        self.assertFalse('<GeoTransform>' in vrt.xml)
        self.assertEqual(vrt.dataset.GetGCPs(), ())

    def test_set_gcps_geolocation_geotransform_with_gcps(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt.create_band({'SourceFilename': vrt.geolocation.x_vrt.filename})
        vrt._remove_geolocation()
        vrt._set_gcps_geolocation_geotransform()
        self.assertFalse('<GeoTransform>' in vrt.xml)
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
        self.assertIn('<WorkingDataType>Float32</WorkingDataType>', warped_vrt.xml)

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

    def test_create_band_name_no_wkv(self):
        self.mock_pti['get_wkv_variable'].side_effect = IndexError
        vrt = VRT()
        self.assertEqual(vrt._create_band_name({'name': 'name1'}), ('name1', {}))

    def test_create_band_name_wkv(self):
        short_name='sigma0'
        wkv = dict(short_name=short_name)
        self.mock_pti['get_wkv_variable'].return_value=wkv
        vrt = VRT()
        self.assertEqual(vrt._create_band_name({'wkv': short_name}), (short_name, wkv))
        self.assertEqual(vrt._create_band_name({'wkv': short_name, 'suffix': 'HH'}),
                         (short_name + '_HH', wkv))

    def test_create_band_name_existing_name(self):
        self.mock_pti['get_wkv_variable'].side_effect = IndexError
        vrt = VRT.from_array(np.zeros((10,10)))
        vrt.dataset.GetRasterBand(1).SetMetadata({'name':'band1'})
        self.assertEqual(vrt._create_band_name({'name': 'band1'}), ('band1_0000', {}))

    def test_create_band_name_wkv_and_name(self):
        name = 'some_name'
        wkv = dict(short_name='sigma0')
        self.mock_pti['get_wkv_variable'].return_value=wkv
        vrt = VRT()
        self.assertEqual(vrt._create_band_name({'wkv': 'sigma0', 'name': name}), (name, wkv))

    def test_leave_few_bands(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt = VRT.copy_dataset(ds)
        vrt.leave_few_bands([1, 'L_469'])
        self.assertEqual(vrt.dataset.RasterCount,2)
        self.assertEqual(vrt.dataset.GetRasterBand(1).GetMetadataItem(str('name')), 'L_645')
        self.assertEqual(vrt.dataset.GetRasterBand(2).GetMetadataItem(str('name')), 'L_469')

    def test_find_complex_band(self):
        a = np.random.randn(100,100)
        vrt1 = VRT.from_array(a)
        vrt2 = VRT.from_array(a.astype(np.complex64))

        vrt3 = VRT.from_gdal_dataset(vrt1.dataset)
        vrt3.create_bands([{'src': {'SourceFilename': vrt1.filename}},
                           {'src': {'SourceFilename': vrt2.filename}}])

        self.assertEqual(vrt1._find_complex_band(), None)
        self.assertEqual(vrt2._find_complex_band(), 1)
        self.assertEqual(vrt3._find_complex_band(), 2)

    def test_split_complex_bands(self):
        a = np.random.randn(100,100)
        vrt1 = VRT.from_array(a.astype(np.complex64))
        vrt2 = VRT.from_array(a)
        vrt3 = VRT.from_array(a.astype(np.complex64))

        vrt4 = VRT.from_gdal_dataset(vrt1.dataset)
        vrt4.create_bands([{'src': {'SourceFilename': vrt1.filename}, 'dst': {'name': 'vrt1'}},
                           {'src': {'SourceFilename': vrt2.filename}, 'dst': {'name': 'vrt2'}},
                           {'src': {'SourceFilename': vrt3.filename}, 'dst': {'name': 'vrt3'}}])

        vrt4.split_complex_bands()

        self.assertEqual(vrt4.dataset.RasterCount,5)
        self.assertEqual(vrt4.dataset.GetRasterBand(1).GetMetadataItem(str('name')), 'vrt2')
        self.assertEqual(vrt4.dataset.GetRasterBand(2).GetMetadataItem(str('name')), 'vrt1_real')
        self.assertEqual(vrt4.dataset.GetRasterBand(3).GetMetadataItem(str('name')), 'vrt1_imag')
        self.assertEqual(vrt4.dataset.GetRasterBand(4).GetMetadataItem(str('name')), 'vrt3_real')
        self.assertEqual(vrt4.dataset.GetRasterBand(5).GetMetadataItem(str('name')), 'vrt3_imag')

    def test_create_geolocation_bands(self):
        lon, lat = np.meshgrid(np.linspace(0,5,10), np.linspace(10,20,30))
        vrt = VRT.from_lonlat(lon, lat)
        vrt.create_geolocation_bands()

        self.assertEqual(vrt.dataset.RasterCount, 2)
        self.assertEqual(vrt.dataset.GetRasterBand(1).GetMetadataItem(str('name')), 'longitude')
        self.assertEqual(vrt.dataset.GetRasterBand(2).GetMetadataItem(str('name')), 'latitude')
        self.assertTrue(np.allclose(vrt.dataset.GetRasterBand(1).ReadAsArray(), lon))
        self.assertTrue(np.allclose(vrt.dataset.GetRasterBand(2).ReadAsArray(), lat))

    def test_fix_band_metadata(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt = VRT.copy_dataset(ds)
        self.assertIn('standard_name', vrt.dataset.GetRasterBand(1).GetMetadata())
        self.assertIn('time', vrt.dataset.GetRasterBand(1).GetMetadata())
        vrt.fix_band_metadata(['standard_name', 'time'])
        self.assertNotIn('standard_name', vrt.dataset.GetRasterBand(1).GetMetadata())
        self.assertNotIn('time', vrt.dataset.GetRasterBand(1).GetMetadata())

    def test_fix_global_metadata(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt = VRT.copy_dataset(ds)
        vrt.dataset.SetMetadataItem(str('test'), str('"test"'))
        vrt.fix_global_metadata(['AREA_OR_POINT'])
        self.assertNotIn('AREA_OR_POINT', vrt.dataset.GetMetadata())
        self.assertEqual('&quot;test&quot;', vrt.dataset.GetMetadataItem(str('test')))

    def test_hardcopy_bands(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt = VRT.copy_dataset(ds)
        vrt.hardcopy_bands()

        self.assertTrue(np.allclose(vrt.dataset.ReadAsArray(), ds.ReadAsArray()))
        band_nodes = Node.create(str(vrt.xml)).nodeList('VRTRasterBand')
        self.assertEqual(band_nodes[0].node('SourceFilename').value, vrt.band_vrts[1].filename)
        self.assertEqual(band_nodes[1].node('SourceFilename').value, vrt.band_vrts[2].filename)
        self.assertEqual(band_nodes[2].node('SourceFilename').value, vrt.band_vrts[3].filename)

    @patch.multiple(VRT, dataset=DEFAULT, __init__=Mock(return_value=None))
    def test_get_projection_dataset(self, dataset):
        proj = 'SOME_PROJECTION'
        dataset.GetProjection.return_value = proj
        dataset.GetGCPProjection.return_value = ''
        dataset.GetMetadata.return_value = {}

        vrt = VRT()
        proj_src = vrt.get_projection()
        self.assertEqual(proj_src, (proj, 'dataset'))

    @patch.multiple(VRT, dataset=DEFAULT, __init__=Mock(return_value=None))
    def test_get_projection_gcps(self, dataset):
        proj = 'SOME_PROJECTION'
        dataset.GetProjection.return_value = ''
        dataset.GetGCPProjection.return_value = proj
        dataset.GetMetadata.return_value = {}

        vrt = VRT()
        proj_src = vrt.get_projection()
        self.assertEqual(proj_src, (proj, 'gcps'))

    @patch.multiple(VRT, dataset=DEFAULT, __init__=Mock(return_value=None))
    def test_get_projection_geolocation(self, dataset):
        proj = 'SOME_PROJECTION'
        dataset.GetProjection.return_value = ''
        dataset.GetGCPProjection.return_value = ''
        dataset.GetMetadata.return_value = {'SRS': proj}

        vrt = VRT()
        proj_src = vrt.get_projection()
        self.assertEqual(proj_src, (proj, 'geolocation'))

    @patch.multiple(VRT, dataset=DEFAULT, __init__=Mock(return_value=None))
    def test_get_projection_raises_NansatProjectionError(self, dataset):
        dataset.GetProjection.return_value = ''
        dataset.GetGCPProjection.return_value = ''
        dataset.GetMetadata.return_value = {}

        vrt = VRT()
        with self.assertRaises(NansatProjectionError):
            proj = vrt.get_projection()

    def test_repr(self):
        # we mock the entire class [MagicMock(VRT)] and set instance attributes/methods
        mock_vrt1 = MagicMock(VRT, filename='aaa', vrt=None, __repr__=VRT.__repr__)
        mock_vrt2 = MagicMock(VRT, filename='bbb', vrt=mock_vrt1, __repr__=VRT.__repr__)

        # str(mock_vrt1) will call mock_vrt1.__repr__ and, hence, VRT.__repr__ which we test
        self.assertEqual(str(mock_vrt1), 'aaa')
        self.assertEqual(str(mock_vrt2), 'bbb=>aaa')

    @patch.multiple(VRT, create_band=DEFAULT, __init__=Mock(return_value=None))
    def test_add_swath_mask_band(self, create_band):
        vrt = VRT()
        vrt.filename = '/temp/filename.vrt'
        vrt._add_swath_mask_band()
        src = [{'SourceFilename': '/temp/filename.vrt',
                'SourceBand':  1,
                'DataType': 1}]
        dst ={'dataType': 1,
                'wkv': 'swath_binary_mask',
                'PixelFunctionType': 'OnesPixelFunc'}
        create_band.assert_called_once_with(src=src, dst=dst)

    def test_remove_strings_in_metadata_keys(self):
        gdal_metadata = {'aaa': 'bbb', 'NC_GLOBAL#ccc': 'ddd', 'NANSAT_eee': 'fff'}
        rm_strings = ['NC_GLOBAL#', 'NANSAT_']
        new_metadata = VRT._remove_strings_in_metadata_keys(gdal_metadata, rm_strings)
        self.assertEqual(new_metadata, {'aaa': 'bbb', 'ccc': 'ddd', 'eee': 'fff'})

    def test_super_vrt_of_geolocation_bands(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt1 = VRT.from_lonlat(lon, lat)
        vrt1.create_geolocation_bands()
        vrt2 = vrt1.get_super_vrt()
        self.assertTrue(hasattr(vrt2.vrt, 'vrt'))

    def test_get_shifted_vrt(self):
        deg = -10
        vrt1 = VRT.from_array(np.zeros((180,360)))
        vrt1.dataset.SetProjection(NSR().wkt)
        vrt2 = vrt1.get_shifted_vrt(deg)
        self.assertEqual(vrt1.dataset.GetGeoTransform()[0]+deg, vrt2.dataset.GetGeoTransform()[0])

    def test_get_super_vrt(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt1 = VRT.from_gdal_dataset(ds, metadata=ds.GetMetadata())
        vrt2 = vrt1.get_super_vrt()
        self.assertIsInstance(vrt2.vrt, VRT)
        self.assertEqual(vrt2.dataset.GetMetadataItem(str('AREA_OR_POINT')), 'Area')

    def test_get_super_vrt_geolocation(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt1 = VRT.from_lonlat(lon, lat)
        vrt2 = vrt1.get_super_vrt()
        vrt1 = None
        self.assertTrue(vrt2.geolocation.x_vrt is not None)
        self.assertTrue(vrt2.geolocation.y_vrt is not None)

    def test_get_super_vrt_and_copy(self):
        array = np.zeros((10,10))
        vrt = VRT.from_array(array)
        vrt = vrt.get_super_vrt()
        vrt = vrt.copy()
        data = vrt.dataset.ReadAsArray()

        self.assertFalse(data is None)
        self.assertTrue(np.all(data == array))

    def test_get_sub_vrt0(self):
        vrt1 = VRT()
        vrt2 = vrt1.get_sub_vrt()
        self.assertEqual(vrt1, vrt2)

    def test_get_sub_vrt3(self):
        vrt1 = VRT().get_super_vrt().get_super_vrt().get_super_vrt()
        vrt2 = vrt1.get_sub_vrt(3)
        self.assertEqual(vrt2.vrt, None)

    def test_get_sub_vrt_steps_0(self):
        vrt1 = VRT().get_super_vrt()
        vrt2 = vrt1.get_sub_vrt(steps=0)
        self.assertEqual(vrt1, vrt2)

    def test_transform_points(self):
        ds = gdal.Open(self.test_file_gcps)
        vrt1 = VRT.from_gdal_dataset(ds, metadata=ds.GetMetadata())
        vrt1.tps = True
        lon, lat = vrt1.transform_points([1, 2, 3], [4, 5, 6])
        self.assertTrue(np.allclose(lon, np.array([28.23549571, 28.24337106, 28.25126129])))
        self.assertTrue(np.allclose(lat, np.array([71.52509848, 71.51913744, 71.51317568])))
        lon, lat = vrt1.transform_points([], [])
        self.assertTrue(np.allclose(lon, np.array([])))
        self.assertTrue(np.allclose(lat, np.array([])))

    def test_make_filename(self):
        filename1 = VRT._make_filename()
        filename2 = VRT._make_filename(extention='smth')
        filename3 = VRT._make_filename(nomem=True)
        self.assertTrue(filename1.startswith('/vsimem/'))
        self.assertTrue(filename2.startswith('/vsimem/'))
        self.assertTrue(filename2.endswith('.smth'))
        self.assertTrue(os.path.exists(filename3))

    def test_transform_coordinates_list(self):
        src_srs = NSR()
        dst_srs = NSR(str('+proj=stere'))
        src_points = ([1,2,3,4],[5,6,7,8])
        dst_x, dst_y, dst_z = VRT.transform_coordinates(src_srs, src_points, dst_srs)
        # check if shape of the result matches the expected shape (list with four points)
        self.assertEqual(dst_x.shape, (4,))
        self.assertEqual(dst_y.shape, (4,))
        self.assertEqual(dst_z.shape, (4,))

    def test_transform_coordinates_1d_array(self):
        src_srs = NSR()
        dst_srs = NSR(str('+proj=stere'))
        src_points = (np.array([1,2,3,4]), np.array([5,6,7,8]), np.array([5,6,7,8]))
        dst_x, dst_y, dst_z = VRT.transform_coordinates(src_srs, src_points, dst_srs)
        # check if shape of the result matches the expected shape (list with four points)
        self.assertEqual(dst_x.shape, (4,))
        self.assertEqual(dst_y.shape, (4,))
        self.assertEqual(dst_z.shape, (4,))

    def test_transform_coordinates_2d_array(self):
        src_srs = NSR()
        dst_srs = NSR(str('+proj=stere'))
        src_points = (np.array([[1,2,3,4],[1,2,3,4]]),
                      np.array([[5,6,7,8],[5,6,7,8]]),
                      np.array([[5,6,7,8],[5,6,7,8]]),)
        dst_x, dst_y, dst_z = VRT.transform_coordinates(src_srs, src_points, dst_srs)
        # check if shape of the result matches the expected shape (2x4 array)
        self.assertEqual(dst_x.shape, (2,4))
        self.assertEqual(dst_y.shape, (2,4))
        self.assertEqual(dst_z.shape, (2,4))

    def test_reproject_gcps(self):
        lon, lat = np.meshgrid(np.linspace(0, 5, 10), np.linspace(10, 20, 30))
        vrt1 = VRT.from_lonlat(lon, lat)
        vrt1.reproject_gcps(str('+proj=cea'))
        self.assertIn('Cylindrical_Equal_Area', vrt1.dataset.GetGCPProjection())
        self.assertEqual(round(vrt1.dataset.GetGCPs()[0].GCPX), 0)
        self.assertEqual(round(vrt1.dataset.GetGCPs()[0].GCPY), 1100286)
        self.assertEqual(round(vrt1.dataset.GetGCPs()[-1].GCPX), 556597)
        self.assertEqual(round(vrt1.dataset.GetGCPs()[-1].GCPY), 2096057)

if __name__ == "__main__":
    unittest.main()
