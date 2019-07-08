# coding=utf-8
#------------------------------------------------------------------------------
# Name:         test_mapper_opendap.py
# Purpose:      Test the Opendap class
#
# Author:       Artem Moiseev
#
# Created:      15.08.2018
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
import sys
import unittest
from netCDF4 import Dataset
import tempfile
from mock import patch
from nansat.exceptions import WrongMapperError
from collections import OrderedDict
from nansat.mappers.opendap import Opendap
import numpy as np
import warnings


class OpenDAPTests(unittest.TestCase):

    def setUp(self):
        _, tmp_filename = tempfile.mkstemp(suffix='.nc')
        ds = Dataset(tmp_filename, 'w')
        lat_sz = 30
        lon_sz = 20
        values = np.random.random_sample((lat_sz, lon_sz))

        # Set dimensions
        ds.createDimension('lat', lat_sz)
        ds.createDimension('lon', lon_sz)
        ds.createDimension('time', 3)
        ds.createDimension('depth', 10)
        # Set variables
        # 1d "dimensional" variables i.e lats, times, etc.
        times = ds.createVariable('var2', 'i4', ('time'))
        lats = ds.createVariable('lat', 'i4', ('lat'))
        lats[:] = np.linspace(0, 60, lat_sz)
        lons = ds.createVariable('lon', 'i4', ('lon'))
        lons[:] = np.linspace(0, 20, lon_sz)
        # Spatial variables 2d, 3d, and 4d
        ds.createVariable('var2d', 'i4', ('lat', 'lon'))
        ds.createVariable('var3d', 'i4', ('time', 'lat', 'lon'))
        ds.createVariable('var4d', 'f4', ('time', 'depth', 'lat', 'lon'))
        # Initiate the Opendap object
        self.od = Opendap()
        self.od.baseURLs = ['http://first.no', 'http://second.com']
        self.od.xName = 'lon'
        self.od.yName = 'lat'
        self.od.timeVarName = 'time'
        self.od.ds = ds
        self.ds = ds

    def test_test_mapper(self):
        res_ok = self.od.test_mapper(filename='http://first.no/path/to/the/file.nc')
        self.assertIsNone(res_ok)
        with self.assertRaises(WrongMapperError):
            self.od.test_mapper(filename='http://not-in-base-urls.net/path/to/the/file.nc')

    """ Not sure how to test get_dataset.. 
    
    According toe the followingm, mocking C modules seems not to be possible.. Can that be the
    reason the below doesn't work?

    See https://stackoverflow.com/questions/192649/can-you-monkey-patch-methods-on-core-types-in-python/192857#192857
    """
    ##@patch('nansat.mappers.opendap.Dataset.__init__')
    #@patch.object('nansat.mappers.opendap.Dataset', '__init__', return_value = None)
    #def test_get_dataset(self, mock_dataset):
    #    dd = Dataset()
    #    self.od.filename = 'http://something.that.is/mocked'
    #    ds1 = self.od.get_dataset(None)
    #    self.assertIsInstance(ds1, Dataset)

    def test_get_dataset_which_does_not_exist(self):
        wrong_filename = '/path/which/does/not/exist/file.nc'
        self.od.filename = wrong_filename
        if sys.version_info.major==2:
            with self.assertRaises(IOError) as ve:
                self.od.get_dataset(None)
        else:
            with self.assertRaises(FileNotFoundError) as ve:
                self.od.get_dataset(None)
        self.assertEqual(ve.exception.errno, 2)

        with self.assertRaises(ValueError) as ve:
            self.od.get_dataset([])
        self.assertEqual(ve.exception.args[0], 'Input ds is not netCDF.Dataset!')

    def test_get_geospatial_variable_names(self):
        ds_vars = self.od.get_geospatial_variable_names()
        self.assertEqual(len(ds_vars), 3)
        self.assertIn('var2d', ds_vars)
        self.assertIn('var3d', ds_vars)
        self.assertIn('var4d', ds_vars)

    def test_get_dataset_time(self):
        pass

    def test_get_layer_datetime(self):
        date1 = '2010-01-02'
        datetimes_1 = np.arange(np.datetime64('2010-01-01'), np.datetime64('2010-01-02'))
        res_layer_num1, res_layer_date1 = Opendap.get_layer_datetime(date1, datetimes_1)

        self.assertEqual(res_layer_num1, 0)
        self.assertIsInstance(res_layer_date1, np.datetime64)

        datetimes_2 = np.arange(np.datetime64('2009-01-01'), np.datetime64('2009-01-08'))

        with self.assertRaises(ValueError):
            Opendap.get_layer_datetime(date1, datetimes_2)

        datetimes_3 = np.arange(np.datetime64('2009-12-31'), np.datetime64('2010-01-08'))
        res_layer_num3, res_layer_date3 = Opendap.get_layer_datetime(date1, datetimes_3)
        self.assertEqual(res_layer_date3, np.datetime64(date1))
        self.assertEqual(res_layer_num3, 2)

    @patch('nansat.mappers.opendap.gdal.Open')
    def test_metaitem(self, mock_open):
        mock_open.return_value = None
        self.od.ds.variables['var3d'].setncattr('test_attr', 'test_val')
        res1 = self.od.get_metaitem('https://file-url.nc', 'var3d', (0, 'y', 'x'))
        self.assertIsInstance(res1, dict)
        self.assertEqual(len(res1.keys()), 2)
        self.assertIn('dst', res1.keys())
        self.assertIn('src', res1.keys())
        self.assertIn('SourceFilename', res1['src'].keys())
        self.assertIn('SourceBand', res1['src'].keys())
        self.assertIn('name', res1['dst'].keys())
        self.assertIn('dataType', res1['dst'].keys())
        self.assertIn('test_attr', res1['dst'].keys())
        self.od.ds.variables['var3d'].delncattr('test_attr')

    @patch('nansat.mappers.opendap.gdal.Open')
    def test_get_metaitem_spec_attrs(self, mock_open):
        mock_open.return_value = None
        test_cases = [
            {'key': 'offset', 'val': 'offset_value', 'meta_key': 'ScaleOffset'},
            {'key': 'add_offset', 'val': 'add_offset_value', 'meta_key': 'ScaleOffset'},
            {'key': 'scale', 'val': 'scale_value', 'meta_key': 'ScaleRatio'},
            {'key': 'scale_factor', 'val': 'scale_factor_value', 'meta_key': 'ScaleRatio'},
        ]

        for test_case in test_cases:
            self.od.ds.variables['var3d'].setncattr(test_case['key'], test_case['val'])
            res1 = self.od.get_metaitem('https://file-url.nc', 'var3d', (0, 'y', 'x'))
            self.assertIn((test_case['meta_key'], test_case['val']), res1['src'].items())
            self.od.ds.variables['var3d'].delncattr(test_case['key'])

    def test_fix_encoding(self):
        self.assertEqual(Opendap._fix_encoding(u'åsnes'), 'snes')
        self.assertEqual(Opendap._fix_encoding('asnes'), 'asnes')

    def test_filter_dimensions(self):
        res1 = list(filter(self.od._filter_dimensions, ['time', 'var1', 'lat', 'lon']))
        self.assertIsInstance(res1, list)
        self.assertEqual(len(res1), 1)
        self.assertEqual(res1[0], 'var1')
        res2 = list(filter(self.od._filter_dimensions, ['time', 'lat', 'lon']))
        self.assertEqual(len(res2), 0)

    @patch('nansat.mappers.opendap.gdal.Open')
    def test_create_metadict(self, mock_open):
        mock_open.return_value = None
        res1 = self.od.create_metadict('test.nc', ['var4d'], 0)
        self.assertIsInstance(res1, list)
        self.assertEqual(len(res1), 10)

        for i in range(len(res1)):
            source_filename = res1[i]['src']['SourceFilename'][19:31]
            self.assertEqual(source_filename, '[0][%s][y][x]' % i)

        res2 = self.od.create_metadict('test.nc', ['var2d', 'var3d'], 0)
        self.assertEqual(len(res2), 2)
        self.assertEqual(res2[0]['src']['SourceFilename'][19:25], '[y][x]')
        self.assertEqual(res2[1]['src']['SourceFilename'][19:28], '[0][y][x]')

    def test_get_time_coverage_resolution(self):
        res1 = self.od.get_time_coverage_resolution()
        self.assertIsInstance(res1, int)
        self.assertEqual(res1, 0)
        self.od.ds.setncattr('time_coverage_resolution', 'P2D')
        res2 = self.od.get_time_coverage_resolution()
        self.assertEqual(res2, 172800)
        self.od.ds.setncattr('time_coverage_resolution', 'wrong_value')
        with warnings.catch_warnings():
            _ = self.od.get_time_coverage_resolution()
        self.od.ds.delncattr('time_coverage_resolution')

    def test_get_shape(self):
        res = self.od.get_shape()
        self.assertIsInstance(res, tuple)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0], 20)
        self.assertEqual(res[1], 30)

    def test_get_geotransfort(self):
        res = self.od.get_geotransform()
        self.assertIsInstance(res, tuple)
        self.assertEqual(len(res), 6)
        self.assertEqual(res, (0, 1, 0, 0, 0, 2))
