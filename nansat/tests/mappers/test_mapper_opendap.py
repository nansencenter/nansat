import unittest
from netCDF4 import Dataset
import tempfile
from mock import patch
from nansat.exceptions import WrongMapperError
from collections import OrderedDict
from nansat.mappers.opendap import Opendap
import numpy as np


class OpenDAPTests(unittest.TestCase):

    def setUp(self):
        ds_file = tempfile.NamedTemporaryFile(suffix='.nc', mode='w')
        ds = Dataset(ds_file.name, 'w')
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
        var2d = ds.createVariable('var2d', 'i4', ('lat', 'lon'))
        var2d.setncattr('scale', 'scale_value')
        var2d.setncattr('scale_factor', 'scale_factor_value')
        var2d.setncattr('useless_attr', 'useless_value')
        var3d = ds.createVariable('var3d', 'i4', ('time', 'lat', 'lon'))
        var3d.setncattr('offset', 'offset_value')
        var3d.setncattr('add_offset', 'add_offset_value')
        var4d = ds.createVariable('var4d', 'f4', ('time', 'depth', 'lat', 'lon'))
        self.ds = ds

    def test_test_mapper(self):
        od = Opendap()
        od.baseURLs = ['http://first.no', 'http://second.com']
        res_ok = od.test_mapper(filename='http://first.no/path/to/the/file.nc')
        self.assertIsNone(res_ok)
        with self.assertRaises(WrongMapperError):
            od.test_mapper(filename='http://not-in-base-urls.net/path/to/the/file.nc')

    def test_get_dataset(self):
        od = Opendap()
        od.filename = 'http://www.ifremer.fr/opendap/cerdap1/globcurrent/v2.0/global_025_deg/' \
                      'total_hs/2010/002/20100102000000-GLOBCURRENT-L4-CUReul_hs-ALT_' \
                      'SUM-v02.0-fv01.0.nc'
        ds1 = od.get_dataset(None)
        self.assertIsInstance(ds1, Dataset)

        wrong_filename = '/path/which/does/not/exist/file.nc'
        od.filename = wrong_filename
        with self.assertRaises(ValueError) as ve:
            ds1 = od.get_dataset(None)
            self.assertEqual(ve.args[0], 'Cannot open /path/which/does/not/exist/file.nc')

        with self.assertRaises(ValueError) as ve:
            ds1 = od.get_dataset([])
            self.assertEqual(ve.args[0], 'Input ds is not netCDF.Dataset!')

    def test_get_geospatial_variable_names(self):
        od = Opendap()
        od.xName = 'lon'
        od.yName = 'lat'
        od.ds = self.ds
        ds_vars = od.get_geospatial_variable_names()
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

    def test_get_metaitem(self):
        url = 'https://file-url.nc'
        dims1 = (0, 'y', 'x')
        od = Opendap()
        od.ds = self.ds
        res1 = od.get_metaitem(url, 'var3d', dims1)

        self.assertIn((b'offset', "b'offset_value'"), res1['dst'].items())
        self.assertIn((b'add_offset', "b'add_offset_value'"), res1['dst'].items())
        self.assertIsInstance(res1, dict)
        self.assertEqual(len(res1.keys()), 2)
        self.assertIn('dst', res1.keys())
        self.assertIn('src', res1.keys())
        self.assertEqual(res1['src']['SourceFilename'],
                         'https://file-url.nc?var3d.var3d[%s][%s][%s]' %
                         (dims1[0], dims1[1], dims1[2]))

        res2 = od.get_metaitem(url, 'var2d', dims1)
        self.assertIn((b'scale', "b'scale_value'"), res2['dst'].items())
        self.assertIn((b'scale_factor', "b'scale_factor_value'"), res2['dst'].items())
        self.assertIn((b'useless_attr', "b'useless_value'"), res2['dst'].items())

    def test_create_vrt(self):
        pass

    def test_filter_dimensions(self):
        od = Opendap()
        od.timeVarName = 'time'
        od.yName = 'y'
        od.xName = 'x'
        res1 = list(filter(od._filter_dimensions, ['time', 'var1', 'y', 'x']))
        self.assertIsInstance(res1, list)
        self.assertEqual(len(res1), 1)
        self.assertEqual(res1[0], 'var1')
        res2 = list(filter(od._filter_dimensions, ['time', 'y', 'x']))
        self.assertEqual(len(res2), 0)

    def test_create_metadict(self):
        od = Opendap()
        od.ds = self.ds
        od.timeVarName = 'time'
        od.yName = 'lat'
        od.xName = 'lon'
        res1 = od.create_metadict('test.nc', ['var4d'], 0)
        self.assertIsInstance(res1, list)
        self.assertEqual(len(res1), 10)

        for i in range(len(res1)):
            source_filename = res1[i]['src']['SourceFilename'][19:31]
            self.assertEqual(source_filename, '[0][%s][y][x]' % i)

        res2 = od.create_metadict('test.nc', ['var2d', 'var3d'], 0)
        self.assertEqual(len(res2), 2)
        self.assertEqual(res2[0]['src']['SourceFilename'][19:25], '[y][x]')
        self.assertEqual(res2[1]['src']['SourceFilename'][19:28], '[0][y][x]')

    def test_get_time_coverage_resolution(self):
        od = Opendap()
        od.ds = self.ds
        res = od.get_time_coverage_resolution()
        self.assertIsInstance(res, int)
        self.assertEqual(res, 0)

    def test_get_shape(self):
        od = Opendap()
        od.ds = self.ds
        od.yName = 'lat'
        od.xName = 'lon'
        res = od.get_shape()
        self.assertIsInstance(res, tuple)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0], 20)
        self.assertEqual(res[1], 30)

    def test_get_geotransfort(self):
        od = Opendap()
        od.ds = self.ds
        od.yName = 'lat'
        od.xName = 'lon'
        res = od.get_geotransform()
        self.assertIsInstance(res, tuple)
        self.assertEqual(len(res), 6)
        self.assertEqual(res, (0, 1, 0, 0, 0, 2))
