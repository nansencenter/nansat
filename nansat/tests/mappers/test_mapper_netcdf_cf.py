import os
import unittest
import datetime
import tempfile

import numpy as np
import gdal
from netCDF4 import Dataset

from nansat.mappers.mapper_netcdf_cf import Mapper

from mock import patch, Mock, DEFAULT

class NetCDF_CF_Tests(unittest.TestCase):

    def setUp(self):
        fd, self.tmp_filename = tempfile.mkstemp(suffix='.nc')
        ds = Dataset(self.tmp_filename, 'w')
        lat_sz = 30
        lon_sz = 20
        height_sz = 10

        # Set dimensions
        ds.createDimension('latitude', lat_sz)
        ds.createDimension('longitude', lon_sz)
        ds.createDimension('time', 3)
        ds.createDimension('pressure', 7)
        ds.createDimension('height', height_sz)
        ds.createDimension('dimension_rgb', 3) # intentionally no variable added..
        # Set variables
        # 1d "dimensional" variables i.e lats, times, etc.
        times = ds.createVariable('time', 'i4', ('time'))
        times.units = 'seconds since 1970-01-01 00:00'
        times.standard_name = 'time'
        times.long_name = 'time'
        times[0] = np.ma.MaskedArray(data = 1560610800.0, mask =
                False, fill_value = 1e+20)
        times[1] = np.ma.MaskedArray(data = 1560621600.0, mask =
                False, fill_value = 1e+20)
        times[2] = np.ma.MaskedArray(data = 1560632400.0, mask =
                False, fill_value = 1e+20)

        heights = ds.createVariable('height', 'i4', ('height'))
        heights[:] = np.linspace(10, 100, height_sz)

        lats = ds.createVariable('latitude', 'i4', ('latitude'))
        lats[:] = np.linspace(0, 60, lat_sz)

        lons = ds.createVariable('longitude', 'i4', ('longitude'))
        lons[:] = np.linspace(0, 20, lon_sz)

        # Spatial variables 2d, 3d, and 4d
        var2d = ds.createVariable('var2d', 'i4', ('latitude', 'longitude'))
        var3d = ds.createVariable('var3d', 'i4', ('time', 'latitude', 'longitude'))
        var3d.standard_name = 'x_wind'
        var4d = ds.createVariable('var4d', 'f4', ('time', 'pressure', 'latitude', 'longitude'))
        var4d.standard_name = 'x_wind'
        var5d = ds.createVariable('var5d', 'f4', ('time', 'pressure', 'height', 'latitude', 'longitude'))
        var5d.standard_name = 'x_wind'
        # gdal should read this as several bands of shape (longitude,pressure)=(20, 7)
        buggy_var = ds.createVariable('buggy_var', 'f4', ('time', 'latitude', 'longitude', 'pressure'))
        buggy_var.standard_name = 'x_wind'
        # A variable with a dimension that is not itself added as a variable
        rgb_var = ds.createVariable('rgb_var', 'f4', ('dimension_rgb', 'latitude', 'longitude'))

        pressures = ds.createVariable('pressure', 'i4', ('pressure'))
        pressures.standard_name = 'air_pressure'
        pressures.description = 'pressure'
        pressures.long_name = 'pressure'
        pressures.positive = 'down'
        pressures.units = 'hPa'
        pressures[0] = np.ma.MaskedArray(data = 200., mask = False, fill_value = 1e+20)
        pressures[1] = np.ma.MaskedArray(data = 250., mask = False, fill_value = 1e+20)
        pressures[2] = np.ma.MaskedArray(data = 300., mask = False, fill_value = 1e+20)
        pressures[3] = np.ma.MaskedArray(data = 400., mask = False, fill_value = 1e+20)
        pressures[4] = np.ma.MaskedArray(data = 500., mask = False, fill_value = 1e+20)
        pressures[5] = np.ma.MaskedArray(data = 700., mask = False, fill_value = 1e+20)
        pressures[6] = np.ma.MaskedArray(data = 800., mask = False, fill_value = 1e+20)
        ds.close()
        os.close(fd) # Just in case - see https://www.logilab.org/blogentry/17873

        fd, self.tmp_filename_xy = tempfile.mkstemp(suffix='.nc')
        ds = Dataset(self.tmp_filename_xy, 'w')
        y_sz = 30
        x_sz = 20

        # Set dimensions
        ds.createDimension('y', y_sz)
        ds.createDimension('x', x_sz)
        ds.createDimension('some_times', 3)
        # Set variables
        # 1d "dimensional" variables i.e lats, times, etc.
        times = ds.createVariable('some_times', 'i4', ('some_times'))
        times.units = 'seconds since 1970-01-01 00:00'
        times.standard_name = 'time'
        times.long_name = 'time'
        times[0] = np.ma.MaskedArray(data = 1560610800.0, mask =
                False, fill_value = 1e+20)
        times[1] = np.ma.MaskedArray(data = 1560621600.0, mask =
                False, fill_value = 1e+20)
        times[2] = np.ma.MaskedArray(data = 1560632400.0, mask =
                False, fill_value = 1e+20)

        ys = ds.createVariable('y', 'i4', ('y'))
        ys[:] = np.linspace(0, 60, y_sz)

        xs = ds.createVariable('x', 'i4', ('x'))
        xs[:] = np.linspace(0, 20, x_sz)

        # Spatial variables 2d and 3d
        var2d = ds.createVariable('var2d', 'i4', ('y', 'x'))
        var3d = ds.createVariable('var3d', 'i4', ('some_times', 'y', 'x'))
        var3d.standard_name = 'x_wind'

        ds.close()
        os.close(fd) # Just in case - see https://www.logilab.org/blogentry/17873

        fd, self.tmp_filename_no_time_var = tempfile.mkstemp(suffix='.nc')
        ds = Dataset(self.tmp_filename_no_time_var, 'w')
        y_sz = 30
        x_sz = 20

        # Set dimensions
        ds.createDimension('y', y_sz)
        ds.createDimension('x', x_sz)
        # Set variables
        # 1d "dimensional" variables i.e lats, times, etc.
        ys = ds.createVariable('y', 'i4', ('y'))
        ys.standard_name = 'projection_y_coordinate'
        ys[:] = np.linspace(0, 60, y_sz)

        xs = ds.createVariable('x', 'i4', ('x'))
        xs.standard_name = 'projection_x_coordinate'
        xs[:] = np.linspace(0, 20, x_sz)

        # Spatial variables 2d and 3d
        var2d = ds.createVariable('var2d', 'i4', ('y', 'x'))
        var2d.standard_name = 'x_wind'

        ds.close()
        os.close(fd) # Just in case - see https://www.logilab.org/blogentry/17873

    def tearDown(self):
        os.unlink(self.tmp_filename)
        os.unlink(self.tmp_filename_xy)
        os.unlink(self.tmp_filename_no_time_var)

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test__timevarname(self, mock_init):
        mock_init.return_value = None
        mm = Mapper()
        mm.input_filename = self.tmp_filename_no_time_var
        timevar_name = mm._timevarname()
        self.assertEqual(timevar_name, '')

        mm = Mapper()
        mm.input_filename = self.tmp_filename_xy
        timevar_name = mm._timevarname()
        self.assertEqual(timevar_name, 'some_times')

        mm = Mapper()
        mm.input_filename = self.tmp_filename
        timevar_name = mm._timevarname()
        self.assertEqual(timevar_name, 'time')
        

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper._time_reference')
    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test_time_count_to_np_datetime64(self, mock_init, mock_units):
        mock_init.return_value = None
        tu = (datetime.datetime(1900, 1, 1, 0, 0), 'days since 1900-1-1 0:0:0 +0')
        mock_units.return_value = tu
        mm = Mapper()
        time_count = '43648.22734953704'

        # TEST DAYS
        tt = mm._time_count_to_np_datetime64(time_count)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)

        tt = mm._time_count_to_np_datetime64(time_count, time_reference=tu)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)
        

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper._time_reference')
    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test_time_count_to_np_datetime64(self, mock_init, mock_units):
        mock_init.return_value = None
        tu = (datetime.datetime(1900, 1, 1, 0, 0), 'days since 1900-1-1 0:0:0 +0')
        mock_units.return_value = tu
        mm = Mapper()
        time_count = '43648.22734953704'

        # TEST DAYS
        tt = mm._time_count_to_np_datetime64(time_count)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)

        tt = mm._time_count_to_np_datetime64(time_count, time_reference=tu)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)

        time_count = '43648'
        tt = mm._time_count_to_np_datetime64(time_count)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)
        self.assertEqual(tt, np.datetime64('2019-07-04'))

        # TEST HOURS
        tu = (datetime.datetime(1900, 1, 1, 0, 0), 'hours since 1900-1-1 0:0:0 +0')
        tt = mm._time_count_to_np_datetime64(time_count, time_reference=tu)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)
        self.assertEqual(tt, np.datetime64('1904-12-24T16:00:00.000000'))

        time_count = '43648.22734953704'
        tt = mm._time_count_to_np_datetime64(time_count, time_reference=tu)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)
        self.assertEqual(tt, np.datetime64('1904-12-24T16:13:38.458333'))

        # TEST MINUTES
        tu = (datetime.datetime(1900, 1, 1, 0, 0), 'minutes since 1900-1-1 0:0:0 +0')
        tt = mm._time_count_to_np_datetime64(time_count, time_reference=tu)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)
        self.assertEqual(tt, np.datetime64('1900-01-31T07:28:13.640972'))

        # TEST SECONDS
        tu = (datetime.datetime(1900, 1, 1, 0, 0), 'seconds since 1900-1-1 0:0:0 +0')
        tt = mm._time_count_to_np_datetime64(time_count, time_reference=tu)
        # Assert data type of tt is np.datetime64
        self.assertEqual(type(tt), np.datetime64)
        self.assertEqual(tt, np.datetime64('1900-01-01T12:07:28.227350'))

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test_variable_with_a_dimension_that_is_not_itself_added_as_a_variable(self, mock_init):
        mock_init.return_value = None
        mm = Mapper()
        mm.input_filename = self.tmp_filename
        fn = 'NETCDF:"' + self.tmp_filename + '":rgb_var'
        bdict = mm._get_band_from_subfile(fn)
        self.assertEqual(bdict['src']['SourceBand'], 1)

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test_with_xy_dims(self, mock_init):
        mock_init.return_value = None
        mm = Mapper()
        mm.input_filename = self.tmp_filename_xy
        fn = 'NETCDF:"' + self.tmp_filename + '":var3d'
        bdict = mm._get_band_from_subfile(fn, 
                bands=['x_wind'])
        self.assertEqual(bdict['src']['SourceBand'], 1)

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test_buggy_var(self, mock_init):
        """ The last band dimensions should be latitude and longitude - otherwise gdal will fail in
        reading the data correctly.
        
        This is to confirm that this understanding is correct..

        The shape of buggy_var is ('time', 'latitude', 'longitude', 'pressure')
        """
        mock_init.return_value = None
        mm = Mapper()
        mm.input_filename = self.tmp_filename
        fn = 'NETCDF:"' + self.tmp_filename + '":buggy_var'
        bdict = mm._get_band_from_subfile(fn, 
                bands=['x_wind'])
        self.assertEqual(bdict['src']['SourceBand'], 1)
        self.assertEqual(bdict['dst']['NETCDF_DIM_latitude'], '0')
        self.assertEqual(bdict['dst']['time_iso_8601'], np.datetime64('2019-06-15T15:00:00.000000'))
        subds = gdal.Open(fn)
        self.assertEqual(subds.RasterXSize, 7) # size of pressure dimension
        self.assertEqual(subds.RasterYSize, 20) # size of longitude dimension

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test__get_band_from_subfile__var4d(self, mock_init):
        mock_init.return_value = None
        mm = Mapper()
        mm.input_filename = self.tmp_filename
        fn = 'NETCDF:"' + self.tmp_filename + '":var4d'
        bdict200 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'time': np.datetime64('2019-06-15T18:00'), # 2nd band
                    'pressure': 200}, 
                bands=['x_wind'])
        bdict500 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'time': np.datetime64('2019-06-15T18:00'), # 2nd band
                    'pressure': 500}, 
                bands=['x_wind'])

        self.assertEqual(bdict200['src']['SourceBand'], 8)
        self.assertEqual(bdict200['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict200['dst']['time_iso_8601'], np.datetime64('2019-06-15T18:00:00.000000'))
        self.assertEqual(bdict500['src']['SourceBand'], 12)
        self.assertEqual(bdict500['dst']['NETCDF_DIM_pressure'], '500')
        self.assertEqual(bdict500['dst']['time_iso_8601'], np.datetime64('2019-06-15T18:00:00.000000'))

    @patch('nansat.mappers.mapper_netcdf_cf.Mapper.__init__')
    def test__get_band_from_subfile__var5d(self, mock_init):
        mock_init.return_value = None
        mm = Mapper()
        mm.input_filename = self.tmp_filename
        fn = 'NETCDF:"' + self.tmp_filename + '":var5d'
        # should give 1 band when time, pressure and height is in netcdf_dim
        bdict1 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'time': np.datetime64('2019-06-15T18:00'),
                    'pressure': 200,
                    'height': 20}, 
                bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '20')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T18:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 72)
        # should give first height_sz band when only time and pressure is in netcdf_dim
        bdict1 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'time': np.datetime64('2019-06-15T18:00'),
                    'pressure': 200},
                bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '10')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T18:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 71)
        # should give first height_sz and pressure bands when only time is in netcdf_dim
        bdict1 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'time': np.datetime64('2019-06-15T18:00'),
                },
                bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '10')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T18:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 71)

        bdict1 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'pressure': 300, # 3rd band
                },
                bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '10')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '300')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T15:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 21)

        bdict1 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'height': 30, # 3rd band
                },
                bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '30')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T15:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 3)

        bdict1 = mm._get_band_from_subfile(fn, 
                netcdf_dim = {
                    'time': np.datetime64('2019-06-15T15:00'),
                },
                bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '10')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T15:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 1)

        # should give first band when netcdf_dim is empty
        bdict1 = mm._get_band_from_subfile(fn, bands=['x_wind'])
        self.assertEqual(bdict1['dst']['NETCDF_DIM_height'], '10')
        self.assertEqual(bdict1['dst']['NETCDF_DIM_pressure'], '200')
        self.assertEqual(bdict1['dst']['time_iso_8601'], np.datetime64('2019-06-15T15:00:00.000000'))
        self.assertEqual(bdict1['src']['SourceBand'], 1)
