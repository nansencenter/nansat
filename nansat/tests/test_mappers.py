import unittest 

import os, datetime
import numpy as np

from nansat.nansat import Nansat

import nansat_test_data as ntd

class NetCDFCFMapperTests(unittest.TestCase):

    def setUp(self):
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.test_file_arome_metcoop = '/vagrant/shared/test_data/generic/arome_metcoop_default2_5km_20170227T00Z.nc'
        self.test_file_arome_arctic = '/vagrant/shared/test_data/generic/arome_arctic_pp_2_5km_20170227T00Z.nc'
        self.test_file_ecmwf = '/vagrant/shared/test_data/generic/ec_atmo_0_1deg_20170227T000000Z_1h.nc'

    def test_netcdfcf_mapper_is_used(self):
        n = Nansat(self.test_file_arctic, mapperName='netcdfcf')
        self.assertEqual(n.mapper, 'netcdfcf')

    def test_open_netcdf_cf(self):
        n = Nansat(self.test_file_arctic, mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)

    def test_open_arome_metcoop(self):
        n = Nansat(self.test_file_arome_metcoop, mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)

    def test_open_arome_metcoop_at_given_time(self):
        n = Nansat(self.test_file_arome_metcoop, netcdf_dim={'time': '1488153600'},
                mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)

    def test_open_arome_metcoop_at_given_height(self):
        n = Nansat(self.test_file_arome_metcoop, netcdf_dim={'height0': '0'},
                mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)

    def test_open_arome_metcoop_y_wind_and_x_wind_at_given_time(self):
        n = Nansat(self.test_file_arome_metcoop, bands=['y_wind', 'x_wind'],
                netcdf_dim={'time': '1488153600'}, mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))

    def test_open_arome_arctic_y_wind_and_x_wind_at_given_time(self):
        n = Nansat(self.test_file_arome_arctic, bands=['y_wind', 'x_wind'],
                netcdf_dim={'time': '1488387600'}, mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))

    def test_open_arome_arctic_y_wind_and_x_wind_at_given_datetime(self):
        n = Nansat(self.test_file_arome_arctic, bands=['y_wind', 'x_wind'],
            netcdf_dim={'time':
                np.datetime64(datetime.datetime(2017,2,28,15,30,0))},
                mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))

    def test_open_ecmwf_y_wind_and_x_wind_at_given_time(self):
        n = Nansat(self.test_file_ecmwf, bands=['y_wind', 'x_wind'],
                netcdf_dim={'time': '1488409200'}, mapperName='netcdfcf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))

if __name__ == "__main__":
    unittest.main()
