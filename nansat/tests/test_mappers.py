import unittest 

import os

from nansat.nansat import Nansat

import nansat_test_data as ntd

class MapperTests(unittest.TestCase):

    def setUp(self):
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.test_file_arome = '/vagrant/shared/test_data/generic/arome_metcoop_default2_5km_20170227T00Z.nc'
        self.tmpfilename = os.path.join(ntd.tmp_data_path, 'test.nc')

    def test_open_netcdf_cf(self):
        n = Nansat(self.test_file_arctic, mapperName='netcdfCF')

    def test_open_arome_netcdf_cf_file(self):
        n = Nansat(self.test_file_arome, mapperName='netcdfcf')

if __name__ == "__main__":
    unittest.main()
