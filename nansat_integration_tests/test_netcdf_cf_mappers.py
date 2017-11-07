import unittest 

import os, datetime
import numpy as np

from nansat.nansat import Nansat

from nansat.tests import nansat_test_data as ntd

class NetCDFCFMapperTests(unittest.TestCase):

    def setUp(self):
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.test_file_arome_metcoop = '/vagrant/shared/test_data/generic/arome_metcoop_default2_5km_20170227T00Z.nc'
        self.test_file_arome_arctic = '/vagrant/shared/test_data/generic/arome_arctic_pp_2_5km_20170227T00Z.nc'
        self.test_file_ecmwf = '/vagrant/shared/test_data/generic/ec_atmo_0_1deg_20170227T000000Z_1h.nc'
        self.test_file_arome_opendap = 'http://thredds.met.no/thredds/catalog/arome25/catalog.html?dataset=arome25/arome_metcoop_default2_5km_latest.nc'
        self.s1aEW = '/vagrant/shared/test_data/sentinel1_l1/S1A_EW_GRDM_1SDH_20170227T065537_20170227T065637_015466_019652_DD9C.SAFE'
        self.s1bIW = '/vagrant/shared/test_data/sentinel1_l1/S1B_IW_GRDM_1SDV_20170227T061040_20170227T061113_004482_007CD7_3205.SAFE'

    def test_netcdf_cf_mapper_is_used(self):
        n = Nansat(self.test_file_arctic)
        self.assertEqual(n.mapper, 'netcdf_cf')

    def test_open_netcdf_cf(self):
        n = Nansat(self.test_file_arctic, mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)

    def test_open_arome_metcoop(self):
        n = Nansat(self.test_file_arome_metcoop, mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_open_arome_metcoop_at_given_time(self):
        n = Nansat(self.test_file_arome_metcoop, netcdf_dim={'time': '1488153600'},
                mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_open_arome_metcoop_at_given_height(self):
        n = Nansat(self.test_file_arome_metcoop, netcdf_dim={'height0': '0'},
                mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertTrue(n['surface_air_pressure'].any())

    def test_open_arome_metcoop_y_wind_and_x_wind_at_given_time(self):
        n = Nansat(self.test_file_arome_metcoop, bands=['y_wind', 'x_wind'],
                netcdf_dim={'time': '1488153600'}, mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_open_arome_arctic_y_wind_and_x_wind_at_given_time(self):
        n = Nansat(self.test_file_arome_arctic, bands=['y_wind', 'x_wind'],
                netcdf_dim={'time': '1488387600'}, mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_open_arome_arctic_y_wind_and_x_wind_at_given_datetime(self):
        n = Nansat(self.test_file_arome_arctic, bands=['y_wind', 'x_wind'],
            netcdf_dim={'time':
                np.datetime64(datetime.datetime(2017,2,28,15,30,0))},
                mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_open_ecmwf_y_wind_and_x_wind_at_given_time(self):
        n = Nansat(self.test_file_ecmwf, bands=['y_wind', 'x_wind'],
                netcdf_dim={'time': '1488409200'}, mapperName='netcdf_cf')
        self.assertIsInstance(n, Nansat)
        self.assertEqual(2, len(n.bands()))
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_arome_mapper_is_used(self):
        n = Nansat(self.test_file_arome_arctic)
        self.assertEqual(n.mapper, 'arome')
        n = Nansat(self.test_file_arome_metcoop)
        self.assertEqual(n.mapper, 'arome')
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_ecmwf_mapper_is_used(self):
        n = Nansat(self.test_file_ecmwf)
        self.assertEqual(n.mapper, 'ecmwf_metno')
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_mapper_opendap_arome(self):
        n = Nansat(self.test_file_arome_opendap, mapperName='opendap_arome')
        self.assertEqual(n.mapper, 'opendap_arome')
        self.assertTrue(n['x_wind_10m'].any())
        self.assertTrue(n['y_wind_10m'].any())

    def test_reproject_arome_to_SAR(self):
        sar = Nansat(self.s1bIW)
        wind = Nansat(self.test_file_arome_metcoop, netcdf_dim={'time':
            np.datetime64(sar.time_coverage_start)},
            bands=['y_wind','x_wind'])
        self.assertTrue(wind['x_wind_10m'].any())
        self.assertTrue(wind['y_wind_10m'].any())
        wind.reproject(sar, addmask=False)
        self.assertTrue(wind['x_wind_10m'].any())
        self.assertTrue(wind['y_wind_10m'].any())

    def test_reproject_ecmwf_to_SAR(self):
        sar = Nansat(self.s1aEW)
        wind = Nansat(self.test_file_ecmwf, netcdf_dim={'time':
            np.datetime64(sar.time_coverage_start)},
            bands=['y_wind','x_wind'])
        self.assertTrue(wind['x_wind_10m'].any())
        self.assertTrue(wind['y_wind_10m'].any())
        wind.reproject(sar, addmask=False)
        self.assertTrue(wind['x_wind_10m'].any())
        self.assertTrue(wind['y_wind_10m'].any())


    def test_issue_193(self):
        fn = [
            '/vagrant/shared/test_data/cmems/GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS-x10-X30-y55-Y73-201705181200-201705271200.nc',
            '/vagrant/shared/test_data/cmems/ARC-METNO-ARC-TOPAZ4_2_PHYS-FOR-TDS-x10-X30-y55-Y73-20170518-20170526.nc',
            '/vagrant/shared/test_data/cmems/GLOBAL_ANALYSIS_FORECAST_BIO_001_014-TDS-x-180-X179.5-y-89-Y90-20170520-20170527.nc',
        ]
        for f in fn:
            n = Nansat(f)
            self.assertTrue(n.get_metadata().has_key('time_coverage_start'))
            self.assertTrue(n.get_metadata().has_key('time_coverage_end'))
            self.assertTrue(n.get_metadata().has_key('instrument'))
            self.assertTrue(n.get_metadata().has_key('platform'))
            self.assertEqual(n.mapper, 'cmems')

if __name__ == "__main__":
    unittest.main()
