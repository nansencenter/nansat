#-------------------------------------------------------------------------------
# Name:     test_radarsat2.py
# Purpose:      End to end testing of Nansat for Radarsat-2
#
# Author:       Morten Wergeland Hansen
# Modified: Morten Wergeland Hansen
#
# Created:  17.04.2015
# Last modified:17.04.2015 10:24
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import os
import unittest
import numpy as np

from nansat.nansat import Nansat
from mapper_test_archive import DataForTestingMappers

class TestRadarsat(object):

    def test_all_rs2_files(self):
        testData = DataForTestingMappers()
        testData.download_all_test_data()
        for rsfile in testData.mapperData['radarsat2']:
            # OBS: do not yield functions that have the word 'test' in
            # their names - these are run automatically by nose...
            yield self.incidence_angle, rsfile
            #yield self.export2thredds, rsfile
            yield self.export, rsfile

    def export2thredds(self, rsfile):
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        orig.export2thredds(ncfile, bands = {'incidence_angle': {}})
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        np.testing.assert_allclose(inc0, inc1)
        os.unlink(ncfile)

    def export(self, rsfile):
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        orig.export(ncfile)
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        lon0, lat0 = orig.get_geolocation_grids()
        lon1, lat1 = copy.get_geolocation_grids()
        sigma0_0 = orig['sigma0_HH']
        sigma0_1 = copy['sigma0_HH']
        np.testing.assert_allclose(lon0, lon1)
        np.testing.assert_allclose(lat0, lat1)
        # If the next tests fail, it could indicate that the data is flipped
        # check by pyplot.imshow orig vs copy...
        np.testing.assert_allclose(inc0, inc1, rtol=1e-3)
        np.testing.assert_allclose(sigma0_0, sigma0_1)
        os.unlink(ncfile)

    def incidence_angle(self, rsfile):
        n = Nansat(rsfile)
        inc_min = float(n.get_metadata()['NEAR_RANGE_INCIDENCE_ANGLE'])-0.5
        inc_max = float(n.get_metadata()['FAR_RANGE_INCIDENCE_ANGLE'])+0.5
        inc = n['incidence_angle']
        assert np.all(np.greater_equal(inc[np.isnan(inc)==False], inc_min))
        assert np.all(np.less_equal(inc[np.isnan(inc)==False], inc_max))

    #def test_export_netcdf(self):


