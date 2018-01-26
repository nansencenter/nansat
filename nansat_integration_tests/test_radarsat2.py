# ------------------------------------------------------------------------------
# Name:     test_radarsat2.py
# Purpose:      End to end testing of Nansat for Radarsat-2
#
# Author:   Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen, Aleksander Vines
#
# Created:  17.04.2015
# Last modified:23.12.2015 13:21
# Copyright:    (c) NERSC
# License:
# ------------------------------------------------------------------------------
import os
import sys
import numpy as np

from nansat.nansat import Nansat
from nansat_integration_tests.mapper_test_archive import DataForTestingMappers


class TestRadarsat(object):

    def test_all_rs2_files(self):
        sys.stderr.write('\ntest_all_rs2_files\n')
        testData = DataForTestingMappers()
        rs2Index = [i for i,
                    v in enumerate(testData.mapperData) if v['mapperName'] == 'radarsat2']
        for index in rs2Index:
            rsfile = testData.mapperData[index][0]
            # yield self.export2thredds, rsfile
            yield self.export, rsfile
            yield self.incidence_angle, rsfile
            yield self.export_band, rsfile
            yield self.resize, rsfile

    def export2thredds(self, rsfile):
        sys.stderr.write('\nexport2thredds:'+rsfile)
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        orig.export2thredds(ncfile, bands={'incidence_angle': {}})
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        np.testing.assert_allclose(inc0, inc1)
        os.unlink(ncfile)

    def export(self, rsfile):
        sys.stderr.write('\nexport:'+rsfile+'\n')
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        sys.stderr.write('\nExporting\n')
        orig.export(ncfile)
        sys.stderr.write('\nOpening Copy\n')
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        sys.stderr.write('\nGet orig grids\n')
        lon0, lat0 = orig.get_geolocation_grids()
        sys.stderr.write('\nGet copy grids\n')
        lon1, lat1 = copy.get_geolocation_grids()
        sys.stderr.write('\nGet orig sigma0_HH\n')
        sigma0_0 = orig['sigma0_HH']
        sys.stderr.write('\nGet copy sigma0_HH\n')
        sigma0_1 = copy['sigma0_HH']
        sys.stderr.write('\nAsserting\n')
        np.testing.assert_allclose(lon0, lon1)
        np.testing.assert_allclose(lat0, lat1)
        # If the next tests fail, it could indicate that the data is flipped
        # check by pyplot.imshow orig vs copy...
        np.testing.assert_allclose(inc0, inc1)
        np.testing.assert_allclose(sigma0_0, sigma0_1)
        os.unlink(ncfile)

    def incidence_angle(self, rsfile):
        sys.stderr.write('\nincidence_angle:'+rsfile+'\n')
        n = Nansat(rsfile)
        inc_min = float(n.get_metadata()['NEAR_RANGE_INCIDENCE_ANGLE'])-0.5
        inc_max = float(n.get_metadata()['FAR_RANGE_INCIDENCE_ANGLE'])+0.5
        inc = n['incidence_angle']
        assert np.all(np.greater_equal(inc[np.isnan(inc) == False], inc_min))
        assert np.all(np.less_equal(inc[np.isnan(inc) == False], inc_max))

    def export_band(self, rsfile):
        sys.stderr.write('\nexport_band:'+rsfile+'\n')
        orig = Nansat(rsfile)
        ncfile = 'test.nc'
        orig.export(ncfile, bands=[orig.get_band_number('incidence_angle')])
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        np.testing.assert_allclose(inc0, inc1)
        os.unlink(ncfile)

    def resize(self, rsfile):
        sys.stderr.write('\nresize:'+rsfile+'\n')
        n = Nansat(rsfile)
        inc_max = float(n.get_metadata()['FAR_RANGE_INCIDENCE_ANGLE'])+0.5
        n.resize(0.5, eResampleAlg=0)
        assert (np.nanmax(n['incidence_angle']) <= inc_max)
        n.undo()
        n.resize(0.5, eResampleAlg=1)
        assert (np.nanmax(n['incidence_angle']) <= inc_max)
        n.undo()
        n.resize(0.5, eResampleAlg=2)
        assert (np.nanmax(n['incidence_angle']) <= inc_max)
        n.undo()
        n.resize(0.5, eResampleAlg=3)
        assert (np.nanmax(n['incidence_angle']) <= inc_max)
        n.undo()
        n.resize(0.5, eResampleAlg=4)
        assert (np.nanmax(n['incidence_angle']) <= inc_max)
        n.undo()
