#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa, Anton Korosov
# Modified:     Anton Korosov
#
# Created:      18.06.2014
# Last modified:17.04.2015 15:23
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Nansat, Domain
from nansat.nansat import _import_mappers
from mapper_test_archive import DataForTestingMappers

nansatMappers = _import_mappers()

class TestDataForTestingMappers(unittest.TestCase):
    def test_create_test_data(self):
        ''' should create TestData instance '''
        t = DataForTestingMappers()
        self.assertTrue(hasattr(t, 'mapperData'))
        self.assertTrue(hasattr(t, 'testDataDir'))

    def test_testDataDir_from_env(self):
        ''' should create TestData instance '''
        fakeDir = '/fake/dir/to/test/data'
        os.environ['MAPPER_TEST_DATA_DIR'] = fakeDir
        t = DataForTestingMappers()
        self.assertEqual(t.testDataDir, fakeDir)

    def test_testDataDir_exists(self):
        ''' should create TestData instance '''
        t = DataForTestingMappers()
        self.assertTrue(os.path.exists(t.testDataDir))

    def test_download_file(self):
        ''' Should download the selected file and put into mapperData'''
        t = DataForTestingMappers()
        t.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs20120328.t00z.master.grbf00',
                'ncep')
        self.assertTrue('ncep' in t.mapperData)
        self.assertEqual(type(t.mapperData['ncep']), list)
        for ifile, kwa in t.mapperData['ncep']:
            self.assertTrue(os.path.exists(ifile))

# https://nose.readthedocs.org/en/latest/writing_tests.html#test-generators
class TestAllMappers(object):

    @classmethod
    def setup_class(cls):
        ''' Download testing data '''
        cls.testData = DataForTestingMappers()
        cls.testData.download_all_test_data()

    def test_mappers(self):
        ''' Run similar tests for all mappers '''
        for mapperName in self.testData.mapperData:
            mapperParams = self.testData.mapperData[mapperName]
            for fileName, kwargs in mapperParams:
                print mapperName, '->', fileName
                # Test call to Nansat, mapper not specified
                yield self.open_with_automatic_mapper, fileName, kwargs
                # Test call to Nansat, mapper specified
                yield self.open_with_specific_mapper, fileName, mapperName, kwargs

    def open_with_automatic_mapper(self, mapperFile, kwargs):
        ''' Perform call to Nansat with each file as a separate test '''
        n = Nansat(mapperFile, **kwargs)
        n.logger.error('Generic mapper for %s ' % mapperFile)
        assert type(n) == Nansat

    def open_with_specific_mapper(self, mapperFile, mapperName, kwargs):
        ''' Perform call to Nansat with each file as a separate test '''
        n = Nansat(mapperFile, mapperName=mapperName, **kwargs)
        n.logger.error('Mapper %s for %s ' % (mapperName, mapperFile))
        assert type(n) == Nansat


class TestRadarsat(object):

    def test_all_rs2_files(self):
        ''' Run tests for all Radarsat2 data '''
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
        orig.resize(0.05)
        orig.export2thredds(ncfile, bands = {'incidence_angle': {}})
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        np.testing.assert_allclose(inc0, inc1)
        os.unlink(ncfile)

    def export(self, rsfile):
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        orig.resize(0.05)
        orig.export(ncfile)
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        lon0, lat0 = orig.get_geolocation_grids()
        lon1, lat1 = copy.get_geolocation_grids()
        sigma0_0 = orig['sigma0_HH']
        sigma0_1 = copy['sigma0_HH']
        np.testing.assert_allclose(inc0, inc1, rtol=1e-3)
        # Make sure data is not flipped
        # OBS: tolerance is reduced to test flipping - see other tests
        # regarding geolocation of exported vs original
        np.testing.assert_allclose(lon0, lon1, rtol=1e-2)
        np.testing.assert_allclose(lat0, lat1, rtol=1e-2)
        np.testing.assert_allclose(sigma0_0, sigma0_1)
        os.unlink(ncfile)

    def incidence_angle(self, rsfile):
        n = Nansat(rsfile)
        inc_min = float(n.get_metadata()['NEAR_RANGE_INCIDENCE_ANGLE'])-0.5
        inc_max = float(n.get_metadata()['FAR_RANGE_INCIDENCE_ANGLE'])+0.5
        inc = n['incidence_angle']
        assert np.all(np.greater_equal(inc[np.isnan(inc)==False], inc_min))
        assert np.all(np.less_equal(inc[np.isnan(inc)==False], inc_max))

## Test Generator with unittests:
## http://stackoverflow.com/questions/32899/how-to-generate-dynamic-parametrized-unit-tests-in-python
#class TestMapper(unittest.TestCase):
#    def setUp(self):
#        self.testData = DataForTestingMappers()
#        self.testData.download_all_test_data()
#
#def test_generator():
#    def test_import(self):
#

if __name__=='__main__':
    #for mapper in nansatMappers:
    #    test_name = 'test_%s'%mapper
    unittest.main()




