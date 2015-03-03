#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa, Anton Korosov
# Modified:	Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:03.03.2015 15:10
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
        for ifile in t.mapperData['ncep']:
            self.assertTrue(os.path.exists(ifile))


class TestAllMappers(object):
    def test_automatic_mapper(self):
        ''' Should open all downloaded files with automatically selected mapper '''
        testData = DataForTestingMappers()
        testData.download_all_test_data()
        for mapper in testData.mapperData:
            mapperFiles = testData.mapperData[mapper]
            for mapperFile in mapperFiles:
                print mapperFile
                yield self.open_with_automatic_mapper, mapperFile

    def open_with_automatic_mapper(self, mapperFile):
        ''' Perform call to Nansat with each file as a separate test '''
        n = Nansat(mapperFile)
        assert type(n) == Nansat

    def test_specific_mapper(self):
        ''' Should open all downloaded files with automatically selected mapper '''
        testData = DataForTestingMappers()
        testData.download_all_test_data()
        for mapperName in testData.mapperData:
            mapperFiles = testData.mapperData[mapperName]
            for mapperFile in mapperFiles:
                print mapperName, '->', mapperFile
                yield self.open_with_specific_mapper, mapperFile, mapperName

    def open_with_specific_mapper(self, mapperFile, mapperName):
        ''' Perform call to Nansat with each file as a separate test '''
        n = Nansat(mapperFile, mapperName=mapperName)
        assert type(n) == Nansat

class TestRadarsat(object):

    def test_all_rs2_files(self):
        testData = DataForTestingMappers()
        testData.download_all_test_data()
        for rsfile in testData.mapperData['radarsat2']:
            yield self.test_incidence_angle, rsfile
            #yield self.test_export2thredds, rsfile
            yield self.test_export, rsfile

    def test_export2thredds(self, rsfile):
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        orig.export2thredds(ncfile, bands = {'incidence_angle': {}})
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        np.testing.assert_allclose(inc0, inc1)
        os.unlink(ncfile)

    def test_export(self, rsfile):
        ncfile = 'test.nc'
        orig = Nansat(rsfile)
        orig.export(ncfile)
        copy = Nansat(ncfile)
        inc0 = orig['incidence_angle']
        inc1 = copy['incidence_angle']
        np.testing.assert_allclose(inc0, inc1)
        os.unlink(ncfile)
        
    def test_incidence_angle(self, rsfile):
        n = Nansat(rsfile)
        inc_min = float(n.get_metadata()['NEAR_RANGE_INCIDENCE_ANGLE'])
        inc_max = float(n.get_metadata()['FAR_RANGE_INCIDENCE_ANGLE'])
        inc = n['incidence_angle']
        assert np.all(np.greater_equal(inc[np.isnan(inc)==False], inc_min))
        assert np.all(np.less_equal(inc[np.isnan(inc)==False], inc_max))

    #def test_export_netcdf(self):

if __name__=='__main__':
    unittest.main()





