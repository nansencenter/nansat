#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa, Anton Korosov
# Modified:	Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:08.12.2014 14:49
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

# https://nose.readthedocs.org/en/latest/writing_tests.html#test-generators
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

    # Test that methods for returning required metadata are working

    def test_complex_data(self):
        ''' Should open all downloaded files with automatically selected mapper '''
        testData = DataForTestingMappers()
        testData.download_all_test_data()
        for mapperName in testData.mapperData:
            mapperFiles = testData.mapperData[mapperName]
            for mapperFile in mapperFiles:
                if mapperName  == 'radarsat2' or mapperName == 'asar':
                    yield self.exist_intensity_band, mapperFile, mapperName

    def exist_intensity_band(self, mapperFile, mapperName):
        ''' test if intensity bands exist for complex data '''
        n = Nansat(mapperFile, mapperName=mapperName)
        allBandNames = []
        complexBandNames = []
        for iBand in range(n.vrt.dataset.RasterCount):
            iBandName = n.get_metadata(bandID=iBand + 1)['name']
            allBandNames.append(iBandName)
            if '_complex' in iBandName:
                complexBandNames.append(iBandName)

        for iComplexName in complexBandNames:
            assert iComplexName.replace('_complex', '') in allBandNames


if __name__=='__main__':
    unittest.main()





