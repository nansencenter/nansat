#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:  18.06.2014
# Last modified:27.08.2014 10:58
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Nansat, Domain
from nansat.nansat import _import_mappers
from mapper_test_archive import DataForTestingMappers

nansatMappers = _import_mappers()

class TestDataTest(unittest.TestCase):
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
                'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs.t00z.master.grbf00',
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
        

if __name__=='__main__':
    unittest.main()





