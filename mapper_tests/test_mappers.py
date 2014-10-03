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
from mapper_test_archive import TestData

nansatMappers = _import_mappers()

class TestDataTest(unittest.TestCase):
    def test_create_test_data(self):
        ''' should create TestData instance '''
        t = TestData()
        self.assertTrue(hasattr(t, 'mapperData'))
        self.assertTrue(hasattr(t, 'testDataDir'))

    def test_testDataDir_from_env(self):
        ''' should create TestData instance '''
        fakeDir = '/fake/dir/to/test/data'
        os.environ['MAPPER_TEST_DATA_DIR'] = fakeDir
        t = TestData()
        self.assertEqual(t.testDataDir, fakeDir)

    def test_testDataDir_exists(self):
        ''' should create TestData instance '''
        t = TestData()
        self.assertTrue(os.path.exists(t.testDataDir))

    def test_download_file(self):
        ''' Should download the selected file and put into mapperData'''
        t = TestData()
        t.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs/gfs20120328/gfs.t00z.master.grbf00',
                'ncep')
        self.assertTrue('ncep' in t.mapperData)
        self.assertEqual(type(t.mapperData['ncep']), list)
        for ifile in t.mapperData['ncep']:
            self.assertTrue(os.path.exists(ifile))


class AllMappersTest(unittest.TestCase):
    def setUp(self):
        # download all data
        self.testData = TestData()
        self.testData.download_all_test_data()

    def test_automatic_mapper(self):
        ''' Should open all downloaded files with automatically selected mapper '''
        for mapper in self.testData.mapperData:
            mapperFiles = self.testData.mapperData[mapper]
            for mapperFile in mapperFiles:
                n = Nansat(mapperFile)

    def test_specific_mapper(self):
        ''' Should open all downloaded files with automatically selected mapper '''
        for mapper in self.testData.mapperData:
            mapperFiles = self.testData.mapperData[mapper]
            for mapperFile in mapperFiles:
                n = Nansat(mapperFile, mapperName=mapper)

if __name__=='__main__':
    unittest.main()





