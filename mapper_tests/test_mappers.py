#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa, Anton Korosov
# Modified:	Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:17.12.2014 15:17
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

    @classmethod
    def setup_class(cls):
        cls.testData = DataForTestingMappers()
        cls.testData.download_all_test_data()

    def test_mappers(self):
        for mapper in self.testData.mapperData:
            mfiles = self.testData.mapperData[mapper]
            for f in mfiles:
                print mapper, '->', f
                # Test call to Nansat, mapper not specified
                yield self.open_with_nansat, f
                # Test call to Nansat, mapper specified
                yield self.open_with_nansat, f, mapper
                # Test nansat.mapper()
                # Test nansat.start_time()
                # Test nansat.end_time()
                # Test nansat.source() (returns, e.g., Envisat/ASAR)
                # Test that SAR objects have sigma0 intensity bands in addition
                # to complex bands
                if mapper == 'radarsat2' or mapper == 'asar':
                    yield self.exist_intensity_band, f, mapper


    def open_with_nansat(self, file, mapper=None):
        ''' Perform call to Nansat and check that it returns a Nansat object '''
        if mapper:
            n = Nansat(file, mapperName=mapper)
        else:
            n = Nansat(file)
        assert type(n) == Nansat

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




