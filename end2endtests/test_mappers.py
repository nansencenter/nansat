#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa, Anton Korosov
# Modified:	Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:02.07.2015 16:05
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob, datetime
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
# The x-flag results in the test stopping at first failure or error - use it
# for easier debugging:
# nosetests -v -x end2endtests.test_mappers:TestAllMappers.test_mappers
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
                sys.stderr.write('\nMapper '+mapperName+' -> '+fileName+'\n')
                # Test call to Nansat, mapper not specified
                #yield self.open_with_automatic_mapper, fileName, kwargs
                yield self.open_with_nansat, fileName, kwargs
                # Test call to Nansat, mapper specified
                #yield self.open_with_specific_mapper, fileName, mapperName, kwargs
                yield self.open_with_nansat, fileName, mapperName, kwargs
                n = Nansat(fileName, mapperName=mapperName)
                # Test nansat.mapper
                yield self.is_correct_mapper, n, mapperName
                # Test nansat.start_time() and nansat.end_time()
                yield self.has_time, n
                # Test nansat.source() (returns, e.g., Envisat/ASAR)
                yield self.has_source, n
                # Test that SAR objects have sigma0 intensity bands in addition
                # to complex bands
                if n.has_band(
                    'surface_backwards_scattering_coefficient_of_radar_wave'
                        ):
                    yield self.exist_intensity_band, n

    def has_time(self, n):
        assert type(n.start_time())==datetime.datetime
        assert type(n.stop_time())==datetime.datetime

    def has_source(self, n):
        assert type(n.source())==str

    def is_correct_mapper(self, n, mapper):
        assert n.mapper==mapper

    #def open_with_automatic_mapper(self, mapperFile, kwargs):
    #    ''' Perform call to Nansat with each file as a separate test '''
    #    n = Nansat(mapperFile, **kwargs)
    #    n.logger.error('Generic mapper for %s ' % mapperFile)
    #    assert type(n) == Nansat

    #def open_with_specific_mapper(self, mapperFile, mapperName, kwargs):
    #    ''' Perform call to Nansat with each file as a separate test '''
    #    n = Nansat(mapperFile, mapperName=mapperName, **kwargs)
    #    n.logger.error('Mapper %s for %s ' % (mapperName, mapperFile))
    #    assert type(n) == Nansat

    def open_with_nansat(self, file, mapper=None, kwargs={}):
        ''' Perform call to Nansat and check that it returns a Nansat object '''
        if mapper:
            n = Nansat(file, mapperName=mapper, **kwargs)
        else:
            n = Nansat(file, **kwargs)
        # Why this?:
        n.logger.error('Mapper %s for %s ' % (n.mapper, file))
        assert type(n) == Nansat

    def exist_intensity_band(self, n):
        ''' test if intensity bands exist for complex data '''
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
    #for mapper in nansatMappers:
    #    test_name = 'test_%s'%mapper
    unittest.main()




