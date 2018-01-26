from __future__ import absolute_import, unicode_literals, division
import unittest
import sys
import datetime
import json

import pythesint as pti

from nansat import Nansat
from nansat.nansat import _import_mappers
from nansat_integration_tests.mapper_test_archive import DataForTestingMappers, DataForTestingOnlineMappers


class TestDataForTestingMappers(unittest.TestCase):
    def test_create_test_data(self):
        ''' should create TestData instance '''
        t = DataForTestingMappers()
        self.assertTrue(hasattr(t, 'mapperData'))


# https://nose.readthedocs.org/en/latest/writing_tests.html#test-generators
# The x-flag results in the test stopping at first failure or error - use it
# for easier debugging:
# nosetests -v -x integration_tests.test_mappers:TestAllMappers.test_mappers_basic
class TestAllMappers(object):

    @classmethod
    def setup_class(cls):
        ''' Download testing data '''
        cls.testData = DataForTestingMappers()

    def test_mappers_basic(self):
        ''' Basic mapper test '''
        for kwargs in self.testData.mapperData:
            fileName = kwargs.pop('fileName')
            mapperName = kwargs.pop('mapperName')
            # Test call to Nansat, mapper not specified
            sys.stderr.write('\n -> '+fileName+'\n')
            yield self.open_with_nansat, fileName, None, kwargs
            # Test call to Nansat, mapper specified
            sys.stderr.write('\n'+mapperName+' -> '+fileName+'\n')
            yield self.open_with_nansat, fileName, mapperName, kwargs

    def open_with_nansat(self, filePath, mapper=None, kwargs=None):
        ''' Ensures that you can open the filePath as a Nansat object '''
        if kwargs is None:
            kwargs = {}

        try:
            if mapper:
                n = Nansat(filePath, mapperName=mapper, **kwargs)
            else:
                n = Nansat(filePath, **kwargs)
        except Exception as e:
            raise Exception('%s: %s'%(filePath, e.message))
        assert type(n) == Nansat

    def test_mappers_advanced(self):
        ''' Run tests to check DIF/GCMD metadata content in all mappers'''
        for dd in self.testData.mapperData:
            sys.stderr.write('\nMapper '+dd['mapperName']+' -> '+dd['fileName']+'\n')
            n = Nansat(dd['fileName'], mapperName=dd['mapperName'])
            yield self.is_correct_mapper, n, dd['mapperName']
            yield self.has_metadata_time_coverage_start, n
            yield self.has_metadata_time_coverage_end, n
            yield self.has_metadata_platform, n
            yield self.has_metadata_instrument, n

            # Test that SAR objects have sigma0 intensity bands in addition
            # to complex bands
            if n.has_band(
                'surface_backwards_scattering_coefficient_of_radar_wave'
                    ):
                yield self.exist_intensity_band, n

    def has_metadata_time_coverage_start(self, n):
        ''' Has start time '''
        assert type(n.time_coverage_start) == datetime.datetime

    def has_metadata_time_coverage_end(self, n):
        assert type(n.time_coverage_end) == datetime.datetime

    def has_metadata_platform(self, n):
        meta1 = json.loads(n.get_metadata('platform'))
        meta1ShortName = meta1['Short_Name']
        meta2 = pti.get_gcmd_platform(meta1ShortName)

        assert type(meta1) == dict
        assert meta1 == meta2

    def has_metadata_instrument(self, n):
        meta1 = json.loads(n.get_metadata('instrument'))
        meta1ShortName = meta1['Short_Name']
        meta2 = pti.get_gcmd_instrument(meta1ShortName)

        assert type(meta1) == dict
        assert meta1 == meta2

    def is_correct_mapper(self, n, mapper):
        assert n.mapper == mapper

    def exist_intensity_band(self, n):
        ''' Test if intensity bands exist for complex data '''
        allBandNames = []
        complexBandNames = []
        for iBand in range(n.vrt.dataset.RasterCount):
            iBandName = n.get_metadata(band=iBand + 1)['name']
            allBandNames.append(iBandName)
            if '_complex' in iBandName:
                complexBandNames.append(iBandName)

        for iComplexName in complexBandNames:
            assert iComplexName.replace('_complex', '') in allBandNames

class TestOnlineMappers(TestAllMappers):
    @classmethod
    def setup_class(cls):
        ''' Download testing data '''
        cls.testData = DataForTestingOnlineMappers()


if __name__ == '__main__':
    # nansatMappers = _import_mappers()
    # for mapper in nansatMappers:
    #     test_name = 'test_%s'%mapper
    unittest.main()

