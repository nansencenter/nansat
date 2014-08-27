#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:27.08.2014 10:58
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Nansat, Domain
from nansat.nansat import import_mappers 
import nansat_test_archive as tna

nansatMappers = import_mappers()

class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noAsarData:
            raise ValueError('No test data available')

    def test_mapper_imports(self):
        for mapper in nansatMappers:
            assert type(__import__('nansat.mappers.'+mapper)) is ModuleType

    def test_pixel_functions(self):
        n = Nansat(self.test_data.asar[0])
        if sys.version_info < (2, 7):
            type(n['sigma0_VV']) == np.ndarray
        else:
            self.assertIsInstance(n['sigma0_VV'], np.ndarray)

    def test_resize_by_pixelsize(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(pixelsize=500, eResampleAlg=1)
        if sys.version_info < (2, 7):
            type(n[1]) == np.ndarray
        else:
            self.assertIsInstance(n[1], np.ndarray)

    def test_resize_by_factor(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(0.5, eResampleAlg=1)
        if sys.version_info < (2, 7):
            type(n[1]) == np.ndarray
        else:
            self.assertIsInstance(n[1], np.ndarray)

    def test_resize_by_width(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(width=100, eResampleAlg=1)
        if sys.version_info < (2, 7):
            type(n[1]) == np.ndarray
        else:
            self.assertIsInstance(n[1], np.ndarray)

    def test_resize_by_height(self):
        n = Nansat(self.test_data.asar[0])
        n.resize(height=500, eResampleAlg=1)
        if sys.version_info < (2, 7):
            type(n[1]) == np.ndarray
        else:
            self.assertIsInstance(n[1], np.ndarray)

    def test_reproject(self):
        n = Nansat(self.test_data.asar[0])
        # resize to avoid memory problems
        n.resize(0.1)
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-lle 15 -35 35 -25 -ts 500 500")
        n.reproject(d)
        if sys.version_info < (2, 7):
            type(n[1]) == np.ndarray
        else:
            self.assertIsInstance(n[1], np.ndarray)

    def tearDown(self):
        # if any test plots are created, they could be deleted here
        pass

class MapperAsterL1aTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noAsterL1aData:
            raise ValueError('No test data available')

    def test_mapper_hirlam(self):
        print self.test_data.asterL1a[0]
        n = Nansat(self.test_data.asterL1a[0])

    def tearDown(self):
        pass

class MapperHirlamTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noHirlamData:
            raise ValueError('No test data available')

    def test_mapper_hirlam(self):
        print self.test_data.hirlam[0]
        n = Nansat(self.test_data.hirlam[0])

    def tearDown(self):
        pass

class MapperCosmoskymedTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noCosmoskymedData:
            raise ValueError('No test data available')

    # Proprietary data cannot be shared and will not be available for all users
    @unittest.skipIf(tna.TestData().noCosmoskymedData, 
            "No Cosmo-Skymed data available (this is proprietary and cannot be shared)")
    def test_mapper_cosmoskymed(self):
        print self.test_data.cosmoskymed[0]
        n = Nansat(self.test_data.cosmoskymed[0])

    def tearDown(self):
        pass

class MapperLandsatTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noLandsatData:
            raise ValueError('No test data available')

    def test_mapper_landsat(self):
        print self.test_data.landsat[0]
        n = Nansat(self.test_data.landsat[0])

    def tearDown(self):
        pass

class MapperMeris2Test(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noMerisData:
            raise ValueError('No test data available')

    def test_mapper_meris(self):
        print self.test_data.meris[0]
        n = Nansat(self.test_data.meris[0])

    def tearDown(self):
        pass

class MapperModisL1Test(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noModisL1Data:
            raise ValueError('No test data available')

    def test_mapper_modisL1(self):
        print self.test_data.modisL1[0]
        n = Nansat(self.test_data.modisL1[0])

    def tearDown(self):
        pass

class MapperNcepTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noNcepData:
            raise ValueError('No test data available')

    def test_mapper_generic(self):
        print self.test_data.ncep[0]
        n = Nansat(self.test_data.ncep[0])

    def tearDown(self):
        pass

class MapperRadarsat2Test(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()

    # Proprietary data cannot be shared and will not be available for all users
    @unittest.skipIf(tna.TestData().noRadarsat2Data, 
            "No Radarsat-2 data available (this is proprietary and cannot be shared)")
    def test_mapper_radarsat2(self):
        print self.test_data.radarsat2[0]
        n = Nansat(self.test_data.radarsat2[0])

    def tearDown(self):
        pass

class MapperGenericTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noGenericData:
            raise ValueError('No test data available')

    def test_mapper_generic(self):
        print self.test_data.generic[0]
        n = Nansat(self.test_data.generic[0])

    def tearDown(self):
        pass

class DomainTest(unittest.TestCase):
    def setUp(self):
        self.test_data = tna.TestData()
        if self.test_data.noAsarData:
            raise ValueError('No test data available')

    def test_get_pixelsize_meters(self):
        n = Nansat(self.test_data.asar[0])
        dX, dY = n.get_pixelsize_meters()
        if sys.version_info < (2, 7):
            type(dX) == FloatType
            type(dY) == FloatType
        else:
            self.assertIsInstance(dX, FloatType)
            self.assertIsInstance(dY, FloatType)

    def tearDown(self):
        pass


if __name__=='__main__':
    unittest.main()





