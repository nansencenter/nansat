#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:26.06.2014 10:36
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import *
import numpy as np

from nansat import *

dirname = os.path.dirname(os.path.abspath(__file__))

'''
    Online datasets
'''
asar_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_WSM_1PNPDE20120327_205532_000002143113_00100_52700_6903.N1'
fname = os.path.basename(asar_agulhas_url)

if not os.path.exists(os.path.join(dirname,fname)):
    os.system('curl -so ' + os.path.join(dirname,fname) + ' ' + asar_agulhas_url )
asar_agulhas = os.path.join(dirname,fname)
if not os.path.isfile(asar_agulhas):
    asar_agulhas = None
    warnings.warn( "Could not access ftp-site with test data - contact " \
            "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )

'''
    Test data should be contained in an instance of the TestData class so we're
    able to correctly remove downloaded files after the tests
'''
class TestData():
    asar = []
    radarsat2 = []
    noData = True

    def __init__(self):
        # OBS: SAR and wind data must be added in pairs for each test
        if asar_agulhas:
            self.asar.append(asar_agulhas)
        # check for data from local archives
        if 'rs2' in globals() and rs2:
            self.radarsat2.append(rs2)
        if 'rs2_quad' in globals() and rs2_quad:
            self.radarsat2.append(rs2_quad)
        if self.asar:
            self.noData = False
        if self.radarsat2:
            self.noData = False

    def __del__(self):
        '''
            Delete any downloaded files
        '''
        if os.path.isfile(str(asar_agulhas)):
            os.unlink(str(asar_agulhas))

test_data = TestData()

class NansatTest(unittest.TestCase):
    def setUp(self):
        if test_data.noData:
            raise ValueError('No test data available')

    def test_mapper_imports(self):
        for folder in sys.path:
            for mapper in glob.glob(folder + '/mapper_*.py'):
                mm = os.path.basename(mapper.replace('.py',''))
                assert type(__import__(mm)) is ModuleType

    def test_pixel_functions(self):
        n = Nansat(test_data.asar[0])
        if sys.version_info < (2, 7):
            type(n['sigma0_VV']) == np.ndarray
        else:
            self.assertIsInstance(n['sigma0_VV'], np.ndarray)

    def tearDown(self):
        test_data.__del__()

if __name__=='__main__':
    unittest.main()




