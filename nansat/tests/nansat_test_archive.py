#-------------------------------------------------------------------------------
# Name:         test_nansat_archive.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:27.06.2014 11:04
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import os

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
    noData = False

    def __init__(self):
        # OBS: SAR and wind data must be added in pairs for each test
        if asar_agulhas:
            self.asar.append(asar_agulhas)
        if not self.asar:
            self.noData = True

    def __del__(self):
        '''
            Delete any downloaded files
        '''
        for a in self.asar:
            if os.path.isfile(a):
                os.unlink(a)
        self.asar = []
        self.noData = True




