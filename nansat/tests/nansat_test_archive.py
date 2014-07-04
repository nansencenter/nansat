#-------------------------------------------------------------------------------
# Name:         test_nansat_archive.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:04.07.2014 14:57
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import os, warnings, timeit

'''
    Test data should be contained in an instance of the TestData class so we're
    able to correctly remove downloaded files after the tests
'''
class TestData():
    asar = []
    noData = False
    dirname = os.path.dirname(os.path.abspath(__file__))

    def __init__(self):
        # OBS: SAR and wind data must be added in pairs for each test
        self._get_online_datasets()
        if not self.asar:
            self.noData = True

    def delete_downloaded(self):
        '''
            Delete any downloaded files
        '''
        for a in self.asar:
            if os.path.isfile(a):
                os.unlink(a)
        self.asar = []
        self.noData = True

    def _get_online_datasets(self):
        '''
            Download and assign online datasets
        '''
        asar_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_WSM_1PNPDE20120327_205532_000002143113_00100_52700_6903.N1'
        fname = os.path.basename(asar_agulhas_url)
        
        if not os.path.exists(os.path.join(self.dirname,fname)):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + os.path.join(self.dirname,fname) + ' ' + asar_agulhas_url )
            end = timeit.timeit()
            print end-start

        asar_agulhas = os.path.join(self.dirname,fname)
        if not os.path.isfile(asar_agulhas):
            asar_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.asar.append(asar_agulhas)
        




