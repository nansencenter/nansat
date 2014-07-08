#-------------------------------------------------------------------------------
# Name:         test_nansat_archive.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:08.07.2014 11:03
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import os, warnings, timeit

dirname_test_data = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'test_data')
if not os.path.exists(dirname_test_data):
    os.mkdir(dirname_test_data)


'''
    Test data should be contained in an instance of the TestData class so we're
    able to correctly remove downloaded files after the tests
'''
class TestData(object):
    asar = []
    noData = False

    def __init__(self):
        # OBS: SAR and wind data must be added in pairs for each test
        self.get_asar_agulhas()
        if not self.asar:
            self.noData = True

    def get_asar_agulhas(self):
        '''
            Download and assign online datasets
        '''
        asar_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_WSM_1PNPDE20120327_205532_000002143113_00100_52700_6903.N1'
        fname = os.path.basename(asar_agulhas_url)

        asar_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(asar_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + asar_agulhas + ' ' + asar_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(asar_agulhas):
            asar_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.asar.append(asar_agulhas)
