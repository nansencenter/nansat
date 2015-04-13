#-------------------------------------------------------------------------------
# Name:         test_nansat_archive.py
# Purpose:      To test nansat
#
# Author:       Anton Korosov, Morten Wergeland Hansen, Asuka Yamakawa
# Modified:     Anton Korosov
#
# Created:      18.06.2014
# Last modified:13.04.2015 15:17
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#-------------------------------------------------------------------------------
import os, warnings, time


class DataForTestingMappers(object):
    ''' Download test data and keep info about each file '''
    mapperData = None
    mustHaveFiles = [
        'ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1',
        ]

    def __init__(self):
        ''' Set directory to store test data

        If MAPPER_TEST_DATA_DIR is in the environment its value will be used
        This is convenient for testing localy and sharing downloaded
        data among several users on the server
        '''
        self.testDataDir = os.getenv('MAPPER_TEST_DATA_DIR')
        if self.testDataDir is None:
            self.testDataDir = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'test_data')

    def download_all_test_data(self):
        ''' Download test data for each mapper '''

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_IMS_1PNIPA20100411_101715_000000162088_00280_42418_0338.N1',
                'asar')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1',
                'asar')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/aster_l1a/AST_L1A_00306192003124632_20120731044546_8073.hdf',
                'aster_l1a')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/csks/CSKS4_SCS_B_PP_11_CO_LA_FF_20111215040251_20111215040257.h5',
                'csks')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/hirlam/DNMI-NEurope.grb',
                'hirlam')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/landsat/LT52280012006208KIS00.tar.gz',
                'landsat')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/landsat/LM52200021984220AAA04.tar.gz',
                'landsat')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/landsat/LE72210032011229EDC00.tar.gz',
                'landsat')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/landsat/LE72210032011229EDC00.tar.gz',
                'landsat',
                resolution='high')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/landsat/LC81750072013176LGN00.tar.gz',
                'landsat')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/landsat/LC81750072013176LGN00.tar.gz',
                'landsat',
                resolution='high')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/meris_l1/MER_FRS_1PNPDK20110503_105638_000001833102_00109_47968_7898.N1',
                'meris_l1')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/modis_l1/MOD021KM.A2010105.2120.005.2010106075131.hdf',
                'modis_l1')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs20120328.t00z.master.grbf00',
                'ncep')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/radarsat2/RS2_20090227_063055_0080_SCWA_HHHV_SCW_30853_0000_1897838',
                'radarsat2')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/radarsat2/RS2_20111109_060616_0045_SCNA_HHHV_SGF_164373_9871_6913894',
                'radarsat2')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/radarsat2/RS2_20140723_161314_0003_U20_VV_SLC_337855_2455_9614320',
                'radarsat2')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/radarsat2/RS2_OK57403_PK539140_DK477416_SCWA_20141022_152035_HH_SGF.ZIP',
                'radarsat2')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/generic/mapperTest_generic.tif',
                'generic')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/obpg_l2/A2014275111000.L2_LAC.NorthNorwegianSeas.hdf',
                'obpg_l2')

        self.download_test_file(
                'ftp://ftp.nersc.no/pub/python_test_data/amsr2_l1r/GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5',
                'amsr2_l1r')

    def download_test_file(self, inputURL, mapperName, **kwargs):
        ''' Download one file for one mapper

        For the given URL and mapper name
        Create local dir with name ./test_data/mapper_name/mapper_file.ext'
        If the downloaded file does not already exist:
            download the file into the dir
        Keep the filepath in self.mapper_data[mapper_name]

        Parameters:
        -----------
            inputUrl : str
                valid URL with the test file to download
            mapperName : str
                name of the mapper for which the data is downloaded

        ModifIes:
        ---------
            self.mapper_data : dict
                adds new <mapper_name> : [<testFileName>]
                or appends <testFileName> to the existing key

        '''
        fileName = os.path.basename(inputURL)
        mapperDir = os.path.split(os.path.split(inputURL)[0])[1]
        mapperDataDir = os.path.join(self.testDataDir, mapperDir)
        mapperFileName = os.path.join(mapperDataDir, fileName)

        if not os.path.exists(mapperDataDir):
            os.makedirs(mapperDataDir)

        if not os.path.exists(mapperFileName):
            print "Downloading %s " % mapperFileName
            t0 = time.time()
            os.system('curl -so ' + mapperFileName + ' ' + inputURL )
            print time.time() - t0

        if not os.path.exists(mapperFileName):
            if os.path.basename(mapperFileName) in self.mustHaveFiles:
                warnings.warn( """
                    Could not access ftp-site with test data - contact
                    morten.stette@nersc.no to get the ftp-server at NERSC restarted""")
            else:
                warnings.warn('%s: file is not available for download' % mapperFName)
        else:
            # create entry for that mapper
            if self.mapperData is None:
                self.mapperData = {}
            # add file and kwargs for that mapper
            if mapperName in self.mapperData:
                self.mapperData[mapperName].append((mapperFileName, kwargs))
            else:
                self.mapperData[mapperName] = [(mapperFileName, kwargs)]

