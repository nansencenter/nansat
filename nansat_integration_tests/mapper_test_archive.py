# ------------------------------------------------------------------------------
# Name:         mapper_test_archive.py
# Purpose:      To discover test data for the integration tests
#
# Author:       Anton Korosov, Morten Wergeland Hansen, Asuka Yamakawa
# Modified:     Morten Wergeland Hansen, Aleksander Vines
#
# Created:      2014-06-18
# Last modified:2015-12-28 13:36
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
import os
import glob
import warnings

class DataForTestingMappers(object):
    def __init__(self):
        ''' Find test files and corresponding mapper names '''
        existingTestFiles = self.find_existing_files()
        self.mapperData = self.identify_mappers(existingTestFiles)

    def find_existing_files(self):
        ''' Find all files for testing inside MAPPER_TEST_DATA_DIR'''
        testFiles = []

        testDataEnv = os.getenv('MAPPER_TEST_DATA_DIR')
        if testDataEnv is None:
            warnings.warn('MAPPER_TEST_DATA_DIR is not defined')
        else:
            testDataDirs = testDataEnv.split(':')
            for testDataDir in testDataDirs:
                if os.path.isdir(testDataDir):
                    testFiles += glob.glob(os.path.join(testDataDir, '*', '*'))

        testFiles = [f for f in testFiles if self.readable(f)]

        return testFiles

    def identify_mappers(self, testFiles):
        ''' Get the name of the mapper from the sub-directory name '''

        return [{'fileName' : testFile,
                 'mapperName' : os.path.split(os.path.split(testFile)[0])[1]}
                for testFile in testFiles]

    def readable(self, testFile):
        ''' Test if file is readable at OS level '''
        if not os.path.exists(testFile):
            return False
        if not os.access(testFile, os.R_OK):
            return False
        if os.stat(testFile).st_size == 0:
            return False
        if os.path.isdir(testFile):
            return False

        return True

class DataForTestingOnlineMappers(object):
    mapperData = [
        {
        'fileName' : 'http://dap.ceda.ac.uk/data/neodc/esacci/sst/data/lt/Analysis/L4/v01.1/2010/05/01/20100501120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.1.nc',
        'mapperName': 'opendap_sstcci',
        'bands': ['analysed_sst'],

        },{
        'fileName' : 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-8DAY',
        'mapperName' : 'opendap_occci',
        'date' : '2010-01-01',
        'bands': ['chlor_a'],
        },{
        'fileName' : 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-MONTHLY',
        'mapperName' : 'opendap_occci',
        'date' : '2010-01-01',
        'bands': ['chlor_a'],


        },{
        'fileName' : 'http://www.ifremer.fr/opendap/cerdap1/globcurrent/v2.0/global_025_deg/total_hs/2010/002/20100102000000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc',
        'mapperName' : 'opendap_globcurrent',
        'bands': ['eastward_eulerian_current_velocity'],

        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh/osisaf-nh_aggregated_ice_concentration_nh_polstere-100_197810010000.nc',
        'mapperName' : 'opendap_osisaf',
        'bands': ['ice_conc_avg'],
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh-agg',
        'mapperName' : 'opendap_osisaf',
        'bands': ['ice_conc_avg'],
        'date': '1980-07-22',
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/conc/2016/04/ice_conc_sh_polstere-100_multi_201604261200.nc',
        'mapperName' : 'opendap_osisaf',
        'bands': ['ice_conc'],
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr/merged/2016/04/ice_drift_nh_polstere-625_multi-oi_201604151200-201604171200.nc',
        'mapperName' : 'opendap_osisaf',
        'bands': ['dX', 'dY'],
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/2016/04/ice_type_nh_polstere-100_multi_201604151200.nc',
        'mapperName' : 'opendap_osisaf',
        'bands': ['ice_type'],
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/edge/2016/04/ice_edge_nh_polstere-100_multi_201604241200.nc',
        'mapperName' : 'opendap_osisaf',
        'bands': ['ice_edge'],

        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03/20121001000000-METNO-L4_GHRSST-SSTfnd-METNO_OI-ARC-v02.0-fv01.0.nc',
        'mapperName' : 'opendap_siwtacsst',
        'bands': ['analysed_sst'],
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03_V1/20120808-METNO-L4UHfnd-ARC-v01-fv01-METNO_OI.nc',
        'mapperName' : 'opendap_siwtacsst',
        'bands': ['analysed_sst'],
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/sea_ice/SST-METNO-ARC-SST_L4-OBS-V2-V1/sst_arctic_aggregated',
        'mapperName' : 'opendap_siwtacsst',
        'bands': ['analysed_sst'],
        'date': '2012-08-08'
        },{
        'fileName' : 'https://thredds.met.no/thredds/dodsC/aromearcticarchive/2017/10/30/arome_arctic_full_2_5km_20171030T21Z.nc',
        'mapperName' : 'opendap_arome',
        'bands': ['x_wind_10m', 'y_wind_10m'],
        'date': '2017-10-30'
        },{'fileName' : 'http://thredds.met.no/thredds/dodsC/meps25epsarchive/2017/10/26/meps_mbr1_pp_2_5km_20171026T06Z.nc',
        'mapperName' : 'opendap_arome',
        'bands': ['x_wind_10m', 'y_wind_10m'],
        'date': '2017-10-26T07:00'
        },{
        'fileName' : 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/L4/GLOB/UKMO/OSTIA/2016/006/20160106-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2',
        'mapperName' : 'opendap_ostia',
        'bands': ['analysed_sst', 'mask'],
        'date': '2016-01-06'
        },{
        'fileName' : 'http://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/2017/10/29/MyWave_wam4_WAVE_20171029T18Z.nc',
        'mapperName' : 'opendap_mywave4km',
        'bands': ['hs'],
        'date': '2017-10-29T18:00Z'
        }
        ]
