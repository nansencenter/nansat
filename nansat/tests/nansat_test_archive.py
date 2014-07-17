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
    asterL1a = []
    cosmoskymed = []
    hirlam = []
    landsat = []
    meris = []
    modisL1 = []
    ncep = []
    generic = []
    radarsat2 = []
    noAsarData = False
    noAsterL1aData = False
    noCosmoskymedData = False
    noLandsatData = False
    noHirlamData = False
    noMerisData = False
    noModisL1Data = False
    noNcepData = False
    noGenericData = False
    noRadarsat2Data =False

    def __init__(self):
        # OBS: SAR and wind data must be added in pairs for each test
        self.get_asar_agulhas()
        if not self.asar:
            self.noAsarData = True

        self.get_asterL1a_agulhas()
        if not self.asterL1a:
            self.noAsterL1aData = True

        self.get_cosmoskymed_agulhas()
        if not self.cosmoskymed:
            self.noCosmoskymedData = True

        self.get_hirlam_agulhas()
        if not self.hirlam:
            self.noHirlamData = True

        self.get_landsat_agulhas()
        if not self.landsat:
            self.noLandsatData = True

        self.get_meris_agulhas()
        if not self.meris:
            self.noMerisData = True

        self.get_modisL1_agulhas()
        if not self.modisL1:
            self.noModisL1Data = True

        self.get_ncep_agulhas()
        if not self.ncep:
            self.noNcepData = True

        self.get_radarsat2_agulhas()
        if not self.radarsat2:
            self.noRadarsat2Data = True

        self.get_generic_agulhas()
        if not self.generic:
            self.noGenericData = True


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

    def get_asterL1a_agulhas(self):
        asterL1a_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/aster_l1a/AST_L1A_00306192003124632_20120731044546_8073.hdf'
        fname = os.path.basename(asterL1a_agulhas_url)

        asterL1a_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(asterL1a_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + asterL1a_agulhas + ' ' + asterL1a_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(asterL1a_agulhas):
            asterL1a_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.asterL1a.append(asterL1a_agulhas)

    def get_cosmoskymed_agulhas(self):
        cosmoskymed_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/cosmoskymed/CSKS4_SCS_B_PP_11_CO_LA_FF_20111215040251_20111215040257.h5'
        fname = os.path.basename(cosmoskymed_agulhas_url)

        cosmoskymed_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(cosmoskymed_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + cosmoskymed_agulhas + ' ' + cosmoskymed_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(cosmoskymed_agulhas):
            cosmoskymed_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.cosmoskymed.append(cosmoskymed_agulhas)

    def get_hirlam_agulhas(self):
        hirlam_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/hirlam/DNMI-NEurope.grb'
        fname = os.path.basename(hirlam_agulhas_url)

        hirlam_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(hirlam_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + hirlam_agulhas + ' ' + hirlam_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(hirlam_agulhas):
            hirlam_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.hirlam.append(hirlam_agulhas)

    def get_landsat_agulhas(self):
        landsat_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/landsat/LC81750072013176LGN00.tar.gz'
        fname = os.path.basename(landsat_agulhas_url)

        landsat_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(landsat_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + landsat_agulhas + ' ' + landsat_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(landsat_agulhas):
            landsat_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.landsat.append(landsat_agulhas)

    def get_meris_agulhas(self):
        meris_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/meris/MER_FRS_1PNUPA20100916_105248_000001012093_00037_44680_8756.N1'
        fname = os.path.basename(meris_agulhas_url)

        meris_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(meris_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + meris_agulhas + ' ' + meris_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(meris_agulhas):
            meris_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.meris.append(meris_agulhas)

    def get_modisL1_agulhas(self):
        modisL1_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/modis_l1/MOD021KM.A2010105.2120.005.2010106075131.hdf'
        fname = os.path.basename(modisL1_agulhas_url)

        modisL1_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(modisL1_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + modisL1_agulhas + ' ' + modisL1_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(modisL1_agulhas):
            modisL1_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.modisL1.append(modisL1_agulhas)

    def get_ncep_agulhas(self):
        ncep_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs/gfs20120328/gfs.t00z.master.grbf00'
        fname = os.path.basename(ncep_agulhas_url)

        ncep_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(ncep_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + ncep_agulhas + ' ' + ncep_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(ncep_agulhas):
            ncep_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.ncep.append(ncep_agulhas)

    def get_radarsat2_agulhas(self):
        radarsat2_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/radarsat2/RS2_20140716_061819_0076_SCWA_HHHV_SGF_336560_2501_9900957.zip'
        fname = os.path.basename(radarsat2_agulhas_url)

        radarsat2_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(radarsat2_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + radarsat2_agulhas + ' ' + radarsat2_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(radarsat2_agulhas):
            radarsat2_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.radarsat2.append(radarsat2_agulhas)

    def get_generic_agulhas(self):
        generic_agulhas_url = 'ftp://ftp.nersc.no/pub/python_test_data/generic/mapperTest_generic.tif'
        fname = os.path.basename(generic_agulhas_url)

        generic_agulhas = os.path.join(dirname_test_data, fname)
        if not os.path.exists(generic_agulhas):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + generic_agulhas + ' ' + generic_agulhas_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(generic_agulhas):
            generic_agulhas = None
            warnings.warn( "Could not access ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.generic.append(generic_agulhas)


