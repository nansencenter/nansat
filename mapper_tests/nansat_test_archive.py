#-------------------------------------------------------------------------------
# Name:         test_nansat_archive.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:  18.06.2014
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
        self.get_asar()
        if not self.asar:
            self.noAsarData = True

        self.get_asterL1a()
        if not self.asterL1a:
            self.noAsterL1aData = True

        self.get_cosmoskymed()
        if not self.cosmoskymed:
            self.noCosmoskymedData = True

        self.get_hirlam()
        if not self.hirlam:
            self.noHirlamData = True

        self.get_landsat()
        if not self.landsat:
            self.noLandsatData = True

        self.get_meris()
        if not self.meris:
            self.noMerisData = True

        self.get_modisL1()
        if not self.modisL1:
            self.noModisL1Data = True

        self.get_ncep()
        if not self.ncep:
            self.noNcepData = True

        self.get_radarsat2()
        if not self.radarsat2:
            self.noRadarsat2Data = True

        self.get_generic()
        if not self.generic:
            self.noGenericData = True


    def get_asar(self):
        '''
            Download and assign online datasets
        '''
        asar_url = 'ftp://ftp.nersc.no/pub/python_test_data/asar/ASA_WSM_1PNPDE20120327_205532_000002143113_00100_52700_6903.N1'
        fname = os.path.basename(asar_url)

        asar = os.path.join(dirname_test_data, fname)
        if not os.path.exists(asar):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + asar + ' ' + asar_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(asar):
            asar = None
            warnings.warn( "Could not access ASAR on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.asar.append(asar)

    def get_asterL1a(self):
        asterL1a_url = 'ftp://ftp.nersc.no/pub/python_test_data/aster_l1a/AST_L1A_00306192003124632_20120731044546_8073.hdf'
        fname = os.path.basename(asterL1a_url)

        asterL1a = os.path.join(dirname_test_data, fname)
        if not os.path.exists(asterL1a):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + asterL1a + ' ' + asterL1a_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(asterL1a):
            asterL1a = None
            warnings.warn( "Could not access ASTER L1A on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.asterL1a.append(asterL1a)

    def get_cosmoskymed(self):
        cosmoskymed_url = 'ftp://ftp.nersc.no/pub/python_test_data/cosmoskymed/CSKS4_SCS_B_PP_11_CO_LA_FF_20111215040251_20111215040257.h5'
        fname = os.path.basename(cosmoskymed_url)

        cosmoskymed = os.path.join(dirname_test_data, fname)
        if not os.path.exists(cosmoskymed):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + cosmoskymed + ' ' + cosmoskymed_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(cosmoskymed):
            cosmoskymed = None
            warnings.warn( "Could not access CosmoSkyMed on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.cosmoskymed.append(cosmoskymed)

    def get_hirlam(self):
        hirlam_url = 'ftp://ftp.nersc.no/pub/python_test_data/hirlam/DNMI-NEurope.grb'
        fname = os.path.basename(hirlam_url)

        hirlam = os.path.join(dirname_test_data, fname)
        if not os.path.exists(hirlam):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + hirlam + ' ' + hirlam_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(hirlam):
            hirlam = None
            warnings.warn( "Could not access HIRLAM on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.hirlam.append(hirlam)

    def get_landsat(self):
        landsat_url = 'ftp://ftp.nersc.no/pub/python_test_data/landsat/LC81750072013176LGN00.tar.gz'
        fname = os.path.basename(landsat_url)

        landsat = os.path.join(dirname_test_data, fname)
        if not os.path.exists(landsat):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + landsat + ' ' + landsat_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(landsat):
            landsat = None
            warnings.warn( "Could not access Landsat ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.landsat.append(landsat)

    def get_meris(self):
        meris_url = 'ftp://ftp.nersc.no/pub/python_test_data/meris/MER_FRS_1PNUPA20100916_105248_000001012093_00037_44680_8756.N1'
        fname = os.path.basename(meris_url)

        meris = os.path.join(dirname_test_data, fname)
        if not os.path.exists(meris):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + meris + ' ' + meris_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(meris):
            meris = None
            warnings.warn( "Could not access MERIS ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.meris.append(meris)

    def get_modisL1(self):
        modisL1_url = 'ftp://ftp.nersc.no/pub/python_test_data/modis_l1/MOD021KM.A2010105.2120.005.2010106075131.hdf'
        fname = os.path.basename(modisL1_url)

        modisL1 = os.path.join(dirname_test_data, fname)
        if not os.path.exists(modisL1):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + modisL1 + ' ' + modisL1_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(modisL1):
            modisL1 = None
            warnings.warn( "Could not access MODIS L1B ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.modisL1.append(modisL1)

    def get_ncep(self):
        ncep_url = 'ftp://ftp.nersc.no/pub/python_test_data/ncep/gfs/gfs20120328/gfs.t00z.master.grbf00'
        fname = os.path.basename(ncep_url)

        ncep = os.path.join(dirname_test_data, fname)
        if not os.path.exists(ncep):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + ncep + ' ' + ncep_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(ncep):
            ncep = None
            warnings.warn( "Could not access NCEP on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.ncep.append(ncep)

    def get_radarsat2(self):
        radarsat2_url = 'ftp://ftp.nersc.no/pub/python_test_data/radarsat2/RS2_20140716_061819_0076_SCWA_HHHV_SGF_336560_2501_9900957.zip'
        fname = os.path.basename(radarsat2_url)

        radarsat2 = os.path.join(dirname_test_data, fname)
        if not os.path.exists(radarsat2):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + radarsat2 + ' ' + radarsat2_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(radarsat2):
            radarsat2 = None
            warnings.warn( "Could not access Radarsat2 on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.radarsat2.append(radarsat2)

    def get_generic(self):
        generic_url = 'ftp://ftp.nersc.no/pub/python_test_data/generic/mapperTest_generic.tif'
        fname = os.path.basename(generic_url)

        generic = os.path.join(dirname_test_data, fname)
        if not os.path.exists(generic):
            print "Downloading test data"
            start = timeit.timeit()
            os.system('curl -so ' + generic + ' ' + generic_url )
            end = timeit.timeit()
            print end-start

        if not os.path.isfile(generic):
            generic = None
            warnings.warn( "Could not access Generic data on ftp-site with test data - contact " \
                    "morten.stette@nersc.no to get the ftp-server at NERSC restarted" )
        else:
            self.generic.append(generic)


