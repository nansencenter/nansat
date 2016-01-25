# Name:         mapper_aapp.py
# Purpose:      Nansat mapper for AAPP output
# Author:       Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

# Description of file format:
# http://research.metoffice.gov.uk/research/interproj/nwpsaf/aapp/NWPSAF-MF-UD-003_Formats.pdf (page 8-)
import sys
import struct
import datetime

from nansat.tools import WrongMapperError
from nansat.vrt import VRT, GeolocationArray

satIDs = {4: 'NOAA-15', 2: 'NOAA-16', 6: 'NOAA-17', 7: 'NOAA-18', 8: 'NOAA-19',
          11: 'Metop-B (Metop-1)', 12: 'Metop-A (Metop-2)',
          13: 'Metop-C (Metop-3)'}
dataFormats = {1: 'LAC', 2: 'GAC', 3: 'HRPT'}
dataSetQualityIndicatorOffset = 114
recordLength = 22016
headerLength = recordLength
imageOffset = headerLength + 1264


class Mapper(VRT):
    ''' VRT with mapping of WKV for AVHRR L1B output from AAPP '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):

        ########################################
        # Read metadata from binary file
        ########################################
        try:
            fp = open(fileName, 'rb')
        except IOError:
            raise WrongMapperError
        fp.seek(72)

        try:
            satNum = int(struct.unpack('<H', fp.read(2))[0])
        except:
            raise WrongMapperError

        if satNum >= 11:
            isMetop = True
        else:
            isMetop = False

        if satNum in satIDs.keys():
            satID = satIDs[satNum]
        else:
            raise WrongMapperError

        fp.seek(76)
        dataFormatNum = int(struct.unpack('<H', fp.read(2))[0])
        if dataFormatNum in dataFormats.keys():
            dataFormat = dataFormats[dataFormatNum]
        else:
            raise WrongMapperError

        fp.seek(dataSetQualityIndicatorOffset + 14)
        numScanLines = int(struct.unpack('<H', fp.read(2))[0])
        numCalibratedScanLines = int(struct.unpack('<H', fp.read(2))[0])
        missingScanLines = int(struct.unpack('<H', fp.read(2))[0])
        if missingScanLines != 0:
            print('WARNING: Missing scanlines: ' + str(missingScanLines))

        ##################
        # Read time
        ##################
        fp.seek(84)
        year = int(struct.unpack('<H', fp.read(2))[0])
        dayofyear = int(struct.unpack('<H', fp.read(2))[0])
        millisecondsOfDay = int(struct.unpack('<l', fp.read(4))[0])
        time = (datetime.datetime(year, 1, 1) +
            datetime.timedelta(dayofyear-1, milliseconds=millisecondsOfDay))

        ##################################
        # Read calibration information
        ##################################
        #IRcalibration = {}
        #fp.seek(202)
        #avh_h_irttcoef=[]
        #for i in range(24):
        #    avh_h_irttcoef.append(int(struct.unpack('<H', fp.read(2))[0]))
        ##print avh_h_irttcoef
        #avh_h_albcnv=[]
        #for i in range(6):
        #    avh_h_albcnv.append(int(struct.unpack('<l', fp.read(4))[0]))
        ##print avh_h_albcnv
        #fp.seek(280)
        #avh_h_radtempcnv = np.zeros((3,3))
        #for IRchannelNo in range(3):
        #    for coeffNo in range(3):
        #        avh_h_radtempcnv[IRchannelNo, coeffNo] = \
        #           int(struct.unpack('<l', fp.read(4))[0])
        #print avh_h_radtempcnv
        #IRcalibration['centralWavenumber'] = (avh_h_radtempcnv[:,0] /
        #                                      [1E2, 1E3, 1E3])
        #IRcalibration['c1'] = avh_h_radtempcnv[:,1] / 1E5
        #IRcalibration['c2'] = avh_h_radtempcnv[:,2] / 1E6

        ########################################################
        # Read visible calibration coefficients per scanline
        # - for channels 1, 2, 3A
        ########################################################
        #for scanline in range(1):
        #    avh_calvis=np.zeros((3,3,5))
        #    fp.seek(headerLength + recordLength*scanline + 48)
        #    for VISchannel in range(3):
        #        for sets in range(3):
        #            for coeff in range(5):
        #                avh_calvis[sets, VISchannel, coeff] = \
        #                    int(struct.unpack('<l', fp.read(4))[0])
        #    print avh_calvis
        #    print '----'

        ########################################################
        # Read IR calibration coefficients per scanline
        # - for channels 3B, 4, 5
        ########################################################
        #for scanline in range(1):
        #    avh_calir=np.zeros((2,3,3))
        #    fp.seek(headerLength + recordLength*scanline + 228)
        #    for IRchannelNo in range(3):
        #        for setNo in range(2):
        #            for coeffNo in range(3):
        #                avh_calir[setNo, IRchannelNo, coeffNo] = \
        #                    int(struct.unpack('<l', fp.read(4))[0])

        #avh_filler2 = np.zeros(3)
        #for fillerNo in range(3):
        #    avh_filler2[fillerNo] = int(struct.unpack('<l', fp.read(4))[0])

        #setNo = 0 # Use operational set (the only available)
        #a = np.zeros((3,3))
        #for IRchannelNo in range(3):
        #    for coeffNo in range(3):
        #        # NB: apparently stored "backwards", therefore 2-coeffNo
        #        a[IRchannelNo, 2-coeffNo] = (avh_calir[setNo, IRchannelNo,
        #                                               coeffNo]
        #                                     / np.power(10,
        #                                                avh_filler2[coeffNo]))
        #
        ###########################
        # Apply calibration
        ###########################
        #C = 410
        #for IRchannelNo in range(3):
        #    Ne = a[IRchannelNo,0] + a[IRchannelNo,1]*C + a[IRchannelNo,2]*C*C
        #    vC = IRcalibration['centralWavenumber'][IRchannelNo]
        #    c1 = -IRcalibration['c1'][IRchannelNo] # Note minus
        #    c2 = IRcalibration['c2'][IRchannelNo]
            #print '-----'
            #print a[IRchannelNo,:]
            #print vC, c1, c2
            #TeStar = c2*vC/np.log(1 + (c1*vC*vC*vC)/Ne)
            #Te = c1 + c2*TeStar
            #print Ne, TeStar, Te
            #print '-----'
        #sys.exit('stop')

        ###########################
        # Make Geolocation Arrays
        ###########################
        srcRasterYSize = numCalibratedScanLines

        # Making VRT with raw (unscaled) lon and lat
        # (smaller bands than full dataset)
        self.bandVRTs = {'RawGeolocVRT': VRT(srcRasterXSize=51,
                                            srcRasterYSize=srcRasterYSize)}
        RawGeolocMetaDict = []
        for lonlatNo in range(1, 3):
            RawGeolocMetaDict.append(
                {'src': {'SourceFilename': fileName,
                         'SourceBand': 0,
                         'SourceType': "RawRasterBand",
                         'DataType': gdal.GDT_Int32,
                         'ImageOffset': (headerLength + 640 +
                                         (lonlatNo - 1) * 4),
                         'PixelOffset': 8,
                         'LineOffset': recordLength,
                         'ByteOrder': 'LSB'},
                 'dst': {}})

        self.bandVRTs['RawGeolocVRT']._create_bands(RawGeolocMetaDict)

        # Make derived GeolocVRT with scaled lon and lat
        self.bandVRTs['GeolocVRT'] = VRT(srcRasterXSize=51,
                                        srcRasterYSize=srcRasterYSize)
        GeolocMetaDict = []
        for lonlatNo in range(1, 3):
            GeolocMetaDict.append(
                {'src': {'SourceFilename': (self.bandVRTs['RawGeolocVRT'].
                                            fileName),
                         'SourceBand': lonlatNo,
                         'ScaleRatio': 0.0001,
                         'ScaleOffset': 0,
                         'DataType': gdal.GDT_Int32},
                 'dst': {}})

        self.bandVRTs['GeolocVRT']._create_bands(GeolocMetaDict)

        GeolocObject = GeolocationArray(xVRT=self.bandVRTs['GeolocVRT'],
                                        yVRT=self.bandVRTs['GeolocVRT'],
                                        xBand=2, yBand=1,  # x = lon, y = lat
                                        lineOffset=0, pixelOffset=25,
                                        lineStep=1, pixelStep=40)

        #######################
        # Initialize dataset
        #######################
        # create empty VRT dataset with geolocation only
        # (from Geolocation Array)
        VRT.__init__(self,
                     srcRasterXSize=2048,
                     srcRasterYSize=numCalibratedScanLines,
                     geolocationArray=GeolocObject,
                     srcProjection=GeolocObject.d['SRS'])

        # Since warping quality is horrible using geolocation arrays
        # which are much smaller than raster bands (due to a bug in GDAL:
        # http://trac.osgeo.org/gdal/ticket/4907), the geolocation arrays
        # are here converted to GCPs. Only a subset of GCPs is added,
        # significantly increasing speed when using -tps warping
        reductionFactor = 2
        self.convert_GeolocationArray2GPCs(1*reductionFactor,
                                           40*reductionFactor)

        ##################
        # Create bands
        ##################
        metaDict = []
        ch = ({}, {}, {}, {}, {}, {})

        ch[1]['wavelength'] = 0.63
        ch[2]['wavelength'] = 0.86
        ch[3]['wavelength'] = '1.6 or 3.7 mum'
        ch[4]['wavelength'] = 10.8
        ch[5]['wavelength'] = 12.0

        ch[1]['minmax'] = '0 700'
        ch[2]['minmax'] = '0 700'
        ch[3]['minmax'] = '0 800'
        ch[4]['minmax'] = '400 1000'
        ch[5]['minmax'] = '400 1000'

        for bandNo in range(1, 6):
            metaDict.append({'src': {'SourceFilename': fileName,
                                     'SourceBand': 0,
                                     'SourceType': "RawRasterBand",
                                     'dataType': gdal.GDT_UInt16,
                                     'ImageOffset': imageOffset + (bandNo-1)*2,
                                     'PixelOffset': 10,
                                     'LineOffset': recordLength,
                                     'ByteOrder': 'LSB'},
                            'dst': {'dataType': gdal.GDT_UInt16,
                                    'wkv': 'raw_counts',
                                    'colormap': 'gray',
                                    'wavelength': ch[bandNo]['wavelength'],
                                    'minmax': ch[bandNo]['minmax'],
                                    'unit': "1"}})

        self._create_bands(metaDict)

        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start', time.isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', time.isoformat())

        return
