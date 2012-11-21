#-------------------------------------------------------------------------------
# Name:        mapper_aapp.py
# Purpose:     Mapping for AAPP output
#
# Author:      Knut-Frode
#
# Created:     08.11.2012
# Copyright:   
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Description of file format:
# http://research.metoffice.gov.uk/research/interproj/nwpsaf/aapp/NWPSAF-MF-UD-003_Formats.pdf (page 8-)

import sys
import struct
import datetime
from vrt import *

satIDs = {4: 'NOAA-15',2: 'NOAA-16', 6: 'NOAA-17',7: 'NOAA-18', 8: 'NOAA-19', 
          11: 'Metop-B (Metop-1)', 12: 'Metop-A (Metop-2)', 13: 'Metop-C (Metop-3)'}
dataFormats = {1: 'LAC', 2: 'GAC', 3: 'HRPT'}
dataSetQualityIndicatorOffset = 114
recordLength = 22016
headerLength = recordLength
imageOffset = headerLength + 1262


class Mapper(VRT):
    ''' VRT with mapping of WKV for AVHRR L1B output from AAPP '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
    
        ########################################
        # Read metadata from binary file
        ########################################
        fp = open(fileName, 'rb')
        fp.seek(72)
        satNum = int(struct.unpack('<H', fp.read(2))[0])
        satID = satIDs[satNum]
        fp.seek(76)
        dataFormatNum = int(struct.unpack('<H', fp.read(2))[0])
        dataFormat = dataFormats[dataFormatNum]
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
        print year
        dayofyear = int(struct.unpack('<H', fp.read(2))[0])
        millisecondsOfDay = int(struct.unpack('<l', fp.read(4))[0])
        time = datetime.datetime(year,1,1) + datetime.timedelta(dayofyear-1, 
                        milliseconds=millisecondsOfDay)

        ###########################
        # Make Geolocation Arrays
        ###########################

        # Making VRT with raw (unscaled) lon and lat (smaller bands than full dataset)
        RawGeolocVRT = VRT(srcRasterXSize=51, srcRasterYSize=numCalibratedScanLines)
        RawGeolocMetaDict = []
        for lonlatNo in range(1,3):
            RawGeolocMetaDict.append({'src': {
                        'SourceFilename': fileName,
                        'SourceBand': 0,
                        'SourceType': "RawRasterBand",
                        'DataType': gdal.GDT_Int32,
                        'ImageOffset' : headerLength + 640 + (lonlatNo-1)*4,
                        'PixelOffset' : 8,
                        'LineOffset' : recordLength,
                        'ByteOrder' : 'LSB'},
                    'dst': {
                        #'dataType': gdal.GDT_Int32,
                        'wkv': 'latitude'}})

        RawGeolocVRT._create_bands(RawGeolocMetaDict)

#       # Make derived GeolocVRT with scaled lon and lat
        GeolocVRT = VRT(srcRasterXSize=51, srcRasterYSize=numCalibratedScanLines)
        GeolocMetaDict = []
        for lonlatNo in range(1,3):
            GeolocMetaDict.append({'src': {
                        'SourceFilename': RawGeolocVRT.fileName,
                        'SourceBand': lonlatNo,
                        'ScaleRatio': 0.0001,
                        'ScaleOffset': 0,
                        'DataType': gdal.GDT_Int32},
                    'dst': {
                        #'dataType': gdal.GDT_Int32,
                        'wkv': 'latitude'}})

        GeolocVRT._create_bands(GeolocMetaDict)

        GeolocObject = Geolocation(xVRT=GeolocVRT, yVRT=GeolocVRT,
                    xBand=2, yBand=1, # x = lon, y = lat
                    lineOffset=0, pixelOffset=25, lineStep=1, pixelStep=40)

        #######################
        # Initialize dataset
        #######################
        # create empty VRT dataset with geolocation only (from Geolocation Array)
        VRT.__init__(self, srcRasterXSize=2048, srcRasterYSize=numCalibratedScanLines, 
                        geolocation=GeolocObject, srcProjection=GeolocObject.d['SRS'])

        ##################
        # Create bands
        ##################
        metaDict = []

        for bandNo in range(1,6):
            metaDict.append({'src': {
                                'SourceFilename': fileName,
                                'SourceBand': 0,
                                'SourceType': "RawRasterBand",
                                'dataType': gdal.GDT_UInt16,
                                'ImageOffset' : imageOffset + (bandNo-1)*2,
                                'PixelOffset' : 10,
                                'LineOffset' : recordLength,
                                'ByteOrder' : 'LSB'},
                            'dst': {
                                'dataType': gdal.GDT_UInt16,
                                'wkv': 'latitude',
                                'minmax': '500 1100',
                                "unit": "1"}})

        self._create_bands(metaDict)

        # Adding valid time to dataset
        self._set_time(time)
        
        return
