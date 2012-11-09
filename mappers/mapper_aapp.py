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

import struct
import datetime
from vrt import *

satIDs = {4: 'NOAA-15',2: 'NOAA-16', 6: 'NOAA-17',7: 'NOAA-18', 8: 'NOAA-19', 
          11: 'Metop-B (Metop-1)', 12: 'Metop-A (Metop-2)', 13: 'Metop-C (Metop-3)'}
dataFormats = {1: 'LAC', 2: 'GAC', 3: 'HRPT'}
dataSetQualityIndicatorOffset = 114

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

        ##################
        # Read time
        ##################
        fp.seek(84)
        year = int(struct.unpack('<H', fp.read(2))[0])
        dayofyear = int(struct.unpack('<H', fp.read(2))[0])
        millisecondsOfDay = int(struct.unpack('<L', fp.read(4))[0])
        time = datetime.datetime(year,1,1) + datetime.timedelta(dayofyear-1, 
                        milliseconds=millisecondsOfDay)
        fp.close()

        #######################
        # Initialize dataset
        #######################
        # To be updated
        projection = 'GEOGCS["Coordinate System imported from GRIB file",DATUM["unknown",SPHEROID["Sphere",6367470,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, srcRasterXSize=2048, srcRasterYSize=numCalibratedScanLines, 
                        srcProjection=projection, 
                        srcGeoTransform=(-12.1, 0.2, 0.0, 81.95, 0.0, -0.1))

        ##################
        # Create bands
        ##################
        recordLength = 22016
        imageOffset = recordLength + 1262
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
                                'dataType': gdal.GDT_Float32,
                                'wkv': 'albedo',
                                'minmax': '500 1100',
                                "unit": "1"}})

        self._create_bands(metaDict)

        # Adding valid time to dataset
        self._set_time(time)
        
        return
