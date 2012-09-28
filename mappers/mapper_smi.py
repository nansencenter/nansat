#-------------------------------------------------------------------------------
# Name:        mapper_smi
# Purpose:     Mapping for HMODISA Level-3 Standard Mapped Image from OBPG
#
# Author:      antonk
#
# Created:     29.11.2011
# Copyright:   (c) NERSC 2012
#-------------------------------------------------------------------------------
import datetime

import numpy as np

from vrt import VRT, Geolocation

class Mapper(VRT):
    ''' Mapper for HMODISA Level-3 Standard Mapped Image from OBPG'''
    
    # detect wkv from metadata 'Parameter'
    wkvs = {
    'Chlorophyll a concentration': 'mass_concentration_of_chlorophyll_a_in_sea_water'
    }
    
    bandNames = {
    'Chlorophyll a concentration': 'algal_1'
    }

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' SMI VRT '''

        print "=>%s<=" % gdalMetadata['Title']

        if gdalMetadata['Title'] != 'HMODISA Level-3 Standard Mapped Image':
            raise AttributeError("SMI BAD MAPPER")
        
        # get wkv
        wkv = self.wkvs[gdalMetadata['Parameter']]
        
        #get array with data and make 'mask'
        a = gdalDataset.ReadAsArray()
        mask = np.zeros(a.shape, 'uint8') + 128
        mask[a < -32000] = 1
        self.maskVRT = VRT(array=mask)

        metaDict = [
        {'src': {'SourceFilename': fileName,
                 'SourceBand':  1,
                 'ScaleRatio': float(gdalMetadata['Slope']),
                 'ScaleOffset': float(gdalMetadata['Intercept'])},
         'dst': {'wkv': wkv,
                 'BandName': self.bandNames[gdalMetadata['Parameter']]}},
        {'src': {'SourceFilename': self.maskVRT.fileName, 'SourceBand':  1},
         'dst': {'BandName': 'mask'}}
        ]

        # create empty VRT dataset with geolocation only
        latitudeStep = float(gdalMetadata['Latitude Step'])
        longitudeStep = float(gdalMetadata['Longitude Step'])
        VRT.__init__(self, srcGeoTransform=(-180.0, longitudeStep, 0.0, 90.0, 0.0, -longitudeStep),
                           srcProjection='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]',
                           srcRasterXSize=int(gdalMetadata['Number of Columns']),
                           srcRasterYSize=int(gdalMetadata['Number of Lines'])
                    )
                           

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # Add valid time
        startYear = int(gdalMetadata['Start Year'])
        startDay = int(gdalMetadata['Start Day'])
        self._set_time(datetime.datetime(startYear, 1, 1) + datetime.timedelta(startDay))
