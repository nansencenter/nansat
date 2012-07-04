#-------------------------------------------------------------------------------
# Name:        mapper_hirlam.py
# Purpose:     Mapping for Hirlam model data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:   
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

class Mapper(VRT):
    ''' VRT with mapping of WKV for HIRLAM '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):

        if gdalDataset.GetGeoTransform()[0:5] != (-12.1, 0.2, 0.0, 81.95, 0.0):
            raise AttributeError("HIRLAM BAD MAPPER");

        metaDict = [\
                    {'source': fileName, 'sourceBand': 2, 'wkv': 'eastward_wind', 
                        'parameters':{'band_name': 'east_wind', 'height': '10 m'}, 'NODATA': 9999},
                    {'source': fileName, 'sourceBand': 3, 'wkv': 'northward_wind',
                        'parameters':{'band_name': 'north_wind', 'height': '10 m'}, 'NODATA': 9999},
                    {'source': fileName, 'sourceBand': [2,3], 'wkv': 'wind_speed', 
                        'parameters':{'band_name': 'windspeed', 'height': '10 m', 'pixel_function': 'UVToMagnitude'}},
                    {'source': fileName, 'sourceBand': [2,3], 'wkv': 'wind_from_direction', 
                        'parameters':{'band_name': 'winddirection', 'height': '10 m', 'pixel_function': 'UVToDirectionFrom'}}                    
                    ];
  
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # Create bands
        self._create_bands(metaDict)

        # Adding valid time from the GRIB file to dataset
        validTime = gdalDataset.GetRasterBand(2).GetMetadata()['GRIB_VALID_TIME']
        self._set_time(datetime.datetime.utcfromtimestamp(
                                int(validTime.strip().split(' ')[0])))
        
        return
