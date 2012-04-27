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
                    {'source': fileName, 'sourceBand': 2, 'wkv': 'eastward_wind_velocity', 'parameters':{'band_name': 'east_wind', 'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 3, 'wkv': 'northward_wind_velocity', 'parameters':{'band_name': 'north_wind', 'height': '10 m'}} \
                    ];
  
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._add_all_bands(metaDict)


        ##############################################################
        # Adding derived bands (wind speed and "wind_from_direction") 
        # calculated with pixel functions 
        ##############################################################        
        self._add_pixel_function('UVToMagnitude', [2, 3], fileName, \
                              {'wkv': 'wind_speed', 'parameters': {'band_name': 'speed', 'height': '10 m'}})
        self._add_pixel_function('UVToDirectionFrom', [2, 3], fileName, \
                              {'wkv': 'wind_from_direction', 'parameters': {'band_name': 'direction', 'height': '10 m'}})
        return
