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

    def __init__(self, rawVRTFileName, fileName, dataset, metadata, vrtBandList):
        ''' Create HIRLAM VRT '''
        VRT.__init__(self, dataset, metadata, rawVRTFileName);

        if dataset.GetGeoTransform()[0:5] != (-12.1, 0.2, 0.0, 81.95, 0.0):
            raise AttributeError("HIRLAM BAD MAPPER");

        metaDict = [\
                    {'source': fileName, 'sourceBand': 2, 'wkv': 'eastward_wind_velocity', 'parameters':{'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 3, 'wkv': 'northward_wind_velocity', 'parameters':{'height': '10 m'}} \
                    ];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self._createVRT(metaDict, vrtBandList);

        ##############################################################
        # Adding derived bands (wind speed and "wind_from_direction") 
        # calculated with pixel functions 
        ##############################################################        
        self._add_pixel_function('UVToMagnitude', [2, 3], fileName, \
                              {'longname': 'wind_speed', 'height': '10 m', 'unit': 'm/s'})
        self._add_pixel_function('UVToDirectionFrom', [2, 3], fileName, \
                              {'longname': 'wind_from_direction', 'height': '10 m', 'unit': 'degrees'})
        return
