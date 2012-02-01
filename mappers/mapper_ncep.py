#-------------------------------------------------------------------------------
# Name:        nansat_mapper_asar
# Purpose:     Mapping for NCEP GFS model data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:   
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

class Mapper(VRT):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, rawVRTFileName, fileName, dataset, metadata, vrtBandList):
        ''' Create NCEP VRT '''
        VRT.__init__(self, dataset, metadata, rawVRTFileName);

        if dataset.GetGeoTransform() != (-0.25, 0.5, 0.0, 90.25, 0.0, -0.5) or dataset.RasterCount != 9: # Not water proof
            raise AttributeError("NCEP BAD MAPPER");

        metaDict = [\
                    {'source': fileName, 'sourceBand': 8, 'wkv': 'eastward_wind_velocity', 'parameters':{'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 9, 'wkv': 'northward_wind_velocity', 'parameters':{'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 6, 'wkv': 'air_temperature', 'parameters':{'height': '2 m'}}
                    ];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self._createVRT(metaDict, vrtBandList);
        
        ##############################################################
        # Adding derived bands (wind speed and "wind_from_direction") 
        # calculated with pixel functions 
        ##############################################################        
        self._add_pixel_function('UVToMagnitude', [8, 9], fileName, \
                              {'longname': 'wind_speed', 'height': '10 m', 'unit': 'm/s'})
        self._add_pixel_function('UVToDirectionFrom', [8, 9], fileName, \
                              {'longname': 'wind_from_direction', 'height': '10 m', 'unit': 'degrees'})
        return
