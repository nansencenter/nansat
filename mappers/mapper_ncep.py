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

    def __init__(self, ds, fileName, metadata, vrtBandList, rawVRTName):
        ''' Create NCEP VRT '''
        VRT.__init__(self, metadata, rawVRTName);

        if ds.GetGeoTransform() != (-0.25, 0.5, 0.0, 90.25, 0.0, -0.5) or ds.RasterCount != 9: # Not water proof
            raise AttributeError("NCEP BAD MAPPER");

        metaDict = [\
                    {'source': fileName, 'sourceBand': 8, 'wkv': 'eastward_wind_velocity', 'parameters':{'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 9, 'wkv': 'northward_wind_velocity', 'parameters':{'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 6, 'wkv': 'air_temperature', 'parameters':{'height': '2 m'}}
                    ];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self.createVRT_and_add_bands(ds, metaDict, vrtBandList);
        
        ##############################################################
        # Adding derived band (wind speed) calculated
        # using pixel function "UVToMagnitude":
        #      wind_speed = sqrt(u*2+v**2) 
        ##############################################################        
        self.addPixelFunction('UVToMagnitude', [8, 9], fileName, \
                              {'longname': 'wind_speed', 'height': '10 m', 'unit': 'm/s'})
        self.addPixelFunction('UVToDirectionFrom', [8, 9], fileName, \
                              {'longname': 'wind_from_direction', 'height': '10 m', 'unit': 'm/s'})
        return
