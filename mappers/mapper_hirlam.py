#-------------------------------------------------------------------------------
# Name:        nansat_mapper_asar
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

    def __init__(self, ds, fileName, metadata, vrtBandList, rawVRTName):
        ''' Create HIRLAM VRT '''
        VRT.__init__(self, metadata, rawVRTName);

        if ds.GetGeoTransform()[0:5] != (-12.1, 0.2, 0.0, 81.95, 0.0):
            raise AttributeError("HIRLAM BAD MAPPER");

        metaDict = [\
                    {'source': fileName, 'sourceBand': 2, 'wkv': 'eastward_wind_velocity', 'parameters':{'height': '10 m'}}, \
                    {'source': fileName, 'sourceBand': 3, 'wkv': 'northward_wind_velocity', 'parameters':{'height': '10 m'}} \
                    ];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self.createVRT_and_add_bands(ds, metaDict, vrtBandList);

        return
