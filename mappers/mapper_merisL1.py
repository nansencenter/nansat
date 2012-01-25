#-------------------------------------------------------------------------------
# Name:        nansat_mapper_merisL1
# Purpose:     Mapping for Meris-L1 data
#
# Author:      antonk
#
# Created:     29.11.2011
# Copyright:   (c) asumak 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

class Mapper(VRT):
    ''' VRT with mapping of WKV for MERIS Level 1 (FR or RR) '''

    def __init__(self, ds, fileName, metadata, vrtBandList, rawVRTName):
        ''' Create MER1 VRT '''
        VRT.__init__(self, metadata, rawVRTName);

        product = metadata.get("MPH_PRODUCT", "Not_MERIS")

        if product[0:9] != "MER_FRS_1" and product[0:9] != "MER_RR__1":
            raise AttributeError("MERIS_L1 BAD MAPPER");

        metaDict = [\
        {'source': fileName, 'sourceBand': 1, 'wkv': 'radiance', 'parameters':{'wavelength': '412'}},\
        {'source': fileName, 'sourceBand': 2, 'wkv': 'radiance', 'parameters':{'wavelength': '443'}},\
        {'source': fileName, 'sourceBand': 3, 'wkv': 'radiance', 'parameters':{'wavelength': '490'}},\
        {'source': fileName, 'sourceBand': 4, 'wkv': 'radiance', 'parameters':{'wavelength': '510'}},\
        {'source': fileName, 'sourceBand': 5, 'wkv': 'radiance', 'parameters':{'wavelength': '560'}},\
        {'source': fileName, 'sourceBand': 6, 'wkv': 'radiance', 'parameters':{'wavelength': '620'}},\
        {'source': fileName, 'sourceBand': 7, 'wkv': 'radiance', 'parameters':{'wavelength': '665'}},\
        {'source': fileName, 'sourceBand': 8, 'wkv': 'radiance', 'parameters':{'wavelength': '680'}},\
        {'source': fileName, 'sourceBand': 9, 'wkv': 'radiance', 'parameters':{'wavelength': '708'}},\
        {'source': fileName, 'sourceBand': 10, 'wkv': 'radiance', 'parameters':{'wavelength': '753'}},\
        {'source': fileName, 'sourceBand': 11, 'wkv': 'radiance', 'parameters':{'wavelength': '761'}},\
        {'source': fileName, 'sourceBand': 12, 'wkv': 'radiance', 'parameters':{'wavelength': '778'}},\
        {'source': fileName, 'sourceBand': 13, 'wkv': 'radiance', 'parameters':{'wavelength': '864'}},\
        {'source': fileName, 'sourceBand': 14, 'wkv': 'radiance', 'parameters':{'wavelength': '849'}},\
        {'source': fileName, 'sourceBand': 15, 'wkv': 'radiance', 'parameters':{'wavelength': '900'}},\
        {'source': fileName, 'sourceBand': 16, 'wkv': 'flags'}
        ];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self.createVRT_and_add_bands(ds, metaDict, vrtBandList);

        return
