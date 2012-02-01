#-------------------------------------------------------------------------------
# Name:        nansat_mapper_merisL2
# Purpose:     Mapping for Meris-L2 data
#
# Author:      antonk
#
# Created:     29.11.2011
# Copyright:   (c) asumak 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

class Mapper(VRT):
    ''' Create VRT with mapping of WKV for MERIS Level 2 (FR or RR) '''

    def __init__(self, rawVRTFileName, fileName, dataset, metadata, vrtBandList):
        ''' Create MER2 VRT '''
        VRT.__init__(self, dataset, metadata, rawVRTFileName);

        product = metadata.get("MPH_PRODUCT", "Not_MERIS")

        if product[0:9] != "MER_FRS_2" and product[0:9] != "MER_RR__2":
            raise AttributeError("MERIS_L2 BAD MAPPER");

        metaDict = [\
        {'source': fileName, 'sourceBand': 1,  'wkv': 'reflectance', 'parameters': {'wavelength': '412'} },\
        {'source': fileName, 'sourceBand': 2, 'wkv': 'reflectance', 'parameters':{'wavelength': '443'}},\
        {'source': fileName, 'sourceBand': 3, 'wkv': 'reflectance', 'parameters':{'wavelength': '490'}},\
        {'source': fileName, 'sourceBand': 4, 'wkv': 'reflectance', 'parameters':{'wavelength': '510'}},\
        {'source': fileName, 'sourceBand': 5, 'wkv': 'reflectance', 'parameters':{'wavelength': '560'}},\
        {'source': fileName, 'sourceBand': 6, 'wkv': 'reflectance', 'parameters':{'wavelength': '620'}},\
        {'source': fileName, 'sourceBand': 7, 'wkv': 'reflectance', 'parameters':{'wavelength': '665'}},\
        {'source': fileName, 'sourceBand': 8, 'wkv': 'reflectance', 'parameters':{'wavelength': '680'}},\
        {'source': fileName, 'sourceBand': 9, 'wkv': 'reflectance', 'parameters':{'wavelength': '708'}},\
        {'source': fileName, 'sourceBand': 10, 'wkv': 'reflectance', 'parameters':{'wavelength': '753'}},\
        {'source': fileName, 'sourceBand': 11, 'wkv': 'reflectance', 'parameters':{'wavelength': '761'}},\
        {'source': fileName, 'sourceBand': 12, 'wkv': 'reflectance', 'parameters':{'wavelength': '778'}},\
        {'source': fileName, 'sourceBand': 13, 'wkv': 'reflectance', 'parameters':{'wavelength': '864'}},\
        {'source': fileName, 'sourceBand': 15,  'wkv': 'chlor_a', 'parameters': {'case': 'I'} },\
        {'source': fileName, 'sourceBand': 16,  'wkv': 'a_doc', 'parameters': {'case': 'II'} },\
        {'source': fileName, 'sourceBand': 17,  'wkv': 'spm', 'parameters': {'case': 'II'} },\
        {'source': fileName, 'sourceBand': 18,  'wkv': 'chlor_a', 'parameters': {'case': 'II'} },\
        {'source': fileName, 'sourceBand': 22, 'wkv': 'flags'},\
        ];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);

        self._createVRT(metaDict, vrtBandList);
        
        return
