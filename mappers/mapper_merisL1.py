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

    def __init__(self, rawVRTFileName, fileName, dataset, metadata, vrtBandList):
        ''' Create MER1 VRT '''
        VRT.__init__(self, dataset, metadata, rawVRTFileName);

        product = metadata.get("MPH_PRODUCT", "Not_MERIS")

        if product[0:9] != "MER_FRS_1" and product[0:9] != "MER_RR__1":
            raise AttributeError("MERIS_L1 BAD MAPPER");

        metaDict = [
        {'source': fileName, 'sourceBand':  1, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '412'}},
        {'source': fileName, 'sourceBand':  2, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '443'}},
        {'source': fileName, 'sourceBand':  3, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '490'}},
        {'source': fileName, 'sourceBand':  4, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '510'}},
        {'source': fileName, 'sourceBand':  5, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '560'}},
        {'source': fileName, 'sourceBand':  6, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '620'}},
        {'source': fileName, 'sourceBand':  7, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '665'}},
        {'source': fileName, 'sourceBand':  8, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '680'}},
        {'source': fileName, 'sourceBand':  9, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '708'}},
        {'source': fileName, 'sourceBand': 10, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '753'}},
        {'source': fileName, 'sourceBand': 11, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '761'}},
        {'source': fileName, 'sourceBand': 12, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '778'}},
        {'source': fileName, 'sourceBand': 13, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '864'}},
        {'source': fileName, 'sourceBand': 14, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '849'}},
        {'source': fileName, 'sourceBand': 15, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'parameters': {'wavelength': '900'}},
        {'source': fileName, 'sourceBand': 16, 'wkv': 'quality_flags', 'parameters': {'band_name': 'l1_flags'}}
        ];

        # add 'band_name' to 'parameters'
        for bandDict in metaDict:
            if bandDict['parameters'].has_key('wavelength'):
                bandDict['parameters']['band_name'] = 'radiance_' + bandDict['parameters']['wavelength']

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self._createVRT(metaDict, vrtBandList);
        
        return
