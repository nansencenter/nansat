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
        {'source': fileName, 'sourceBand':  1, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '412'} },\
        {'source': fileName, 'sourceBand':  2, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '443'}},\
        {'source': fileName, 'sourceBand':  3, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '490'}},\
        {'source': fileName, 'sourceBand':  4, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '510'}},\
        {'source': fileName, 'sourceBand':  5, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '560'}},\
        {'source': fileName, 'sourceBand':  6, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '620'}},\
        {'source': fileName, 'sourceBand':  7, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '665'}},\
        {'source': fileName, 'sourceBand':  8, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '680'}},\
        {'source': fileName, 'sourceBand':  9, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '708'}},\
        {'source': fileName, 'sourceBand': 10, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '753'}},\
        {'source': fileName, 'sourceBand': 11, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '761'}},\
        {'source': fileName, 'sourceBand': 12, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '778'}},\
        {'source': fileName, 'sourceBand': 13, 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '864'}},\
        {'source': fileName, 'sourceBand': 15, 'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'parameters': {'band_name': 'algal_1', 'case': 'I'} },\
        {'source': fileName, 'sourceBand': 16, 'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'parameters': {'band_name': 'yellow_subs', 'case': 'II'} },\
        {'source': fileName, 'sourceBand': 17, 'wkv': 'mass_concentration_of_suspended_matter_in_sea_water', 'parameters': {'band_name': 'total_susp', 'case': 'II'} },\
        {'source': fileName, 'sourceBand': 18, 'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'parameters': {'band_name': 'algal_2', 'case': 'II'} },\
        {'source': fileName, 'sourceBand': 22, 'wkv': 'quality_flags', 'parameters': {'band_name': 'l2_flags'} },\
        ];

        # add 'band_name' to 'parameters'
        for bandDict in metaDict:
            if bandDict['parameters'].has_key('wavelength'):
                bandDict['parameters']['band_name'] = 'reflectance_' + bandDict['parameters']['wavelength']

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);

        self._createVRT(metaDict, vrtBandList);
        
        return
