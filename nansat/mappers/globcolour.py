# Name:         globcolour
# Purpose:      Data and methods shared by GLOBCOLOUR mappers
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from copy import deepcopy


class Globcolour():
    ''' Mapper for GLOBCOLOR L3M products'''

    # detect wkv from metadata 'Parameter'
    varname2wkv = {'CHL1_mean': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                   'CHL2_mean': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                   'KD490_mean': 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water',
                   'L412_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L443_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L490_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L510_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L531_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L555_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L620_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L670_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L681_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'L709_mean': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water',
                   'CDM_mean': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter',
                   'BBP_mean': 'volume_backscattering_coefficient_of_radiative_flux_in_sea_water_due_to_suspended_particles',
                   'PAR_mean': 'surface_downwelling_photosynthetic_radiative_flux_in_air',
                       }

    def make_rrsw_meta_entry(self, nlwMetaEntry):
        '''Make metaEntry for calculation of Rrsw'''
        iWKV = nlwMetaEntry['dst']['wkv']
        if 'solar_irradiance' not in nlwMetaEntry['dst']:
            return None
        if iWKV == 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water':
            solarIrradiance = nlwMetaEntry['dst']['solar_irradiance']
            wavelength = nlwMetaEntry['dst']['wavelength']
            metaEntry = deepcopy(nlwMetaEntry)
            metaEntry['dst']['wkv'] = 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_water'
            metaEntry['dst']['expression'] = ('self["nLw_%s"] / %s'
                                              % (wavelength, solarIrradiance))
        else:
            metaEntry = None

        return metaEntry
