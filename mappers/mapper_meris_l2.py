# Name:        nansat_mapper_merisL2
# Purpose:     Mapping for Meris-L2 data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from vrt import VRT
from envisat import Envisat

class Mapper(VRT, Envisat):
    ''' Create VRT with mapping of WKV for MERIS Level 2 (FR or RR)'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create MER2 VRT

        Parameters (**kwargs)
        ---------------------
        geolocation : bool (default True)
            if True, add 'geoloection'
        '''

        product = gdalMetadata["MPH_PRODUCT"]

        if product[0:9] != "MER_FRS_2" and product[0:9] != "MER_RR__2":
            raise AttributeError("MERIS_L2 BAD MAPPER")

        kwDict = {'geolocation' : True}
        # choose kwargs for envisat and asar and change keyname
        for key in kwargs:
            if key.startswith('envisat') or key.startswith('meris'):
                keyName = key.replace('envisat_', '').replace('meris_', '')
                kwDict[keyName] = kwargs[key]
            else:
                kwDict[key] = kwargs[key]

        Envisat.__init__(self, fileName, product[0:4], **kwDict)

        metaDict = [
        {'src': {'SourceFilename': fileName, 'SourceBand':  1}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '412'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  2}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '443'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  3}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '490'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  4}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '510'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  5}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '560'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  6}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '620'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  7}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '665'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  8}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '680'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  9}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '708'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 10}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '753'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 11}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '761'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 12}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '778'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 13}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '864'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 15}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'suffix': '1_log', 'case': 'I'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 16}, 'dst': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'suffix': '2_log', 'case': 'II'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 17}, 'dst': {'wkv': 'mass_concentration_of_suspended_matter_in_sea_water', 'suffix': '2_log', 'case': 'II'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 18}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'suffix': '2_log', 'case': 'II'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 22}, 'dst': {'wkv': 'quality_flags', 'suffix': 'l2'}},
        ]

        # add 'name' to 'parameters'
        for bandDict in metaDict:
            if 'wavelength' in bandDict['dst']:
                bandDict['dst']['suffix'] = bandDict['dst']['wavelength']

        #get GADS from header
        scales = self.read_scaling_gads(range(7, 20) + [20, 21, 22, 20])
        offsets = self.read_scaling_gads(range(33, 46) + [46, 47, 48, 46])
        # set scale/offset to the band metadata (only reflectance)
        for i, bandDict in enumerate(metaDict[:-1]):
            bandDict['src']['ScaleRatio'] = str(scales[i])
            bandDict['src']['ScaleOffset'] = str(offsets[i])

        # add log10-scaled variables
        metaDict += [
            {'src': {'SourceFilename': fileName, 'SourceBand': 1}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'suffix': '1', 'case': 'I', 'expression': 'np.power(10., self["chlor_a_1_log"])'}},
            {'src': {'SourceFilename': fileName, 'SourceBand': 1}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'suffix': '2', 'case': 'II', 'expression': 'np.power(10., self["chlor_a_2_log"])'}},
            {'src': {'SourceFilename': fileName, 'SourceBand': 1}, 'dst': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'suffix': '2', 'case': 'II', 'expression': 'np.power(10., self["cdom_a_2_log"])'}},
            {'src': {'SourceFilename': fileName, 'SourceBand': 1}, 'dst': {'wkv': 'mass_concentration_of_suspended_matter_in_sea_water', 'suffix': '2', 'case': 'II', 'expression': 'np.power(10., self["tsm_2_log"])'}},
        ]

        # get list with resized VRTs from ADS
        self.adsVRTs = []
        self.adsVRTs = self.get_ads_vrts(gdalDataset, ['sun zenith angles', "sun azimuth angles", "zonal winds", "meridional winds"])
        # add bands from the ADS VRTs
        for adsVRT in self.adsVRTs:
            metaDict.append({'src': {'SourceFilename': adsVRT.fileName,
                                     'SourceBand': 1},
                             'dst': {'name':  adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name'),
                                     'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}
                            })

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, **kwDict)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays
        if self.d['geolocation']:
            self.add_geolocation_from_ads(gdalDataset)
