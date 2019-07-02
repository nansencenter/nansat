# Name:        nansat_mapper_merisL2
# Purpose:     Mapping for Meris-L2 data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse

from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError
from nansat.mappers.envisat import Envisat


class Mapper(VRT, Envisat):
    ''' Create VRT with mapping of WKV for MERIS Level 2 (FR or RR)'''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 geolocation=False, zoomSize=500, step=1, **kwargs):

        ''' Create MER2 VRT

        Parameters
        -----------
        filename : string
        gdalDataset : gdal dataset
        gdalMetadata : gdal metadata
        geolocation : bool (default is False)
            if True, add gdal geolocation
        zoomSize: int (used in envisat.py)
            size, to which the ADS array will be zoomed using scipy
            array of this size will be stored in memory
        step: int (used in envisat.py)
            step of pixel and line in GeolocationArrays. lat/lon grids are
            generated at that step
        '''

        self.setup_ads_parameters(filename, gdalMetadata)

        if self.product[0:9] != "MER_FRS_2" and self.product[0:9] != "MER_RR__2":
            raise WrongMapperError(filename)

        metaDict = [{'src': {'SourceFilename': filename, 'SourceBand': 1},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '412'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 2},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '443'}},
                    {'src': {'SourceFilename': filename, 'SourceBand':  3},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '490'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 4},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '510'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 5},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '560'}},
                    {'src': {'SourceFilename': filename, 'SourceBand':  6},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '620'}},
                    {'src': {'SourceFilename': filename, 'SourceBand':  7},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '665'}},
                    {'src': {'SourceFilename': filename, 'SourceBand':  8},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '680'}},
                    {'src': {'SourceFilename': filename, 'SourceBand':  9},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '708'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 10},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '753'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 11},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '761'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 12},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '778'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 13},
                     'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                             'wavelength': '864'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 15},
                     'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                             'suffix': '1_log', 'case': 'I'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 16},
                     'dst': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter',
                             'suffix': '2_log', 'case': 'II'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 17},
                     'dst': {'wkv': 'mass_concentration_of_suspended_matter_in_sea_water',
                             'suffix': '2_log', 'case': 'II'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 18},
                     'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                             'suffix': '2_log', 'case': 'II'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 22},
                     'dst': {'wkv': 'quality_flags', 'suffix': 'l2'}}
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
        metaDict += [{'src': {'SourceFilename': filename, 'SourceBand': 1},
                      'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                              'suffix': '1', 'case': 'I',
                              'expression': 'np.power(10., self["chlor_a_1_log"])'}},
                     {'src': {'SourceFilename': filename, 'SourceBand': 1},
                      'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                              'suffix': '2', 'case': 'II',
                              'expression': 'np.power(10., self["chlor_a_2_log"])'}},
                     {'src': {'SourceFilename': filename, 'SourceBand': 1},
                      'dst': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter',
                              'suffix': '2', 'case': 'II',
                              'expression': 'np.power(10., self["cdom_a_2_log"])'}},
                     {'src': {'SourceFilename': filename, 'SourceBand': 1},
                      'dst': {'wkv': 'mass_concentration_of_suspended_matter_in_sea_water',
                              'suffix': '2', 'case': 'II',
                              'expression': 'np.power(10., self["tsm_2_log"])'}}
                     ]

        # get list with resized VRTs from ADS
        self.band_vrts = {'adsVRTs': self.get_ads_vrts(gdalDataset,
                                                     ['sun zenith angles',
                                                      'sun azimuth angles',
                                                      'zonal winds',
                                                      'meridional winds'],
                                                     zoomSize=zoomSize,
                                                     step=step)}

        # add bands from the ADS VRTs
        for adsVRT in self.band_vrts['adsVRTs']:
            metaDict.append({'src': {'SourceFilename': adsVRT.filename,
                                     'SourceBand': 1},
                             'dst': {'name': (adsVRT.dataset.GetRasterBand(1).
                                              GetMetadataItem('name')),
                                     'units': (adsVRT.dataset.GetRasterBand(1).
                                               GetMetadataItem('units'))}
                             })

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays
        if geolocation:
            self.add_geolocation_from_ads(gdalDataset,
                                          zoomSize=zoomSize, step=step)
        # set time
        self._set_envisat_time(gdalMetadata)

        self.dataset.SetMetadataItem('sensor', 'MERIS')
        self.dataset.SetMetadataItem('satellite', 'ENVISAT')
