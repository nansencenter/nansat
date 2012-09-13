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
from vrt import VRT, Geolocation
from envisat import Envisat

class Mapper(VRT, Envisat):
    ''' Create VRT with mapping of WKV for MERIS Level 2 (FR or RR) '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        product = gdalMetadata.get("MPH_PRODUCT", "Not_MERIS")

        if product[0:9] != "MER_FRS_2" and product[0:9] != "MER_RR__2":
            raise AttributeError("MERIS_L2 BAD MAPPER")

        # Create VRTdataset with small VRTRawRasterbands
        #geoDataset = self.create_VRT_with_rawbands(fileName, product[0:4], ["latitude", "zonal winds"])
        #
        # Enlarge the band to the underlying data band size
        #self.geoDataset = geoDataset.resized(gdalDataset.RasterXSize, gdalDataset.RasterYSize)

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
        {'src': {'SourceFilename': fileName, 'SourceBand': 15}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'BandName': 'algal_1', 'case': 'I'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 16}, 'dst': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'BandName': 'yellow_subs', 'case': 'II'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 17}, 'dst': {'wkv': 'mass_concentration_of_suspended_matter_in_sea_water', 'BandName': 'total_susp', 'case': 'II'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 18}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'BandName': 'algal_2', 'case': 'II'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 22}, 'dst': {'wkv': 'quality_flags', 'BandName': 'l2_flags'}},
        ]

        # add 'BandName' to 'parameters'
        for bandDict in metaDict:
            if 'wavelength' in bandDict['dst']:
                bandDict['dst']['BandName'] = 'reflectance_' + bandDict['dst']['wavelength']

        #get GADS from header
        scales = self.read_scaling_gads(fileName, range(7, 20) + [20, 21, 22, 20])
        offsets = self.read_scaling_gads(fileName, range(33, 46) + [46, 47, 48, 46])
        # set scale/offset to the band metadata (only reflectance)
        for i, bandDict in enumerate(metaDict[:-1]):
            bandDict['src']['ScaleRatio'] = str(scales[i])
            bandDict['src']['ScaleOffset'] = str(offsets[i])

        #add geolocation dictionary into metaDict
        #for iBand in range(self.geoDataset.dataset.RasterCount):
        #    bandMetadata = self.geoDataset.dataset.GetRasterBand(iBand+1).GetMetadata()
        #    metaDict.append({'src': {'SourceFilename': self.geoDataset.fileName, 'SourceBand': iBand+1}, 'dst': {'wkv': '', 'parameters':bandMetadata})

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        ''' Set GeolocationArray '''
        #latlonName = {"latitude":"latitude","longitude":"longitude"}
        #self.add_geoarray_dataset(fileName, product[0:4], gdalDataset.RasterXSize, gdalDataset.RasterYSize, latlonName, gdalDataset.GetGCPProjection())
