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
from vrt import VRT
from envisat import Envisat

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for MERIS Level 1 (FR or RR) '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create MER1 VRT '''
        # get ENVISAT MPH_PRODUCT
        product = gdalMetadata.get("MPH_PRODUCT")
        
        if product[0:9] != "MER_FRS_1" and product[0:9] != "MER_RR__1":
            raise AttributeError("MERIS_L1 BAD MAPPER")
         
        # init ADS parameters
        Envisat.__init__(self, fileName, product[0:4])

        metaDict = [
        {'src': {'SourceFilename': fileName, 'SourceBand':  1}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '412'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  2}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '443'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  3}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '490'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  4}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '510'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  5}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '560'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  6}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '620'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  7}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '665'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  8}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '680'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  9}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '708'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 10}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '753'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 11}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '761'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 12}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '778'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 13}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '864'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 14}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '849'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 15}, 'dst': {'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water', 'wavelength': '900'}},
        {'src': {'SourceFilename': fileName, 'SourceBand': 16, 'DataType': 1}, 'dst': {'wkv': 'quality_flags', 'suffix': 'l1'}}
        ]

        # add DataType into 'src' and suffix into 'dst'
        for bandDict in metaDict:
            if 'DataType' not in bandDict['src']:
                bandDict['src']['DataType'] = 2  # default for meris L1 DataType UInt16
            if bandDict['dst'].has_key('wavelength'):
                bandDict['dst']['suffix'] = bandDict['dst']['wavelength']
        
        # get scaling GADS from header
        scales = self.read_scaling_gads(range(7, 22));
        # set scale factor to the band metadata (only radiances)
        for i, bandDict in enumerate(metaDict[:-1]):
            bandDict['src']['ScaleRatio'] = str(scales[i])

        # get list with resized VRTs from ADS
        self.adsVRTs = self.get_ads_vrts(gdalDataset, ['sun zenith angles', "sun azimuth angles", "zonal winds", "meridional winds"])
        # add bands from the ADS VRTs
        for adsVRT in self.adsVRTs:
            metaDict.append({'src': {'SourceFilename': adsVRT.fileName,
                                     'SourceBand': 1},
                             'dst': {'name':  adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name'),
                                     'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}
                            })

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays
        #self.add_geolocation_from_ads(gdalDataset, step=1)
