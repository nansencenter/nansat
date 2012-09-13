#-------------------------------------------------------------------------------
# Name:        nansat_mapper_asar
# Purpose:     Mapping for NCEP GFS model data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:   
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from vrt import VRT, datetime

class Mapper(VRT):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create NCEP VRT '''

        if gdalDataset.GetGeoTransform() != (-0.25, 0.5, 0.0, 90.25, 0.0, -0.5) or gdalDataset.RasterCount != 9: # Not water proof
            raise AttributeError("NCEP BAD MAPPER");

        metaDict = [
                    {'src': {'SourceFilename': fileName, 'SourceBand': 8}, 'dst': {'wkv': 'eastward_wind', 'BandName': 'east_wind', 'height': '10 m'}},
                    {'src': {'SourceFilename': fileName, 'SourceBand': 9}, 'dst': {'wkv': 'northward_wind', 'BandName': 'north_wind', 'height': '10 m'}},
                    {'src': [{'SourceFilename': fileName, 'SourceBand': 8}, {'SourceFilename': fileName, 'SourceBand': 9}], 'dst': {'wkv': 'wind_speed', 'PixelFunctionType': 'UVToMagnitude', 'BandName': 'windspeed', 'height': '2 m'}},
                    {'src': [{'SourceFilename': fileName, 'SourceBand': 8}, {'SourceFilename': fileName, 'SourceBand': 9}], 'dst': {'wkv': 'wind_from_direction', 'PixelFunctionType': 'UVToDirectionFrom', 'BandName': 'winddirection', 'height': '2 m'}},
                    {'src': {'SourceFilename': fileName, 'SourceBand': 6}, 'dst': {'wkv': 'air_temperature', 'BandName': 'air_t', 'height': '2 m'}}
                    ];

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset);
            
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
        
        # Adding valid time from the GRIB file to dataset
        validTime = gdalDataset.GetRasterBand(8).GetMetadata()['GRIB_VALID_TIME']
        self._set_time(datetime.datetime.utcfromtimestamp(
                                int(validTime.strip().split(' ')[0])))
        
        return
