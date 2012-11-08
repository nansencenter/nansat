#-------------------------------------------------------------------------------
# Name:        mapper_hirlam.py
# Purpose:     Mapping for Hirlam model data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:   
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

class Mapper(VRT):
    ''' VRT with mapping of WKV for HIRLAM '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        if gdalDataset.GetGeoTransform()[0:5] != (-12.1, 0.2, 0.0, 81.95, 0.0):
            raise AttributeError("HIRLAM BAD MAPPER");

        metaDict = [
                    {'src': {'SourceFilename': fileName, 'SourceBand': 2, 'NODATA': 9999},
                     'dst': {'wkv': 'eastward_wind', 'name': 'east_wind', 'height': '10 m'}},
                    {'src': {'SourceFilename': fileName, 'SourceBand': 3, 'NODATA': 9999},
                     'dst': {'wkv': 'northward_wind', 'name': 'north_wind', 'height': '10 m'}},
                    {'src': [{'SourceFilename': fileName, 'SourceBand': 2}, {'SourceFilename': fileName, 'SourceBand': 3}],
                     'dst': {'wkv': 'wind_speed', 'name': 'windspeed', 'height': '10 m', 'PixelFunctionType': 'UVToMagnitude'}},
                    {'src': [{'SourceFilename': fileName, 'SourceBand': 2}, {'SourceFilename': fileName, 'SourceBand': 2}],
                     'dst': {'wkv': 'wind_from_direction', 'name': 'winddirection', 'height': '10 m', 'PixelFunctionType': 'UVToDirectionFrom'}}                    
                    ];
  
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset);

        # Create bands
        self._create_bands(metaDict)

        # Adding valid time from the GRIB file to dataset
        validTime = gdalDataset.GetRasterBand(2).GetMetadata()['GRIB_VALID_TIME']
        self._set_time(datetime.datetime.utcfromtimestamp(
                                int(validTime.strip().split(' ')[0])))
        
        return
