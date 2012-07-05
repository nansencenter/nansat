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
                    {'source': fileName, 'sourceBand': 8, 'wkv': 'eastward_wind', 'parameters':{'band_name': 'east_wind', 'height': '10 m'}},
                    {'source': fileName, 'sourceBand': 9, 'wkv': 'northward_wind', 'parameters':{'band_name': 'north_wind', 'height': '10 m'}},
                    {'source': fileName, 'sourceBand': [8, 9], 'wkv': 'wind_speed', 'parameters':{'pixel_function': 'UVToMagnitude', 'band_name': 'windspeed', 'height': '2 m'}},
                    {'source': fileName, 'sourceBand': [8, 9], 'wkv': 'wind_from_direction', 'parameters':{'pixel_function': 'UVToDirectionFrom', 'band_name': 'winddirection', 'height': '2 m'}},
                    {'source': fileName, 'sourceBand': 6, 'wkv': 'air_temperature', 'parameters':{'band_name': 'air_t', 'height': '2 m'}}
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
