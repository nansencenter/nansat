# Name:         mapper_ncep
# Purpose:      Nansat mapping for NCEP GFS model data
# Author:       Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# Made for GRIB files downloaded from http://nomads.ncep.noaa.gov/
# NB: Band numbers is hardcoded for band subsets extracted at NERSC,
# mapper will not work for other NCEP GFS files before made more generic

from nansat.vrt import VRT, datetime


class Mapper(VRT):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create NCEP VRT '''

        if (gdalDataset.GetGeoTransform() != (-0.25, 0.5, 0.0,
                                              90.25, 0.0, -0.5) or
                gdalDataset.RasterCount != 9):  # Not water proof
            raise AttributeError("NCEP BAD MAPPER")

        metaDict = [{'src': {'SourceFilename': fileName,
                             'SourceBand': 8},
                     'dst': {'wkv': 'eastward_wind',
                             'height': '10 m'}},
                    {'src': {'SourceFilename': fileName,
                             'SourceBand': 9},
                     'dst': {'wkv': 'northward_wind',
                             'height': '10 m'}},
                    {'src': [{'SourceFilename': fileName,
                              'SourceBand': 8,
                              'DataType': gdalDataset.GetRasterBand(8).DataType
                              },
                             {'SourceFilename': fileName,
                              'SourceBand': 9,
                              'DataType': gdalDataset.GetRasterBand(9).DataType
                              }],
                     'dst': {'wkv': 'wind_speed',
                             'PixelFunctionType': 'UVToMagnitude',
                             'name': 'windspeed',
                             'height': '2 m'
                             }},
                    {'src': [{'SourceFilename': fileName,
                              'SourceBand': 8,
                              'DataType': gdalDataset.GetRasterBand(8).DataType
                              },
                             {'SourceFilename': fileName,
                              'SourceBand': 9,
                              'DataType': gdalDataset.GetRasterBand(9).DataType
                              }],
                     'dst': {'wkv': 'wind_from_direction',
                             'PixelFunctionType': 'UVToDirectionFrom',
                             'name': 'winddirection',
                             'height': '2 m'
                             }},
                    {'src': {'SourceFilename': fileName,
                             'SourceBand': 6},
                     'dst': {'wkv': 'air_temperature',
                             'name': 'air_t',
                             'height': '2 m'}
                     }]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # Adding valid time from the GRIB file to dataset
        band = gdalDataset.GetRasterBand(8)
        validTime = band.GetMetadata()['GRIB_VALID_TIME']
        self._set_time(datetime.datetime.
                       utcfromtimestamp(int(validTime.strip().split(' ')[0])))

        return
