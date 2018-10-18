# Name:         mapper_hirlam.py
# Purpose:      Nansat mapping for Hirlam model data
#               (GRIB files from www.yr.no)
# Authors:      Knut-Frode Dagestad, Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import datetime
import json

import numpy

from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError

import pythesint as pti


class Mapper(VRT):
    ''' VRT with mapping of WKV for HIRLAM '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):

        try:
            geo_transform = gdalDataset.GetGeoTransform()[0:5]
        except AttributeError:
            raise WrongMapperError
        if geo_transform != (-12.1, 0.2, 0.0, 81.95, 0.0):
            raise WrongMapperError

        metaDict = [{'src': {'SourceFilename': filename,
                             'SourceBand': 2,
                             'NODATA': 9999},
                     'dst': {'wkv': 'eastward_wind',
                             'height': '10 m'}
                     },
                    {'src': {'SourceFilename': filename,
                             'SourceBand': 3,
                             'NODATA': 9999},
                     'dst': {'wkv': 'northward_wind',
                             'height': '10 m'}
                     },
                    {'src': [{'SourceFilename': filename,
                              'SourceBand': 2,
                              'DataType': gdalDataset.GetRasterBand(2).DataType
                              },
                             {'SourceFilename': filename,
                              'SourceBand': 3,
                              'DataType': gdalDataset.GetRasterBand(3).DataType
                              }],
                     'dst': {'wkv': 'wind_speed',
                             'name': 'windspeed',
                             'height': '10 m',
                             'PixelFunctionType': 'UVToMagnitude',
                             'NODATA': 9999}
                     },
                    {'src': [{'SourceFilename': filename,
                              'SourceBand': 2,
                              'DataType': gdalDataset.GetRasterBand(2).DataType
                              },
                             {'SourceFilename': filename,
                              'SourceBand': 3,
                              'DataType': gdalDataset.GetRasterBand(3).DataType
                              }],
                     'dst': {'wkv': 'wind_from_direction',
                             'name': 'winddirection',
                             'height': '10 m',
                             'PixelFunctionType': 'UVToDirectionFrom',
                             'NODATA': 9999
                             }
                     }]

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalDataset, metadata=gdalMetadata)

        # Create bands
        self.create_bands(metaDict)

        # set source, start_date, stop_date
        self.dataset.SetMetadataItem('source', 'HIRLAM')

        # Adding valid time from the GRIB file to dataset
        start_date = gdalDataset.GetRasterBand(1).GetMetadata()['GRIB_VALID_TIME']
        self.dataset.SetMetadataItem('time_coverage_start',
                datetime.datetime.utcfromtimestamp(
                int(start_date.strip().split(' ')[0])).isoformat() + '+00:00')

        stop_date = gdalDataset.GetRasterBand(gdalDataset.RasterCount).GetMetadata()['GRIB_VALID_TIME']
        self.dataset.SetMetadataItem('time_coverage_end',
                datetime.datetime.utcfromtimestamp(
                int(stop_date.strip().split(' ')[0])).isoformat() + '+00:00')

        mm = pti.get_gcmd_instrument('computer')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        ee = pti.get_gcmd_platform('merged analysis')
        self.dataset.SetMetadataItem('platform', json.dump(ee))
