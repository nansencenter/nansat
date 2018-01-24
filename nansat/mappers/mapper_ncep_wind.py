# Name:         mapper_ncep_wind
# Purpose:      Nansat mapping for NCEP GFS model data subset as downloaded
#               by mapper_ncep_wind_online
# Author:       Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# Made for GRIB files downloaded from http://nomads.ncep.noaa.gov/data/gfs4/
import datetime
import json
import pythesint as pti

from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError


class Mapper(VRT):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' Create NCEP VRT '''

        if not gdalDataset:
            raise WrongMapperError(filename)

        geotransform = gdalDataset.GetGeoTransform()
        if (geotransform != (-0.25, 0.5, 0.0, 90.25, 0.0, -0.5) or
                gdalDataset.RasterCount != 2):  # Not water proof
            raise WrongMapperError(filename)

        metaDict = [{'src': {'SourceFilename': filename,
                             'SourceBand': 1},
                     'dst': {'wkv': 'eastward_wind',
                             'height': '10 m'}},
                    {'src': {'SourceFilename': filename,
                             'SourceBand': 2},
                     'dst': {'wkv': 'northward_wind',
                             'height': '10 m'}},
                    {'src': [{'SourceFilename': filename,
                              'SourceBand': 1,
                              'DataType': gdalDataset.GetRasterBand(1).DataType
                              },
                             {'SourceFilename': filename,
                              'SourceBand': 2,
                              'DataType': gdalDataset.GetRasterBand(2).DataType
                              }],
                     'dst': {'wkv': 'wind_speed',
                             'PixelFunctionType': 'UVToMagnitude',
                             'name': 'windspeed',
                             'height': '2 m'
                             }},
                    {'src': [{'SourceFilename': filename,
                              'SourceBand': 1,
                              'DataType': gdalDataset.GetRasterBand(1).DataType
                              },
                             {'SourceFilename': filename,
                              'SourceBand': 2,
                              'DataType': gdalDataset.GetRasterBand(2).DataType
                              }],
                     'dst': {'wkv': 'wind_from_direction',
                             'PixelFunctionType': 'UVToDirectionFrom',
                             'name': 'winddirection',
                             'height': '2 m'
                             }
                     }]

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # Adding valid time from the GRIB file to dataset
        validTime = gdalDataset.GetRasterBand(1).GetMetadata()['GRIB_VALID_TIME']
        self.dataset.SetMetadataItem('time_coverage_start',
            (datetime.datetime.utcfromtimestamp(
                int(validTime.strip().split(' ')[0])).isoformat()))
        self.dataset.SetMetadataItem('time_coverage_end',
            (datetime.datetime.utcfromtimestamp(
                int(validTime.strip().split(' ')[0])).isoformat()))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('ncep-gfs')

        # TODO: Validate that the found instrument and platform are indeed what we
        # want....

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
