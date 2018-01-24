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
import datetime
import json
import numpy as np
from dateutil.parser import parse
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
        if (geotransform == (-0.25, 0.5, 0.0, 90.25, 0.0, -0.5) or
                geotransform == (-0.5, 1.0, 0.0, 90.5, 0.0, -1.0)):
            if gdalDataset.RasterCount == 4:
                srcBandId = {'temperature': 2,
                             'u-component': 3,
                             'v-component': 4}
            elif gdalDataset.RasterCount == 9:
                srcBandId = {'temperature': 6,
                             'u-component': 8,
                             'v-component': 9}
            else:
                raise WrongMapperError(filename)
        else:
            raise WrongMapperError(filename)  # Not water proof

        # Adding valid time from the GRIB file to dataset
        band = gdalDataset.GetRasterBand(srcBandId['u-component'])
        validTime = band.GetMetadata()['GRIB_VALID_TIME']
        time_isoformat = (datetime.datetime.utcfromtimestamp(
            int(validTime.strip().split(' ')[0])).isoformat())

        # Set band metadata time_iso_8601 for use in OpenWind
        time_iso_8601 = np.datetime64(parse(time_isoformat))
        metaDict = [{'src': {'SourceFilename': filename,
                             'SourceBand': srcBandId['u-component']},
                     'dst': {'wkv': 'eastward_wind',
                             'height': '10 m',
                             'time_iso_8601': time_iso_8601}},
                    {'src': {'SourceFilename': filename,
                             'SourceBand': srcBandId['v-component']},
                     'dst': {'wkv': 'northward_wind',
                             'height': '10 m',
                             'time_iso_8601': time_iso_8601}},
                    {'src': [{'SourceFilename': filename,
                              'SourceBand': srcBandId['u-component'],
                              'DataType': (gdalDataset.GetRasterBand(srcBandId['u-component']).DataType)
                              },
                             {'SourceFilename': filename,
                              'SourceBand': srcBandId['v-component'],
                              'DataType': gdalDataset.GetRasterBand(srcBandId['v-component']).DataType
                              }],
                     'dst': {'wkv': 'wind_speed',
                             'PixelFunctionType': 'UVToMagnitude',
                             'name': 'windspeed',
                             'height': '2 m',
                             'time_iso_8601': time_iso_8601
                             }},
                    {'src': [{'SourceFilename': filename,
                              'SourceBand': srcBandId['u-component'],
                              'DataType': gdalDataset.GetRasterBand(srcBandId['u-component']).DataType
                              },
                             {'SourceFilename': filename,
                              'SourceBand': srcBandId['v-component'],
                              'DataType': gdalDataset.GetRasterBand(srcBandId['v-component']).DataType
                              }],
                     'dst': {'wkv': 'wind_from_direction',
                             'PixelFunctionType': 'UVToDirectionFrom',
                             'name': 'winddirection',
                             'height': '2 m',
                             'time_iso_8601': time_iso_8601
                             }},
                    {'src': {'SourceFilename': filename,
                             'SourceBand': srcBandId['temperature']},
                     'dst': {'wkv': 'air_temperature',
                             'name': 'air_t',
                             'height': '2 m',
                             'time_iso_8601': time_iso_8601}
                     }]

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        self.dataset.SetMetadataItem('time_coverage_start', time_isoformat)

        self.dataset.SetMetadataItem('time_coverage_end', time_isoformat)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('ncep-gfs')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

