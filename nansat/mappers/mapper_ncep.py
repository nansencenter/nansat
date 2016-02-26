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
import pythesint as pti

from nansat.vrt import VRT
from nansat.tools import WrongMapperError


class Mapper(VRT):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create NCEP VRT '''

        if not gdalDataset:
            raise WrongMapperError

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
                raise WrongMapperError
        else:
            raise WrongMapperError  # Not water proof

        metaDict = [{'src': {'SourceFilename': fileName,
                             'SourceBand': srcBandId['u-component']},
                     'dst': {'wkv': 'eastward_wind',
                             'height': '10 m'}},
                    {'src': {'SourceFilename': fileName,
                             'SourceBand': srcBandId['v-component']},
                     'dst': {'wkv': 'northward_wind',
                             'height': '10 m'}},
                    {'src': [{'SourceFilename': fileName,
                              'SourceBand': srcBandId['u-component'],
                              'DataType': (gdalDataset.GetRasterBand(srcBandId['u-component']).DataType)
                              },
                             {'SourceFilename': fileName,
                              'SourceBand': srcBandId['v-component'],
                              'DataType': gdalDataset.GetRasterBand(srcBandId['v-component']).DataType
                              }],
                     'dst': {'wkv': 'wind_speed',
                             'PixelFunctionType': 'UVToMagnitude',
                             'name': 'windspeed',
                             'height': '2 m'
                             }},
                    {'src': [{'SourceFilename': fileName,
                              'SourceBand': srcBandId['u-component'],
                              'DataType': gdalDataset.GetRasterBand(srcBandId['u-component']).DataType
                              },
                             {'SourceFilename': fileName,
                              'SourceBand': srcBandId['v-component'],
                              'DataType': gdalDataset.GetRasterBand(srcBandId['v-component']).DataType
                              }],
                     'dst': {'wkv': 'wind_from_direction',
                             'PixelFunctionType': 'UVToDirectionFrom',
                             'name': 'winddirection',
                             'height': '2 m'
                             }},
                    {'src': {'SourceFilename': fileName,
                             'SourceBand': srcBandId['temperature']},
                     'dst': {'wkv': 'air_temperature',
                             'name': 'air_t',
                             'height': '2 m'}
                     }]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # Adding valid time from the GRIB file to dataset
        band = gdalDataset.GetRasterBand(srcBandId['u-component'])
        validTime = band.GetMetadata()['GRIB_VALID_TIME']

        self.dataset.SetMetadataItem('time_coverage_start',
            (datetime.datetime.utcfromtimestamp(
                int(validTime.strip().split(' ')[0])).isoformat()))

        self.dataset.SetMetadataItem('time_coverage_end',
            ((datetime.datetime.utcfromtimestamp(
                int(validTime.strip().split(' ')[0]))
             + datetime.timedelta(hours=3)).isoformat()))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('ncep-gfs')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
