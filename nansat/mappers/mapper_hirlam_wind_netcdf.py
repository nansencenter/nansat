# Name:         mapper_hirlam2nc.py
# Purpose:      Mapper for Hirlam wind data converted from felt to netCDF
# Authors:      Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import datetime

from nansat.tools import gdal, ogr, WrongMapperError
from nansat.vrt import VRT


class Mapper(VRT):
    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30,
                 **kwargs):

        if not gdalMetadata:
            raise WrongMapperError

        isHirlam = False
        for key in gdalMetadata.keys():
            if 'creation by fimex from file' in gdalMetadata[key]:
                isHirlam = True

        if not isHirlam:
            raise WrongMapperError

        #GeolocMetaDict = [{'src':
        #        {'SourceFilename': 'NETCDF:"' + fileName + '":longitude',
        #         'SourceBand': 1,
        #         'ScaleRatio': 1,
        #         'ScaleOffset': 0},
        #     'dst': {}},
        #           {'src':
        #        {'SourceFilename': 'NETCDF:"' + fileName + '":latitude',
        #         'SourceBand': 1,
        #         'ScaleRatio': 1,
        #         'ScaleOffset': 0},
        #     'dst': {}}]

        subDataset = gdal.Open('NETCDF:"' + fileName + '":x_wind_10m')
        #self.GeolocVRT = VRT(srcRasterXSize=subDataset.RasterXSize,
        #                srcRasterYSize=subDataset.RasterYSize)
        #self.GeolocVRT._create_bands(GeolocMetaDict)

        #GeolocObject = GeolocationArray(xVRT=self.GeolocVRT,
        #                                yVRT=self.GeolocVRT,
        #    xBand=1, yBand=2,
        #    lineOffset=0, pixelOffset=0,
        #    lineStep=1, pixelStep=1)

        ## create empty VRT dataset with geolocation only
        #VRT.__init__(self, srcRasterXSize = subDataset.RasterXSize,
        #                   srcRasterYSize = subDataset.RasterYSize,
        #                geolocationArray = GeolocObject,
        #                srcProjection = GeolocObject.d['SRS'])
        lon = gdal.Open(
            'NETCDF:"' + fileName + '":longitude"').ReadAsArray()
        lat = gdal.Open(
            'NETCDF:"' + fileName + '":latitude"').ReadAsArray()
        VRT.__init__(self, lat=lat, lon=lon)

        # Add bands with wind components
        metaDict = [{'src': {'SourceFilename': ('NETCDF:"' + fileName +
                                                '":x_wind_10m'),
                             'NODATA': -32767},
                     'dst': {'name': 'U',
                             'wkv': 'eastward_wind'}},
                    {'src': {'SourceFilename': ('NETCDF:"' + fileName +
                                                '":y_wind_10m'),
                             'NODATA': -32767},
                     'dst': {'name': 'V',
                             'wkv': 'northward_wind'}}]

        # Add pixel function with wind speed
        metaDict.append({'src': [{'SourceFilename': ('NETCDF:"' + fileName +
                                                     '":x_wind_10m'),
                                  'SourceBand': 1,
                                  'DataType': 6},
                                 {'SourceFilename': ('NETCDF:"' + fileName +
                                                     '":y_wind_10m'),
                                  'SourceBand': 1,
                                  'DataType': 6}],
                         'dst': {'wkv': 'wind_speed',
                                 'name': 'windspeed',
                                 'height': '10 m',
                                 'PixelFunctionType': 'UVToMagnitude',
                                 'NODATA': 9999}})

        # add bands with metadata and corresponding values
        # to the empty VRT
        self._create_bands(metaDict)

        # Add valid time
        validTime = datetime.datetime.utcfromtimestamp(
            int(subDataset.GetRasterBand(1).
                GetMetadata()['NETCDF_DIM_time']))
        self.dataset.SetMetadataItem('time_coverage_start', validTime.isoformat())
