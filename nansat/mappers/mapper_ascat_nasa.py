# Name:        mapper_ascat_nasa
# Purpose:     Mapping for ASCAT scatterometer winds
# Authors:     Knut-Frode Dagestad, Morten W. Hansen
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html
#
# For NetCDF files of ASCAT wind data from the NASA JPL archive:
# ftp://podaac-ftp.jpl.nasa.gov/allData/ascat/preview/L2/metop_a/12km/
import os.path
import datetime
import warnings

import json
import pythesint as pti

from nansat.tools import gdal, ogr
from nansat.vrt import VRT, GeolocationArray
from nansat.tools import WrongMapperError


class Mapper(VRT):
    ''' Create VRT with mapping of WKV '''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 latlonGrid=None, mask='', **kwargs):

        ''' Create VRT

        Parameters
        -----------
        fileName : string
        gdalDataset : gdal dataset
        gdalMetadata : gdal metadata
        latlonGrid : numpy 2 layered 2D array with lat/lons of desired grid
        '''
        # test if input files is ASCAT
        iDir, iFile = os.path.split(fileName)
        iFileName, iFileExt = os.path.splitext(iFile)
        try:
            assert iFileName[0:6] == 'ascat_' and iFileExt == '.nc'
        except:
            raise WrongMapperError

        # Create geolocation
        subDataset = gdal.Open('NETCDF:"' + fileName + '":lat')
        self.GeolocVRT = VRT(srcRasterXSize=subDataset.RasterXSize,
                             srcRasterYSize=subDataset.RasterYSize)

        GeolocMetaDict = [{'src': {'SourceFilename': ('NETCDF:"' + fileName +
                                                      '":lon'),
                                   'SourceBand': 1,
                                   'ScaleRatio': 0.00001,
                                   'ScaleOffset': -360},
                           'dst': {}},
                          {'src': {'SourceFilename': ('NETCDF:"' + fileName +
                                                      '":lat'),
                                   'SourceBand': 1,
                                   'ScaleRatio': 0.00001,
                                   'ScaleOffset': 0},
                           'dst': {}}]

        self.GeolocVRT._create_bands(GeolocMetaDict)

        GeolocObject = GeolocationArray(xVRT=self.GeolocVRT,
                                        yVRT=self.GeolocVRT,
                                        # x = lon, y = lat
                                        xBand=1, yBand=2,
                                        lineOffset=0, pixelOffset=0,
                                        lineStep=1, pixelStep=1)

        # create empty VRT dataset with geolocation only
        VRT.__init__(self,
                     srcRasterXSize=subDataset.RasterXSize,
                     srcRasterYSize=subDataset.RasterYSize,
                     gdalDataset=subDataset,
                     geolocationArray=GeolocObject,
                     srcProjection=GeolocObject.d['SRS'])

        # Scale and NODATA should ideally be taken directly from raw file
        metaDict = [{'src': {'SourceFilename': ('NETCDF:"' + fileName +
                                                '":wind_speed'),
                             'ScaleRatio': 0.01,
                             'NODATA': -32767},
                     'dst': {'name': 'windspeed',
                             'wkv': 'wind_speed'}
                     },
                    {'src': {'SourceFilename': ('NETCDF:"' + fileName +
                                                '":wind_dir'),
                             'ScaleRatio': 0.1,
                             'NODATA': -32767},
                     'dst': {'name': 'winddirection',
                             'wkv': 'wind_from_direction'}}]

        self._create_bands(metaDict)

        # This should not be necessary
        # - should be provided by GeolocationArray!
        self.dataset.SetProjection(GeolocObject.d['SRS'])

        # Add time
        startTime = datetime.datetime(int(iFileName[6:10]),
                                      int(iFileName[10:12]),
                                      int(iFileName[12:14]),
                                      int(iFileName[15:17]),
                                      int(iFileName[17:19]),
                                      int(iFileName[19:21]))
        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start', startTime.isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', startTime.isoformat())

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('ascat')
        ee = pti.get_gcmd_platform('metop-a')

        # TODO: Validate that the found instrument and platform are indeed what
        # we want....

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
