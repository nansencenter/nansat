# Name:         mapper_metno_hires_seaice.py
# Purpose:      Nansat mapping for high resolution sea ice
#                   from met.no Thredds server
# Authors:      Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

# High resolution (1 km) manual ice concentration, based on SAR imagery
# http://thredds.met.no/thredds/catalog/myocean/siw-tac/siw-metno-svalbard/
#
# Mapper may be called with full URL:
#   'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/siw-metno-svalbard/2014/01/ice_conc_svalbard_201401091500.nc
# or with keyword (fake filename):
#   'metno_hires_seaice:20140109'
#
# The latter syntax will construct the URL,
# and will return the closest available data within +/- 3 days
import sys
import urllib2
from datetime import datetime, timedelta

from nansat.tools import gdal, ogr, osr, WrongMapperError
from nansat.vrt import VRT


class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Met.no seaice '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create VRT '''

        ThreddsBase = 'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/siw-metno-svalbard/'
        # First check if mapper is called with keyword syntax:
        # filename = metno_hires_seaice:YYYYmmdd
        keywordBase = 'metno_hires_seaice'
        foundDataset = False
        if fileName[0:len(keywordBase)] == keywordBase:
            keywordTime = fileName[len(keywordBase)+1:]
            requestedTime = datetime.strptime(keywordTime, '%Y%m%d')
            # Search for nearest available file, within the closest 3 days
            for deltaDay in [0, -1, 1, -2, 2, -3, 3]:
                validTime = (requestedTime + timedelta(days=deltaDay) +
                             timedelta(hours=15))
                fileName = (ThreddsBase +
                            validTime.strftime(
                                '%Y/%m/ice_conc_svalbard_%Y%m%d1500.nc'))
                try:
                    urllib2.urlopen(fileName + '.dds')
                    foundDataset = True
                    # Data is found for this day
                    break
                except:
                    # No data for this day
                    pass

        if not foundDataset:
            raise WrongMapperError

        # Then check if a valid OPeNDAP URL is given
        # (or has been constructed from keyword)
        if fileName[0:len(ThreddsBase)] != ThreddsBase:
            AttributeError("Not Met.no Svalbard-ice Thredds URL")
        else:
            timestr = fileName[-15:-3]
            validTime = datetime.strptime(timestr, '%Y%m%d%H%M')

        fileName = fileName + '?ice_concentration[0][y][x]'
        srcProjection = osr.SpatialReference()
        srcProjection.ImportFromProj4('+proj=stere lon_0=0.0 +lat_0=90 +datum=WGS84 +ellps=WGS84 +units=km +no_defs')
        srcProjection = srcProjection.ExportToWkt()

        # From thredds web, with manual shift
        srcGeotransform = (-1243.008 - 1, 1, 0, -3190.026 - 7, 0, 1)

        # create empty VRT dataset with geolocation only
        VRT.__init__(self,
                     srcGeoTransform=srcGeotransform,
                     srcProjection=srcProjection,
                     srcRasterXSize=3812,
                     srcRasterYSize=2980)

        metaDict = [{'src': {'SourceFilename': fileName,
                             'sourceBand': 1},
                     'dst': {'name': 'sea_ice_area_fraction',
                             'wkv': 'sea_ice_area_fraction'}}]

        # Add band
        self._create_bands(metaDict)

        # Set time
        self.logger.info('Valid time: %s', str(validTime))
        self.dataset.SetMetadataItem('time_coverage_start',
                                     validTime.isoformat())
