# Name:         mapper_ncep_wind_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# Usage:
#    w = Nansat('occci_online:1D:chlor_a:2010-01-01')

import os
import datetime
from dateutil.parser import parse
from time import sleep as time_sleep

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    raise ImportError('''
         Cannot import Dataset from netCDF4.
         You cannot access OC CCI data but
         Nansat will work.''')

from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.tools import gdal, WrongMapperError, OptionError

class Mapper(VRT, object):
    ''' VRT with mapping of WKV for NCEP GFS '''
    SST_CCI_URL_FORMAT = 'http://dap.ceda.ac.uk/data/neodc/esacci_sst/data/lt/Analysis/L4/v01.0/%Y/%m/%d/%Y%m%d120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.0.nc'
    cachePrefix = '_'.join(os.path.split(SST_CCI_URL_FORMAT)[1].split('-')[1:])

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 cache='', lons=None, lats=None, **kwargs):
        ''' Create NCEP VRT
        Parameters:
            fileName : str
                sstcci_online:analysed_sst:2010-01-01
                sstcci_online:analysis_error:2010-01-01
                sstcci_online:sea_ice_fraction:2010-01-01
                sstcci_online:sea_ice_fraction_error:2010-01-01
                sstcci_online:mask:2010-01-01
            cache : str or bool
                if str - name of the cahcing directory
                If False or None - no caching
            lon : list
                minimum and maimum values of longitude
            lat : list
                minimum and maimum values of latitude
        '''


        keywordBase = 'sstcci_online'
        if not fileName.startswith(keywordBase):
            raise WrongMapperError

        # create caching directory
        if cache == '':
            cache = os.path.curdir
        if cache and not os.path.exists(cache):
            os.mkdir(cache)

        # Get prod name
        prodName = fileName.split(':')[1]

        # Get date
        iDate = parse(fileName.split(':')[2])

        # create dataset URL
        dsURL = iDate.strftime(self.SST_CCI_URL_FORMAT)

        # get lon, lat, time dimensions from the OC CCI Dataset
        self.lon, self.lat = self.get_lon_lat(cache, dsURL)

        # get rows and cols that contain predefined spatial domain
        self.rows, self.cols, lons, lats, geoTransform = self.get_rows_cols(lons, lats)

        # create VRT with correct lon/lat (geotransform)
        VRT.__init__(self, srcProjection=NSR().wkt,
                     srcRasterXSize=len(self.cols),
                     srcRasterYSize=len(self.rows),
                     srcGeoTransform=geoTransform)

        # Get SourceFilename either from memory array or from cached file
        sourceFilename = self.get_sourcefilename(cache, dsURL, iDate, prodName, lons, lats)

        metaDict = [{'src': {
                        'SourceFilename': sourceFilename,
                        'SourceBand': 1,
                            },
                    'dst': {
                        'name': prodName,
                            }
                    }]

        self._create_bands(metaDict)

        # set time
        self.dataset.SetMetadataItem('time_coverage_start', iDate.isoformat())

    def get_sourcefilename(self, cache, dsURL, iDate, prodName, lons, lats):
        ''' Get SourceFilename either from memory array or from cached file '''
        print 'Get ', iDate, prodName
        # try to find cached layer
        if cache:
            layerFilename = os.path.join(cache,
                                     '%s_%s_%s_%+04d%+04d%+04d%+04d.tif' % (
                                                        self.cachePrefix,
                                                        prodName, iDate.strftime('%Y%m%d'),
                                                        min(lons), max(lons),
                                                        min(lats), max(lats)))
            print 'from ', layerFilename, '...'
            if os.path.exists(layerFilename):
                print 'from ', layerFilename
                return layerFilename

        print 'from THREDDS'
        ### Continue without pre-cached file
        # get product array from remote dataset
        ds = Dataset(dsURL)
        prodArray = ds.variables[prodName][0, min(self.rows):max(self.rows),
                                              min(self.cols):max(self.cols)]
        prodArray.data[prodArray.mask] = np.nan
        # create VRT and add to self.bandVRTs
        vrt = VRT(array=prodArray.data, srcProjection=NSR().wkt,
                    srcRasterXSize=self.dataset.RasterXSize,
                    srcRasterYSize=self.dataset.RasterYSize,
                    srcGeoTransform=self.dataset.GetGeoTransform())
        sourceFilename = vrt.fileName
        self.bandVRTs[os.path.split(sourceFilename)[1]] = vrt
        if cache:
            gdal.GetDriverByName('GTiff').CreateCopy(layerFilename, vrt.dataset)

        return sourceFilename

    def get_lon_lat(self, cache, dsURL):
        ### Get TIME, LAT, LON
        # first try from cache
        print 'Get lon, lat'
        lon, lat  = None, None
        if cache:
            gridFile = os.path.join(cache, self.cachePrefix+'_grid.npz')
            if os.path.exists(gridFile):
                try:
                    lon = np.load(gridFile)['lon']
                    lat = np.load(gridFile)['lat']
                except:
                    time_sleep(0.5)
                    lon = np.load(gridFile)['lon']
                    lat = np.load(gridFile)['lat']

        # if cache does not exist try to fetch from remote dataset
        if lon is None:
            ds = Dataset(dsURL)
            lat = ds.variables['lat'][:]
            lon = ds.variables['lon'][:]

        # cache grid specs
        if cache:
            np.savez_compressed(gridFile, lon=lon, lat=lat)

        return lon, lat

    def get_rows_cols(self, lons, lats):
        ''' Get rows and cols, estimate actual min/max of lat/lon'''

        ### Get min/max lon/lat
        if lons is None:
            lons = [-180, 180]
        if type(lons) in [int, float]:
            lons = [lons]
        if lats is None:
            lats = [-90, 90]
        if type(lats) in [int, float]:
            lats = [lats]

        rows = np.nonzero((self.lat >= min(lats)) * (self.lat <= max(lats)))[0]
        cols = np.nonzero((self.lon >= min(lons)) * (self.lon <= max(lons)))[0]

        lons = [min(self.lon[cols]), max(self.lon[cols])]
        lats = [min(self.lat[rows]), max(self.lat[rows])]

        geoTransform = (self.lon[cols][0], (self.lon[cols][-1] - self.lon[cols][0]) / len(cols), 0,
                        self.lat[rows][0], 0, (self.lat[rows][-1] - self.lat[rows][0]) / len(rows))

        return rows, cols, lons, lats, geoTransform
