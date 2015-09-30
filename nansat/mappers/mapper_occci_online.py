# Name:         mapper_ncep_wind_online.py
# Purpose:      Nansat mapping for NCEP GFS model data, stored online
# Author:       Knut-Frode Dagestad, Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# Mapper searches two online archives for NCEP GFS grib files
# covering the requested time, and downloads using curl, if found:
#    1. ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/ (~last month)
#    2. http://nomads.ncdc.noaa.gov/data/gfs4/ (back to June 2012, with holes)
#
# Usage:
#    w = Nansat('ncep_wind_online:YYYYMMDDHHMM')
# Example:
#    w = Nansat('ncep_wind_online:201405011000')

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

    OC_CCI_URLS = {
        '1d': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-DAILY',
        '1D': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-DAILY',
        '5d': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-5DAY',
        '5D': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-5DAY',
        '8d': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-8DAY',
        '8D': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-8DAY',
        '1m': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-MONTHLY',
        '1M': 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-MONTHLY',
    }

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 cache='', lons=None, lats=None, **kwargs):
        ''' Create NCEP VRT
        Parameters:
            fileName : str
                occci_online:1D:chlor_a:2010-01-01
                occci_online:8D:kd_490:12/22/2014
            cache : str or bool
                if str - name of the cahcing directory
                If False or None - no caching
            lon : list
                minimum and maimum values of longitude
            lat : list
                minimum and maimum values of latitude
        '''


        keywordBase = 'occci_online'
        if not fileName.startswith(keywordBase):
            raise WrongMapperError

        # create caching directory
        if cache == '':
            cache = os.path.curdir
        if cache and not os.path.exists(cache):
            os.mkdir(cache)

        ### Get temporal resolution
        timeStep = fileName.split(':')[1]
        # get dataset URL. If resolution doesn't match, get monthly dataset URL
        dsURL = self.OC_CCI_URLS.get(timeStep, self.OC_CCI_URLS['1M'])

        ### Get OC CCI product name
        prodName = fileName.split(':')[2]

        ### Get Date from input fileName
        iDate = parse(fileName.split(':')[3])

        # get lon, lat, time dimensions from the OC CCI Dataset
        self.lon, self.lat, self.time = self.get_lon_lat_time(cache, dsURL)

        # Get actual cci date and number of layer in the dataset
        self.date, self.layer = self.get_date_layer(iDate)

        # get rows and cols that contain predefined spatial domain
        self.rows, self.cols, lons, lats, geoTransform = self.get_rows_cols(lons, lats)

        # create VRT with correct lon/lat (geotransform)
        VRT.__init__(self, srcProjection=NSR().wkt,
                     srcRasterXSize=len(self.cols),
                     srcRasterYSize=len(self.rows),
                     srcGeoTransform=geoTransform)

        # Get SourceFilename either from memory array or from cached file
        sourceFilename = self.get_sourcefilename(cache, dsURL, timeStep, prodName, lons, lats)

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
        self._set_time(self.date)

    def get_sourcefilename(self, cache, dsURL, timeStep, prodName, lons, lats):
        ''' Get SourceFilename either from memory array or from cached file '''
        print 'Get ', timeStep, prodName
        # try to find cached layer
        if cache:
            layerFilename = os.path.join(cache,
                                     '%s_%s_%s_%s_%+04d%+04d%+04d%+04d.tif' % (
                                                        os.path.split(dsURL)[1],
                                                        timeStep, prodName, self.date.strftime('%Y%m%d'),
                                                        min(lons), max(lons),
                                                        min(lats), max(lats)))
            if os.path.exists(layerFilename):
                print 'from ', layerFilename
                return layerFilename

        print 'from THREDDS'
        ### Continue without pre-cached file
        # get product array from remote dataset
        ds = Dataset(dsURL)
        prodArray = ds.variables[prodName][self.layer,
                                           min(self.rows):max(self.rows),
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

    def timecci2time(self, timeCCI, dsURL):
        '''' Convert time from CCI units to internal calendar '''
        if 'CCI_ALL-v1.0-8DAY' in dsURL:
            time = np.zeros(timeCCI.shape[0])
            for i, t1 in enumerate(timeCCI):
                dt = parse(''.join(t1).replace('Z',''))
                time[i] = (dt - datetime.datetime(1970, 1, 1)).days
        else:
            time = timeCCI

        return time

    def get_lon_lat_time(self, cache, dsURL):
        ### Get TIME, LAT, LON
        # first try from cache
        print 'Get lon, lat, time'
        lon, lat, time  = None, None, None
        if cache:
            gridFile = os.path.join(cache, os.path.split(dsURL)[1]+'_grid.npz')
            if os.path.exists(gridFile):
                try:
                    lon = np.load(gridFile)['lon']
                    lat = np.load(gridFile)['lat']
                    time = np.load(gridFile)['time']
                except:
                    time_sleep(0.5)
                    lon = np.load(gridFile)['lon']
                    lat = np.load(gridFile)['lat']
                    time = np.load(gridFile)['time']

        # if cache does not exist try to fetch from remote dataset
        if lon is None:
            ds = Dataset(dsURL)
            lat = ds.variables['lat'][:]
            lon = ds.variables['lon'][:]
            timeCCI = ds.variables['time'][:]
            time = self.timecci2time(timeCCI, dsURL)

        # cache grid specs
        if cache:
            np.savez_compressed(gridFile, lon=lon, lat=lat, time=time)

        return lon, lat, time

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

    def get_date_layer(self, iDate):
        ''' Get actual cci date and number of layer in the dataset '''
        iDateCCI = (iDate - datetime.datetime(1970,1,1)).days

        # find number of the layer to fecth data from
        timeDiff = np.abs(self.time - iDateCCI)
        if timeDiff.min() > 31:
            raise OptionError('Date is outside OC CCI range')
        cciLayer = np.argmin(timeDiff)
        cciDate = datetime.datetime(1970,1,1) + datetime.timedelta(int(self.time[cciLayer]))

        return cciDate, cciLayer
