# Name:         opendap.py
# Purpose:      Abstract class Opendap is extended by some mappers
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import os
import datetime
from dateutil.parser import parse
from time import sleep as time_sleep
import warnings
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


class Opendap(VRT):
    ''' Methods for all OpenDAP mappers '''

    def get_dataset(self, ds):
        ''' Open Dataset '''
        if ds is None:
            try:
                ds = Dataset(self.fileName)
            except:
                raise OptionError('Cannot open %s' % self.fileName)
        elif type(ds) != Dataset:
            raise OptionError('Input ds is not netCDF.Dataset!')

        return ds

    def get_geospatial_variable_names(self):
        ''' Get names of variables with both spatial dimentions'''
        dsNames = []
        for var in self.ds.variables:
            if     (self.xName in self.ds.variables[var].dimensions and
                    self.yName in self.ds.variables[var].dimensions):
                dsNames.append(var)

        return dsNames

    def get_dataset_time(self):
        ''' Load data from time variable '''
        cachefile = ''
        if (type(self.cachedir) in [str, unicode] and
            os.path.isdir(self.cachedir)):
            # do caching
            cachefile = '%s_%s.npz' % (os.path.join(self.cachedir,
                                       os.path.split(self.fileName)[1]),
                                       self.timeVarName)

        if os.path.exists(cachefile):
            dsTime = np.load(cachefile)[self.timeVarName]
        else:
            warnings.warn('Time consuming loading time from OpenDAP...')
            dsTime = self.ds.variables[self.timeVarName][:]

        if os.path.exists(cachefile):
            np.savez(cachefile, **{self.timeVarName: dsTime})

        return dsTime

    def get_layer_datetime(self, date, datetimes):
        ''' Get datetime of the matching layer and layer number '''

        datetimeResolution = np.abs(datetimes[0] - datetimes[1])
        date = np.datetime64(date).astype('M8[s]')
        matchingDateDiff = np.min(np.abs(datetimes - date))
        if matchingDateDiff > datetimeResolution:
            raise OptionError('Date %s is out of range' % date)
        layerNumber = np.argmin(np.abs(datetimes - date))
        layerStartDate = datetimes[layerNumber]
        layerEndDate = layerStartDate + datetimeResolution

        return layerNumber, layerStartDate, layerEndDate

    def get_metaitem(self, url, varName, layerNo):
        ''' Set metadata for creating band VRT '''

        metaItem = {'src': {
                        'SourceFilename': '%s?%s.%s[%d][y][x]' % (url, varName, varName, layerNo),
                        'SourceBand': 1,
                            },
                    'dst': {
                        'name': varName,
                            }
                    }

        #for attr in self.ds.variables[varName].ncattrs():
        #    metaItem['dst'][str(attr)] = str(self.ds.variables[varName].getncattr(attr))

        return metaItem

    def create_vrt(self, fileName, gdalDataset, gdalMetadata, date, ds, bands, cachedir):
        ''' Create VRT '''
        if not fileName.startswith(self.baseURL):
            raise WrongMapperError

        self.fileName = fileName
        self.cachedir = cachedir
        self.ds = self.get_dataset(ds)

        dsTime = self.get_dataset_time()

        dsDatetimes = self.convert_dstime_datetimes(dsTime)

        (dsLayerNo,
         dsLayerStartDate,
         dsLayerEndDate) = self.get_layer_datetime(date, dsDatetimes)

        if bands is None:
            dsVarNames = self.get_geospatial_variable_names()
        else:
            dsVarNames = bands

        # create VRT with correct lon/lat (geotransform)
        VRT.__init__(self, srcProjection=self.srcDSProjection,
                     srcRasterXSize=self.srcDSRasterXSize,
                     srcRasterYSize=self.srcDSRasterYSize,
                     srcGeoTransform=self.srcDSGeoTransform)

        metaDict = [self.get_metaitem(fileName, dsVarName, dsLayerNo)
                      for dsVarName in dsVarNames]

        self._create_bands(metaDict)

        # set time
        self.dataset.SetMetadataItem('time_coverage_start', str(dsLayerStartDate))
        self.dataset.SetMetadataItem('time_coverage_end', str(dsLayerEndDate))
