# Name:         opendap.py
# Purpose:      Abstract class Opendap is extended by some mappers
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

# http://cfconventions.org/wkt-proj-4.html

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
from nansat.tools import gdal

from nansat.exceptions import WrongMapperError


class Opendap(VRT):
    ''' Methods for all OpenDAP mappers '''

    P2S = {
        'H': 60*60,
        'D': 86400,
        'M': 30*24*60*60,
        'Y': 31536000,
        }

    ### TODOs:
    # add band metadata

    def test_mapper(self, filename):
        ''' Tests if filename fits mapper. May raise WrongMapperError '''
        baseURLmatch = False
        for baseURL in self.baseURLs:
            if filename.startswith(baseURL):
                baseURLmatch = True
                break
        if not baseURLmatch:
            raise WrongMapperError(filename)


    def get_dataset(self, ds):
        ''' Open Dataset '''
        if ds is None:
            try:
                ds = Dataset(self.filename)
            except:
                raise ValueError('Cannot open %s' % self.filename)
        elif type(ds) != Dataset:
            raise ValueError('Input ds is not netCDF.Dataset!')

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
                                       os.path.split(self.filename)[1]),
                                       self.timeVarName)

        if os.path.exists(cachefile):
            dsTime = np.load(cachefile)[self.timeVarName]
        else:
            warnings.warn('Time consuming loading time from OpenDAP...')
            dsTime = self.ds.variables[self.timeVarName][:]
            warnings.warn('Loading time - OK!')

        if os.path.exists(cachefile):
            np.savez(cachefile, **{self.timeVarName: dsTime})

        return dsTime

    def get_layer_datetime(self, date, datetimes):
        ''' Get datetime of the matching layer and layer number '''

        if len(datetimes) == 1 or date is None:
            layerNumber = 0
        else:
            # find closest layer
            datetimeResolution = np.abs(datetimes[0] - datetimes[1])
            date = np.datetime64(date).astype('M8[s]')
            matchingDateDiff = np.min(np.abs(datetimes - date))
            if matchingDateDiff > datetimeResolution:
                raise ValueError('Date %s is out of range' % date)
            layerNumber = np.argmin(np.abs(datetimes - date))

        layerDate = datetimes[layerNumber]

        return layerNumber, layerDate

    def get_metaitem(self, url, varName, layerNo):
        ''' Set metadata for creating band VRT '''

        metaItem = {'src': {
                        'SourceFilename': '%s?%s.%s[%d][y][x]' % (url, varName, varName, layerNo),
                        'SourceBand': 1,
                            },
                    'dst': {
                        'name': varName,
                        'dataType': 6,
                            }
                    }

        for attr in self.ds.variables[varName].ncattrs():
            attrKey = attr.encode('ascii', 'ignore')
            attrVal = self.ds.variables[varName].getncattr(attr)
            if type(attrVal) in [str, unicode]:
                attrVal = attrVal.encode('ascii', 'ignore')
            if attrKey in ['scale', 'scale_factor']:
                metaItem['src']['ScaleRatio'] = attrVal
            elif attrKey in ['offset', 'add_offset']:
                metaItem['src']['ScaleOffset'] = attrVal
            else:
                metaItem['dst'][attrKey] = str(attrVal)

        return metaItem

    def create_vrt(self, filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir):
        ''' Create VRT '''
        if date is None:
            warnings.warn('''
            Date is not specified! Will return the first layer.
            Please add date="YYYY-MM-DD"''')

        self.filename = filename
        self.cachedir = cachedir
        self.ds = self.get_dataset(ds)

        dsTime = self.get_dataset_time()

        dsDatetimes = self.convert_dstime_datetimes(dsTime)

        dsLayerNo, dsLayerDate = self.get_layer_datetime(date, dsDatetimes)

        if bands is None:
            dsVarNames = self.get_geospatial_variable_names()
        else:
            dsVarNames = bands

        # create VRT with correct lon/lat (geotransform)
        srcRasterXSize, srcRasterYSize = self.get_shape()
        srcGeoTransform = self.get_geotransform()
        self._init_from_dataset_params(srcRasterXSize, srcRasterYSize, srcGeoTransform, self.srcDSProjection)

        metaDict = [self.get_metaitem(filename, dsVarName, dsLayerNo)
                      for dsVarName in dsVarNames]

        self.create_bands(metaDict)

        # set time
        timeResSecs = self.get_time_coverage_resolution()
        self.dataset.SetMetadataItem('time_coverage_start', str(dsLayerDate))
        self.dataset.SetMetadataItem('time_coverage_end', str(dsLayerDate + timeResSecs))

    def get_time_coverage_resolution(self):
        ''' Try to fecth time_coverage_resolution and convert to seconds '''
        timeResSecs = 0
        if 'time_coverage_resolution' in self.ds.ncattrs():
            time_res = self.ds.time_coverage_resolution
            try:
                timeResSecs = int(time_res[1]) * self.P2S[time_res[2].upper()]
            except:
                warnings.warn('Cannot get time_coverage_resolution')

        return timeResSecs

    def get_shape(self):
        ''' Get srcRasterXSize and srcRasterYSize from OpenDAP '''
        return (self.ds.variables[self.xName].size,
                self.ds.variables[self.yName].size)


    def get_geotransform(self):
        ''' Get first two values of X,Y variables and create geoTranform '''

        xx = self.ds.variables[self.xName][0:2]
        yy = self.ds.variables[self.yName][0:2]
        return (xx[0], xx[1]-xx[0], 0, yy[0], 0, yy[1]-yy[0])
