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


class Opendap(object):
    ''' Methods for all OpenDAP mappers '''

    def get_dataset(self, url, ds):
        ''' Open Dataset '''
        if ds is None:
            try:
                ds = Dataset(url)
            except:
                raise OptionError('Cannot open %s' % url)
        elif type(ds) != Dataset:
            raise OptionError('Input ds is not netCDF.Dataset!')

        return ds

    def get_geospatial_variable_names(self, ds):
        ''' Get names of variables with both spatial dimentions'''
        dsNames = []
        for var in ds.variables:
            if     (self.xName in ds.variables[var].dimensions and
                    self.yName in ds.variables[var].dimensions):
                dsNames.append(var)

        return dsNames

    def get_dataset_datetimes(self, ds):
        ''' Convert time variable to np.datetime64 '''
        dsTime = ds.variables[self.timeVarName][:]
        dsDatetimes = np.array(
                        [np.datetime64(self.timeCalendarStart) + day
                         for day in dsTime]).astype('M8[s]')

        return dsDatetimes

    def get_layer_datetime(self, date, ds):
        ''' Get datetime of the matching layer and layer number '''

        datetimes = self.get_dataset_datetimes(ds)
        datetimeResolution = np.abs(datetimes[0] - datetimes[1])
        date = np.datetime64(date).astype('M8[s]')
        matchingDateDiff = np.min(np.abs(datetimes - date))
        if matchingDateDiff > datetimeResolution:
            raise OptionError('Date %s is out of range' % date)
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
                            }
                    }

        #for attr in self.ds.variables[varName].ncattrs():
        #    metaItem['dst'][str(attr)] = str(self.ds.variables[varName].getncattr(attr))

        return metaItem

class Mapper(Opendap, VRT):
    ''' VRT with mapping of WKV for NCEP GFS '''

    baseURL = 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0'
    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    timeCalendarStart = '1970-01-01'

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 date='2010-05-01', ds=None, bands=None, **kwargs):
        ''' Create NCEP VRT
        Parameters:
            fileName : URL
            date : str
                2010-05-01
            ds : netCDF.Dataset
                previously opened dataset

        '''

        ### TODOs:
        # add <bands> parameter to reduce number of bands

        # get type of variable from Dataset and add to metaItem

        # add time

        # add metadata

        # CACHING! <cache>[=False] specify where to cache time variable

        # <layerNumber> - integer of layer (instead of time var reading)

        # generic way to generate VRT: read x,y dimensions, first, second values

        if not fileName.startswith(self.baseURL):
            raise WrongMapperError

        ds = self.get_dataset(fileName, ds)

        dsLayerNo, dsLayerDate = self.get_layer_datetime(date, ds)

        if bands is None:
            dsVarNames = self.get_geospatial_variable_names(ds)
        else:
            dsVarNames = bands

        # create VRT with correct lon/lat (geotransform)
        VRT.__init__(self, srcProjection=NSR().wkt,
                     srcRasterXSize=8640,
                     srcRasterYSize=4320,
                     srcGeoTransform=(-179.97920227, 0.04167175, 0, 89.97915649, 0, -0.04166412))

        metaDict = [self.get_metaitem(fileName, dsVarName, dsLayerNo)
                      for dsVarName in dsVarNames]

        self._create_bands(metaDict)

        # set time
        #self.dataset.SetMetadataItem('time_coverage_start', gcDate.isoformat())
