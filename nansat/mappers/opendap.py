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
        if (type(self.cachedir) in [str] and
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

    def get_metaitem(self, url, var_name, var_dimensions):
        """Set metadata for creating band VRT"""
        # assemble dimensions string
        dims = ''.join(['[%s]' % dim for dim in var_dimensions])
        meta_item = {
            'src': {'SourceFilename': '{url}?{var}.{var}{shape}'.format(url=url,
                                                                        var=var_name,
                                                                        shape=dims),
                    'SourceBand': 1},
            'dst': {'name': var_name,
                    'dataType': 6}
        }

        for attr in self.ds.variables[var_name].ncattrs():
            attrKey = attr.encode('ascii', 'ignore')
            attrVal = self.ds.variables[var_name].getncattr(attr)
            if type(attrVal) in [str]:
                attrVal = attrVal.encode('ascii', 'ignore')
            if attrKey in ['scale', 'scale_factor']:
                meta_item['src']['ScaleRatio'] = attrVal
            elif attrKey in ['offset', 'add_offset']:
                meta_item['src']['ScaleOffset'] = attrVal
            else:
                meta_item['dst'][attrKey] = str(attrVal)

        return meta_item

    def create_vrt(self, filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir):
        """Create VRT"""
        if date is None:
            warnings.warn('Date is not specified! Will return the first layer. '
                          'Please add date="YYYY-MM-DD"')

        self.filename = filename
        self.cachedir = cachedir
        self.ds = self.get_dataset(ds)

        ds_time = self.get_dataset_time()
        ds_times = self.convert_dstime_datetimes(ds_time)
        layer_time_id, layer_date = self.get_layer_datetime(date, ds_times)

        if bands is None:
            var_names = self.get_geospatial_variable_names()
        else:
            var_names = bands

        # create VRT with correct lon/lat (geotransform)
        raster_x, raster_y = self.get_shape()
        geotransform = self.get_geotransform()
        self._init_from_dataset_params(int(raster_x), int(raster_y),
                                       geotransform, self.srcDSProjection)

        meta_dict = []
        for var_name in var_names:
            # Get list of dimensions for a variable
            var_dimensions = list(self.ds.variables[var_name].dimensions)
            # get variable specific dimensions
            spec_dimension = list(filter(lambda dim_name: dim_name not in ['time', 'y', 'x'],
                                         var_dimensions))[0]

            var_dimensions[var_dimensions.index('time')] = layer_time_id
            if spec_dimension:
                for i in range(self.ds.dimensions[spec_dimension].size):
                    var_dimensions_copy = var_dimensions.copy()
                    var_dimensions_copy[var_dimensions_copy.index(spec_dimension)] = i
                    meta_dict.append(self.get_metaitem(filename, var_name, var_dimensions_copy))
            else:
                meta_dict.append(self.get_metaitem(filename, var_name, var_dimensions))

        self.create_bands(meta_dict)

        # set time
        timeResSecs = self.get_time_coverage_resolution()
        self.dataset.SetMetadataItem('time_coverage_start', str(layer_date))
        self.dataset.SetMetadataItem('time_coverage_end', str(layer_date + timeResSecs))

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
