# Name:         mapper_occci_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import datetime as dt
import numpy as np
import os

from nansat.nsr import NSR
from nansat.mappers.opendap import Dataset, Opendap

#http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03/20121001000000-METNO-L4_GHRSST-SSTfnd-METNO_OI-ARC-v02.0-fv01.0.nc

class Mapper(Opendap):
    ''' VRT with mapping of WKV for NCEP GFS '''
    baseURLs = ['http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03/']
    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    t0 = dt.datetime(1981,01,01)
    srcDSProjection = NSR().wkt

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 date=None, ds=None, bands=None, cachedir=None,
                 **kwargs):
        ''' Create NCEP VRT
        Parameters:
            fileName : URL
            date : str
                2010-05-01
            ds : netCDF.Dataset
                previously opened dataset

        '''
        self.create_vrt(fileName, gdalDataset, gdalMetadata, date, ds, bands, cachedir)


    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''

        dsDatetimes = np.array([np.datetime64(self.t0 + dt.timedelta(seconds=int(day)))
                                for day in dsTime]).astype('M8[s]')
        return dsDatetimes
