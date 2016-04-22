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
from nansat.mappers.opendap import Opendap

class Mapper(Opendap):
    ''' VRT with mapping of WKV for NCEP GFS '''
            #  http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh-agg
            #  http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh/osisaf-nh_aggregated_ice_concentration_nh_polstere-100_200703010000.nc
    baseURL = 'http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh'
    timeVarName = 'time'
    xName = 'xc'
    yName = 'yc'
    timeCalendarStart = '1978-01-01'

    srcDSProjection = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45 +units=km').wkt
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
        if fileName[-3:] == '.nc':
            fname = os.path.split(fileName)[1]
            fdate = fname.split('.')[0].split('_')[-1]
            date = '%s-%s-%s' % (fdate[0:4], fdate[4:6], fdate[6:8])

        self.create_vrt(fileName, gdalDataset, gdalMetadata, date, ds, bands, cachedir)


    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''
        t0 = dt.datetime.strptime(self.timeCalendarStart, '%Y-%m-%d')
        dsDatetimes = np.array([np.datetime64(t0 + dt.timedelta(seconds=day))
                                for day in dsTime]).astype('M8[s]')
        return dsDatetimes
