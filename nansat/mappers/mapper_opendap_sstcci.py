# Name:         mapper_ncep_wind_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os

from dateutil.parser import parse
import json

import numpy as np

import pythesint as pti

from nansat.nsr import NSR
from nansat.mappers.opendap import Opendap

class Mapper(Opendap):
    ''' VRT with mapping of WKV for NCEP GFS '''
    #http://dap.ceda.ac.uk/data/neodc/esacci/sst/data/lt/Analysis/L4/v01.1/2010/05/01/20100501120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.1.nc
    baseURLs = ['http://dap.ceda.ac.uk/data/neodc/esacci/sst/data/lt/Analysis/L4/v01.1/']
    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    timeCalendarStart = '1981-01-01'

    srcDSProjection = NSR().wkt
    def __init__(self, filename, gdalDataset, gdalMetadata,
                 date=None, ds=None, bands=None, cachedir=None,
                 **kwargs):
        ''' Create NCEP VRT
        Parameters:
            filename : URL
            date : str
                2010-05-01
            ds : netCDF.Dataset
                previously opened dataset

        '''
        self.test_mapper(filename)
        fname = os.path.split(filename)
        date = '%s-%s-%s' % (fname[0:4], fname[4:6], fname[6:8])

        self.create_vrt(filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir)

        # add instrument and platform
        mm = pti.get_gcmd_instrument('Passive Remote Sensing')
        ee = pti.get_gcmd_platform('Earth Observation Satellites')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''
        dsDatetimes = np.array([(np.datetime64(self.timeCalendarStart).astype('M8[s]') +
                                 np.timedelta64(int(day), 'D').astype('m8[s]') +
                                 np.timedelta64(int(24*(day - int(day))), 'h').astype('m8[s]'))
                                for day in dsTime]).astype('M8[s]')

        return dsDatetimes
