# Name:         mapper_occci_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
import json
import datetime as dt

import numpy as np

import pythesint as pti

from nansat.nsr import NSR
from nansat.mappers.opendap import Dataset, Opendap

#http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03/20121001000000-METNO-L4_GHRSST-SSTfnd-METNO_OI-ARC-v02.0-fv01.0.nc
#http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03_V1/20120808-METNO-L4UHfnd-ARC-v01-fv01-METNO_OI.nc
#http://thredds.met.no/thredds/dodsC/sea_ice/SST-METNO-ARC-SST_L4-OBS-V2-V1/sst_arctic_aggregated


class Mapper(Opendap):
    """VRT with mapping of WKV for NCEP GFS"""
    baseURLs = ['http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03/',
                'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/sst-metno-arc-sst03_V1/',
                'http://thredds.met.no/thredds/dodsC/sea_ice/SST-METNO-ARC-SST_L4-OBS-V2-V1/']
    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    t0 = dt.datetime(1981, 1, 1)
    srcDSProjection = NSR().wkt

    def __init__(self, filename, gdalDataset, gdalMetadata, date=None, ds=None, bands=None,
                 cachedir=None, **kwargs):
        ''' Create NCEP VRT
        Parameters:
            filename : URL
            date : str
                2010-05-01
            ds : netCDF.Dataset
                previously opened dataset

        '''
        self.test_mapper(filename)
        self.create_vrt(filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir)

        # add instrument and platform
        mm = pti.get_gcmd_instrument('Passive Remote Sensing')
        ee = pti.get_gcmd_platform('Earth Observation Satellites')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''

        dsDatetimes = np.array([np.datetime64(self.t0 + dt.timedelta(seconds=int(day)))
                                for day in dsTime]).astype('M8[s]')
        return dsDatetimes
