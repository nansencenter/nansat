# Name:         mapper_occci_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os

import json
import numpy as np
import datetime as dt

import pythesint as pti

from nansat.nsr import NSR
from nansat.mappers.opendap import Dataset, Opendap

# http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/conc/2016/04/ice_conc_sh_polstere-100_multi_201604261200.nc ice_conc
# http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr/merged/2016/04/ice_drift_nh_polstere-625_multi-oi_201604151200-201604171200.nc dX dY
# http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/2016/04/ice_type_nh_polstere-100_multi_201604151200.nc ice_type
# http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/edge/2016/04/ice_edge_nh_polstere-100_multi_201604241200.nc ice_edge
# http://thredds.met.no/thredds/dodsC/osisaf_test/met.no/ice/drift_lr/merged/2013/09/ice_drift_nh_polstere-625_multi-oi_201309171200-201309191200.nc dX dY


class Mapper(Opendap):
    ''' VRT with mapping of WKV for NCEP GFS '''
    baseURLs = ['http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh',
                'http://thredds.met.no/thredds/dodsC/osisaf_test/met.no/ice/',
                'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/']
    timeVarName = 'time'
    xName = 'xc'
    yName = 'yc'
    t0 = dt.datetime(1978, 1, 1)
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
        ds = Dataset(filename)
        proj4str = '%s +units=%s' % (ds.variables['Polar_Stereographic_Grid'].proj4_string,
                                     ds.variables['xc'].units)
        self.srcDSProjection = NSR(proj4str).wkt
        if filename[-3:] == '.nc':
            date = self.t0 + dt.timedelta(seconds=ds.variables['time'][0])
            date = date.strftime('%Y-%m-%d')

        self.create_vrt(filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir)

        # add instrument and platform
        mm = pti.get_gcmd_instrument('Passive Remote Sensing')
        ee = pti.get_gcmd_platform('Earth Observation Satellites')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''

        dsDatetimes = np.array([np.datetime64(self.t0 + dt.timedelta(seconds=day))
                                for day in dsTime]).astype('M8[s]')
        return dsDatetimes
