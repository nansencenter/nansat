# Name:         mapper_globcurrent_online.py
# Purpose:      Nansat mapping for GLOBCURRENT data, stored online in HYRAX
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
    #http://www.ifremer.fr/opendap/cerdap1/globcurrent/v2.0/global_025_deg/total_hs/2010/001/20100101000000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc
    #http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V02.0_FULL_TIME_SERIE
    baseURLs = ['http://www.ifremer.fr/opendap/cerdap1/globcurrent/v2.0/']
    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    timeCalendarStart = '1950-01-01'

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
        fname = os.path.split(filename)[1]
        date = '%s-%s-%sT%s:00Z' % (fname[0:4], fname[4:6], fname[6:8], fname[8:10])

        self.create_vrt(filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir)

        # add instrument and platform
        #instr = pti.get_gcmd_instrument('active remote sensing')
        #pltfr = pti.get_gcmd_platform('Earth Observation Satellites')
        pltfr = pti.get_gcmd_platform('JASON-1')
        instr = pti.get_gcmd_instrument('JASON-2 RADAR ALTIMETER')

        self.dataset.SetMetadataItem('instrument',  json.dumps(instr))
        self.dataset.SetMetadataItem('platform',    json.dumps(pltfr))
        self.dataset.SetMetadataItem('Data Center', 'FR/IFREMER/CERSAT')
        self.dataset.SetMetadataItem('Entry Title', 'GLOBCURRENT')

    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''
        dsDatetimes = np.array([(np.datetime64(self.timeCalendarStart).astype('M8[s]') +
                                 np.timedelta64(int(day), 'D').astype('m8[s]') +
                                 np.timedelta64(int(24*(day - int(day))), 'h').astype('m8[s]'))
                                for day in dsTime]).astype('M8[s]')

        return dsDatetimes
