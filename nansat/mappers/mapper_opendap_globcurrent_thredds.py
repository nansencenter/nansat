# Name:         mapper_occci_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import json

import numpy as np

import pythesint as pti

from nansat.nsr import NSR
from nansat.mappers.opendap import Opendap

#https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-8DAY
#https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0-MONTHLY
class Mapper(Opendap):
    ''' VRT with mapping of WKV for NCEP GFS '''

    baseURLs = ['http://tds0.ifremer.fr/thredds/dodsC/CLS-L4']
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
        self.create_vrt(filename, gdalDataset, gdalMetadata, date, ds, bands, cachedir)

        # add instrument and platform
        mm = pti.get_gcmd_instrument('Passive Remote Sensing')
        ee = pti.get_gcmd_platform('Earth Observation Satellites')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
        self.dataset.SetMetadataItem('Data Center', 'FR/IFREMER/CERSAT')
        self.dataset.SetMetadataItem('Entry Title', 'GLOBCURRENT')

    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''
        dsDatetimes = np.array([np.datetime64(self.timeCalendarStart) + int(day)
                                for day in dsTime]).astype('M8[s]')

        return dsDatetimes
