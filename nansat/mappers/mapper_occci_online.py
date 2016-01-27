# Name:         mapper_occci_online.py
# Purpose:      Nansat mapping for OC CCI data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import numpy as np

from nansat.nsr import NSR
from nansat.mappers.opendap import Opendap

class Mapper(Opendap):
    ''' VRT with mapping of WKV for NCEP GFS '''

    baseURL = 'https://rsg.pml.ac.uk/thredds/dodsC/CCI_ALL-v2.0'
    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    timeCalendarStart = '1970-01-01'

    srcDSProjection = NSR().wkt
    srcDSRasterXSize = 8640
    srcDSRasterYSize = 4320
    srcDSGeoTransform = (-179.97920227, 0.04167175, 0, 89.97915649, 0, -0.04166412)

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 date='2010-05-01', ds=None, bands=None, cachedir=None,
                 **kwargs):
        ''' Create NCEP VRT
        Parameters:
            fileName : URL
            date : str
                2010-05-01
            ds : netCDF.Dataset
                previously opened dataset

        '''
        ### TODOs:
        # add metadata

        # generic way to generate VRT: read x,y dimensions, first, second values

        self.create_vrt(fileName, gdalDataset, gdalMetadata, date, ds, bands, cachedir)


    def convert_dstime_datetimes(self, dsTime):
        ''' Convert time variable to np.datetime64 '''
        dsDatetimes = np.array([np.datetime64(self.timeCalendarStart) + day
                                for day in dsTime]).astype('M8[s]')

        return dsDatetimes
