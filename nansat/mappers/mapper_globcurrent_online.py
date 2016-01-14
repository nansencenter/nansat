# Name:         mapper_ncep_wind_online.py
# Purpose:      Nansat mapping for GLOBCURRENT data, stored online in THREDDS
# Author:       Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# Usage:
#    w = Nansat('globcurrent:2010-01-01T12')

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

class Mapper(VRT, object):
    ''' VRT with mapping of WKV for NCEP GFS '''

    GLOBCURR_URLS = {
        'EKMAN_15M'    : 'http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREKM_15M-ERAWS_EEM-V01.0_FULL_TIME_SERIE',
        'EKMAN_HS'     : 'http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREKM_HS-ERAWS_EEM-V01.0_FULL_TIME_SERIE',
        'GEOSTROPHIC'  : 'http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CURGEO_0M-ALT_OI-V01.0_FULL_TIME_SERIE',
        'STOKES_DRIFT' : 'http://tds0.ifremer.fr/thredds/dodsC/GC_MOD_STK_GLO_010_WW3_FULL_TIME_SERIE',
        'TIDAL'        : 'http://tds0.ifremer.fr/thredds/dodsC/GC_MOD_TIDE_GLO_010_FES2012_FULL_TIME_SERIE',
        'TOTAL_15M'    : 'http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_15M-ALT_SUM-V01.0_FULL_TIME_SERIE',
        'TOTAL_HS'     : 'http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V01.0_FULL_TIME_SERIE',
        }

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 product='TOTAL_HS', url='', ds=None,
                 keywords=['northward', 'eastward'], **kwargs):
        ''' Create NCEP VRT
        Parameters:
            fileName : str
                globcurrent:2010-01-01T12
            product : str
                one of the GLOBCURRENT products:
                EKMAN_15M
                EKMAN_HS
                GEOSTROPHIC
                STOKES_DRIFT
                TIDAL
                TOTAL_15M
                TOTAL_HS
            url: str
                absolute url of GLOBCURRENT DATA
        '''

        keywordBase = 'globcurrent'
        if not fileName.startswith(keywordBase):
            raise WrongMapperError

        # get dataset URL
        if url == '' and product.upper() in self.GLOBCURR_URLS:
            url = self.GLOBCURR_URLS[product.upper()]

        if ds is None:
            try:
                self.ds = Dataset(url)
            except:
                raise OptionError('Cannot open %s' % url)
        else:
            self.ds = ds

        ### Get Date from input fileName
        iDate = parse(fileName.split(':')[1])

        # get time variable from GC
        gcTime = self.ds.variables['time'][:]

        # compute date in GC calendar (days since 1950-1-1)
        iDateGC = (iDate - datetime.datetime(1950, 1, 1)).total_seconds() / 60. / 60. / 24.

        # return index of closes time
        gcLayer = np.argmin(np.abs(gcTime - iDateGC))

        # convert closest time to datetime
        gcDate = (datetime.datetime(1950, 1, 1) +
                  datetime.timedelta(float(gcTime[gcLayer])))

        # create VRT with correct lon/lat (geotransform)
        VRT.__init__(self, srcProjection=NSR().wkt,
                     srcRasterXSize=3600,
                     srcRasterYSize=1600,
                     srcGeoTransform=(-179.95, 0.1, 0, -79.95, 0, 0.1))

        varNames = [self.get_var_name(keyword) for keyword in keywords]
        metaDict = [self.get_metaitem(varName, gcLayer, url)
                      for varName in varNames
                      if varName is not None]

        self._create_bands(metaDict)

        # set time
        self._set_time(gcDate)

    def get_var_name(self, keyword):
        ''' Get names of variable based on keyword '''
        for varName in self.ds.variables.keys():
            if keyword in varName:
                return str(varName)

    def get_metaitem(self, varName, gcLayer, url):
        ''' Set metadata for creating band VRT '''

        metaItem = {'src': {
                        'SourceFilename': '%s?%s.%s[%d][y][x]' % (url, varName, varName, gcLayer),
                        'SourceBand': 1,
                            },
                    'dst': {
                        'name': varName,
                            }
                    }

        for attr in self.ds.variables[varName].ncattrs():
            metaItem['dst'][str(attr)] = str(self.ds.variables[varName].getncattr(attr))

        return metaItem
