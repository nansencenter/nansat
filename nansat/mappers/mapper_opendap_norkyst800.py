# Name:         mapper_opendap_norkyst800.py
# Purpose:      Nansat mapping for ROMS Norkyst-800 ocean model data
#               provided by MET Norway
# Author:       Artem Moiseev
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import os
import json
import pythesint as pti
import numpy as np
from datetime import datetime
from nansat.mappers.opendap import Opendap
from nansat.nsr import NSR
from netCDF4 import Dataset


class Mapper(Opendap):

    baseURLs = ['http://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/',
                'https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/']

    timeVarName = 'time'
    xName = 'X'
    yName = 'Y'
    timeCalendarStart = '1970-01-01'

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        ds = Dataset(filename)
        scale_factor = (np.sin(0.9330127018922193) + 1) / 2

        timestamp = date if date else self.get_date(filename)
        self.srcDSProjection = NSR('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=70 +a=6371000 +b=6371000 +units=m +no_defs')
        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)
        self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('Computer')))
        self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform('MODELS')))
        self.dataset.SetMetadataItem('iso_category',
                                     json.dumps(pti.get_iso19115_topic_category('Oceans')))
        self.dataset.SetMetadataItem('Entry Title', str(ds.getncattr('title')))
        self.dataset.SetMetadataItem('gcmd_location',
                                     json.dumps(pti.get_gcmd_location('NORTH ATLANTIC OCEAN')))

    @staticmethod
    def get_date(filename):
        """Extract date and time parameters from filename and return
        it as a formatted string

        Parameters
        ----------

        filename: str
            nn

        Returns
        -------
            str, YYYY-mm-ddThh:MMZ

        Examples
        --------
            >>> Mapper.get_date('/path/to/NorKyst-800m_ZDEPTHS_his.an.2018110100.nc')
            '2018-11-01T00:00Z'
        """
        _, filename = os.path.split(filename)
        t = datetime.strptime(filename.split('.')[-2], '%Y%m%d%H')
        return datetime.strftime(t, '%Y-%m-%dT%H:%MZ')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes
