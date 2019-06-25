# Name:         mapper_arome.py
# Purpose:      Nansat mapping for AROME-Arctic and MEPS (MetCoOp Ensemble
#               Prediction System) data provided by MET.NO
# Author:       Artem Moiseev
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from nansat.mappers.mapper_arome import Mapper as MapperArome
from nansat.mappers.opendap import Opendap
from nansat.exceptions import WrongMapperError
from nansat.nsr import NSR
import pythesint as pti
import os
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import json


class Mapper(Opendap, MapperArome):

    baseURLs = ['http://thredds.met.no/thredds/catalog/arome25/catalog.html',
                'https://thredds.met.no/thredds/dodsC/aromearcticarchive',
                'http://thredds.met.no/thredds/dodsC/aromearcticarchive',
                'https://thredds.met.no/thredds/dodsC/meps25epsarchive',
                'http://thredds.met.no/thredds/dodsC/meps25epsarchive']
    timeVarName = 'time'
    xName = 'x'
    yName = 'y'
    timeCalendarStart = '1970-01-01'

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)
        ds = Dataset(filename)

        try:
            self.srcDSProjection = NSR(ds.variables['projection_lambert'].proj4)
        except KeyError:
            raise WrongMapperError

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

        mm = pti.get_gcmd_instrument('Computer')
        ee = pti.get_gcmd_platform('ecmwfifs')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

        md_item = 'Data Center'
        if not self.dataset.GetMetadataItem(md_item):
            self.dataset.SetMetadataItem(md_item, 'NO/MET')
        md_item = 'Entry Title'
        if not self.dataset.GetMetadataItem(md_item):
            self.dataset.SetMetadataItem(md_item, str(ds.getncattr('title')))
        md_item = 'summary'
        if not self.dataset.GetMetadataItem(md_item):
            summary = """
            AROME_Arctic is a convection-permitting atmosphere model covering parts of the Barents
            Sea and the Nordic Arctic. It has horizontal resolution of 2.5 km and 65 vertical
            levels. AROME_Arctic runs for 66 hours four times a day (00,06,12,18) with three-hourly
            cycling for data assimilation. Boundary data is from ECMWF. Model code based on HARMONIE
            cy40h1.1
            """
            self.dataset.SetMetadataItem(md_item, str(summary))

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
            >>> Mapper.get_date('/path/to/arome_arctic_full_2_5km_20171030T21Z.nc')
            '2017-10-30T21:00Z'
        """
        _, filename = os.path.split(filename)
        t = datetime.strptime(filename.split('_')[-1], '%Y%m%dT%HZ.nc')
        return datetime.strftime(t, '%Y-%m-%dT%H:%MZ')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes
