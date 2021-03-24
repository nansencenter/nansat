# Name:         mapper_opendap_mywave.py
# Purpose:      Nansat mapping for MyWaveWAV data provided by MET.NO
# Author:       Artem Moiseev
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from nansat.mappers.opendap import Opendap
from nansat.exceptions import WrongMapperError
from nansat.nsr import NSR
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import json
import pythesint as pti
import os


class Mapper(Opendap):

    baseURLs = ['http://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive', 
                'https://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be']

    timeVarName = 'time'
    xName = 'rlon'
    yName = 'rlat'
    timeCalendarStart = '1970-01-01'

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)
        ds = Dataset(filename)
        try:
            self.srcDSProjection = NSR(ds.variables['projection_3'].proj4 +
                                       ' +to_meter=0.0174532925199 +wktext')
        except KeyError:
            raise WrongMapperError

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

        self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('Computer')))
        self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform('MODELS')))
        self.dataset.SetMetadataItem('Data Center', json.dumps(pti.get_gcmd_provider('NO/MET')))
        self.dataset.SetMetadataItem('Entry Title', str(ds.getncattr('title')))
        self.dataset.SetMetadataItem('Entry Title',
                                     json.dumps(pti.get_iso19115_topic_category('Oceans')))
        self.dataset.SetMetadataItem('gcmd_location',
                                     json.dumps(pti.get_gcmd_location('sea surface')))

    @staticmethod
    def get_date(filename):
        """Extract date and time parameters from filename and return
        it as a formatted (isoformat) string

        Parameters
        ----------

        filename: str

        Returns
        -------
            str, YYYY-mm-ddThh:MM:00Z

        Examples
        --------
            >>> Mapper.get_date('/path/to/MyWave_wam4_WAVE_20171029T18Z.nc')
            '2017-10-29T18:00:00Z'
        """
        _, filename = os.path.split(filename)
        t = datetime.strptime(filename.split('_')[-1], '%Y%m%dT%HZ.nc')
        return datetime.strftime(t, '%Y-%m-%dT%H:%M:00Z')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes
