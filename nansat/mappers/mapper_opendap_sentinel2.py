# Name:         mapper_opendap_sentinel2.py
# Purpose:      Nansat mapping for ESA Sentinel-2 data from the Norwegian ground segment
# Author:       Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from nansat.mappers.opendap import Opendap
from nansat.nsr import NSR
import pythesint as pti
import os
from datetime import datetime
import numpy as np
import json
from netCDF4 import Dataset


class Mapper(Opendap):

    baseURLs = [
            'http://nbstds.met.no/thredds/dodsC/NBS/S1B',
            'http://nbstds.met.no/thredds/dodsC/NBS/S2B',
    ]

    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    timeCalendarStart = '1981-01-01'
    srcDSProjection = NSR().wkt

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)
        ds = Dataset(filename)
        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)
        self.dataset.SetMetadataItem('entry_title', str(ds.getncattr('title')))
        self.dataset.SetMetadataItem('data_center', json.dumps(pti.get_gcmd_provider('UK/MOD/MET')))
        self.dataset.SetMetadataItem('ISO_topic_category',
                pti.get_iso19115_topic_category('oceans')['iso_topic_category'])
        self.dataset.SetMetadataItem('gcmd_location', json.dumps(pti.get_gcmd_location('sea surface')))

        mm = pti.get_gcmd_instrument('amsr-e')
        ee = pti.get_gcmd_platform('aqua')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    @staticmethod
    def get_date(filename):
        """Extract date and time parameters from filename and return
        it as a formatted (isoformat) string

        Parameters
        ----------

        filename: str
            nn

        Returns
        -------
            str, YYYY-mm-ddThh:MMZ

        """
        _, filename = os.path.split(filename)
        t = datetime.strptime(filename.split('_')[2], '%Y%m%dT%H%M%S')
        return datetime.strftime(t, '%Y-%m-%dT%H:%M:%SZ')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes
