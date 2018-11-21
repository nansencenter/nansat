# Name:         mapper_opendap_sentinel1.py
# Purpose:      Nansat mapping for ESA Sentinel-1 data from the Norwegian ground segment
# Author:       Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
from datetime import datetime
import json
import warnings

import numpy as np
from netCDF4 import Dataset
import gdal

import pythesint as pti

from nansat.mappers.opendap import Opendap
from nansat.vrt import VRT
from nansat.nsr import NSR


class Mapper(Opendap):

    baseURLs = [
            'http://nbstds.met.no/thredds/dodsC/NBS/S1A',
            'http://nbstds.met.no/thredds/dodsC/NBS/S1B',
    ]

    timeVarName = 'time'
    xName = 'x'
    yName = 'y'
    timeCalendarStart = '1981-01-01'
    srcDSProjection = NSR().wkt

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

        self._remove_geotransform()
        self._remove_geolocation()
        self.dataset.SetGCPs(self.get_gcps(), NSR().wkt)

        self.dataset.SetMetadataItem('entry_title', filename)
        self.dataset.SetMetadataItem('data_center', json.dumps(pti.get_gcmd_provider('NO/MET')))
        self.dataset.SetMetadataItem('ISO_topic_category',
                pti.get_iso19115_topic_category('Imagery/Base Maps/Earth Cover')['iso_topic_category'])

        mm = pti.get_gcmd_instrument('sar')
        if self.ds.MISSION_ID=='S1A':
            ee = pti.get_gcmd_platform('sentinel-1a')
        else:
            ee = pti.get_gcmd_platform('sentinel-1b')
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
        t = datetime.strptime(filename.split('_')[4], '%Y%m%dT%H%M%S')
        return datetime.strftime(t, '%Y-%m-%dT%H:%M:%SZ')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes

    def get_gcps(self):

        lon = self.ds.variables['GCP_longitude_'+self.ds.polarisation[:2]]
        lat = self.ds.variables['GCP_latitude_'+self.ds.polarisation[:2]]
        line = self.ds.variables['GCP_line_'+self.ds.polarisation[:2]]
        pixel = self.ds.variables['GCP_pixel_'+self.ds.polarisation[:2]]

        gcps = []
        for i0 in range(0, self.ds.dimensions['gcp_index'].size):
            gcp = gdal.GCP(float(lon[i0].data), float(lat[i0].data), 0, float(pixel[i0].data),
                    float(line[i0].data))
            gcps.append(gcp)

        return gcps

    def get_geotransform(self):
        """ Return fake and temporary geotransform. This will be replaced by gcps at the end of
        __init__ 
        """
        xx = self.ds.variables['lon'][0:100:50, 0].data
        yy = self.ds.variables['lat'][0, 0:100:50].data
        return xx[0], xx[1]-xx[0], 0, yy[0], 0, yy[1]-yy[0]

