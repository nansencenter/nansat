# Name:         mapper_opendap_sentinel2.py
# Purpose:      Nansat mapping for ESA Sentinel-2 data from the Norwegian ground segment
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
from nansat.nsr import NSR


class Mapper(Opendap):

    baseURLs = [
            'http://nbstds.met.no/thredds/dodsC/NBS/S2A',
            'http://nbstds.met.no/thredds/dodsC/NBS/S2B',
    ]

    timeVarName = 'time'
    xName = 'x'
    yName = 'y'
    timeCalendarStart = '1981-01-01'
    GCP_STEP=100

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

        self.dataset.SetMetadataItem('entry_title', str(self.ds.getncattr('product_id')))
        self.dataset.SetMetadataItem('data_center', json.dumps(pti.get_gcmd_provider('NO/MET')))
        self.dataset.SetMetadataItem('ISO_topic_category',
                pti.get_iso19115_topic_category('Imagery/Base Maps/Earth Cover')['iso_topic_category'])

        mm = pti.get_gcmd_instrument('multi-spectral')
        ee = pti.get_gcmd_platform('sentinel-2')
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

    def get_gcps(self):

        lon_grid = self.ds.variables['lon']
        lat_grid = self.ds.variables['lat']

        dx = .5
        dy = .5
        gcps = []
        k = 0
        maxY = 0
        minY = lat_grid.shape[0]
        for i0 in range(0, lat_grid.shape[0], self.GCP_STEP):
            for i1 in range(0, lat_grid.shape[1], self.GCP_STEP):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(lon_grid[i0, i1])
                lat = float(lat_grid[i0, i1])
                #if (lon >= -180 and
                #    lon <= 180 and
                #    lat >= MIN_LAT and
                #    lat <= MAX_LAT):
                gcp = gdal.GCP(lon, lat, 0, i1 + dx, i0 + dy)
                gcps.append(gcp)
                k += 1
                maxY = max(maxY, i0)
                minY = min(minY, i0)
        yOff = minY
        ySize = maxY - minY

        # remove Y-offset from gcps
        for gcp in gcps:
            gcp.GCPLine -= yOff

        return gcps

