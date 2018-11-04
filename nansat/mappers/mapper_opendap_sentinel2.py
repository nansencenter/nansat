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
            'http://nbstds.met.no/thredds/dodsC/NBS/S1B',
            'http://nbstds.met.no/thredds/dodsC/NBS/S2B',
    ]

    timeVarName = 'time'
    xName = 'x'
    yName = 'y'
    timeCalendarStart = '1981-01-01'
    srcDSProjection = NSR().wkt
    GCP_STEP=100

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)

        if date is None:
            warnings.warn('Date is not specified! Will return the first layer. '
                          'Please add date="YYYY-MM-DD"')

        # TODO: <self.filename> will be changed to vrt filename after init vrt
        self.filename = filename
        self.cachedir = cachedir
        self.ds = self.get_dataset(ds)

        ds_time = self.get_dataset_time()
        ds_times = self.convert_dstime_datetimes(ds_time)
        layer_time_id, layer_date = Opendap.get_layer_datetime(date, ds_times)

        if bands is None:
            var_names = self.get_geospatial_variable_names()
        else:
            var_names = bands
        # create VRT with correct lon/lat (geotransform)
        raster_x, raster_y = self.get_shape()
    
        import ipdb
        ipdb.set_trace()

        self._init_from_lonlat(self.ds.variables['lon'], self.ds.variables['lat'], add_gcps=True, **kwargs)

        meta_dict = self.create_metadict(filename, var_names, layer_time_id)

        self.create_bands(meta_dict)

        # set time
        time_res_sec = self.get_time_coverage_resolution()
        self.dataset.SetMetadataItem('time_coverage_start', str(layer_date))
        self.dataset.SetMetadataItem('time_coverage_end', str(layer_date + time_res_sec))

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

        #self.dataset.SetMetadataItem('entry_title', str(ds.getncattr('title')))
        #self.dataset.SetMetadataItem('data_center', json.dumps(pti.get_gcmd_provider('UK/MOD/MET')))
        #self.dataset.SetMetadataItem('ISO_topic_category',
        #        pti.get_iso19115_topic_category('oceans')['iso_topic_category'])
        #self.dataset.SetMetadataItem('gcmd_location', json.dumps(pti.get_gcmd_location('sea surface')))

        #mm = pti.get_gcmd_instrument('amsr-e')
        #ee = pti.get_gcmd_platform('aqua')
        #self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        #self.dataset.SetMetadataItem('platform', json.dumps(ee))

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

