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

try:
    import scipy
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

import pythesint as pti

from nansat.mappers.sentinel1 import Sentinel1
from nansat.mappers.opendap import Opendap
from nansat.vrt import VRT
from nansat.nsr import NSR
from nansat.utils import initial_bearing


class Mapper(Opendap, Sentinel1):

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

        if not IMPORT_SCIPY:
            raise NansatReadError('Sentinel-1 data cannot be read because scipy is not installed')

        timestamp = date if date else self.get_date(filename)

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

        Sentinel1.__init__(self, filename)
        self.add_calibrated_nrcs(filename)
        self.add_nrcs_VV_from_HH(filename)

    def add_calibrated_nrcs(self, filename):
        layer_time_id, layer_date = Opendap.get_layer_datetime(None,
                self.convert_dstime_datetimes(self.get_dataset_time()))
        polarizations = [self.ds.polarisation[i:i+2] for i in range(0,len(self.ds.polarisation),2)]
        for pol in polarizations:
            dims = list(self.ds.variables['dn_%s' %pol].dimensions)
            dims[dims.index(self.timeVarName)] = layer_time_id
            src = [
                    self.get_metaitem(filename, 'Amplitude_%s' %pol, dims)['src'],
                    self.get_metaitem(filename, 'sigmaNought_%s' %pol, dims)['src']
                ]
            dst = {
                    'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                    'PixelFunctionType': 'Sentinel1Calibration',
                    'polarization': pol,
                    'suffix': pol,
                }
            self.create_band(src, dst)
            self.dataset.FlushCache()

    def add_nrcs_VV_from_HH(self, filename):
        if not 'Amplitude_HH' in self.ds.variables.keys():
            return
        layer_time_id, layer_date = Opendap.get_layer_datetime(None,
                self.convert_dstime_datetimes(self.get_dataset_time()))
        dims = list(self.ds.variables['dn_HH'].dimensions)
        dims[dims.index(self.timeVarName)] = layer_time_id
        src = [
                self.get_metaitem(filename, 'Amplitude_HH', dims)['src'],
                self.get_metaitem(filename, 'sigmaNought_HH', dims)['src'],
                {'SourceFilename': self.band_vrts['inciVRT'].filename, 'SourceBand': 1}
            ]
        dst = {
                'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                'PixelFunctionType': 'Sentinel1Sigma0HHToSigma0VV',
                'polarization': 'VV',
                'suffix': 'VV'}
        self.create_band(src, dst)
        self.dataset.FlushCache()

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

    def get_geotransform(self):
        """ Return fake and temporary geotransform. This will be replaced by gcps in
        Sentinel1.__init__ 
        """
        xx = self.ds.variables['lon'][0:100:50, 0].data
        yy = self.ds.variables['lat'][0, 0:100:50].data
        return xx[0], xx[1]-xx[0], 0, yy[0], 0, yy[1]-yy[0]

