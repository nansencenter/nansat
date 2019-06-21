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

from nansat.mappers.opendap import Opendap
from nansat.vrt import VRT
from nansat.nsr import NSR
from nansat.tools import initial_bearing


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

        if not IMPORT_SCIPY:
            raise NansatReadError('Sentinel-1 data cannot be read because scipy is not installed')

        timestamp = date if date else self.get_date(filename)

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

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

        # Add look direction band
        pixel = self.ds['GCP_pixel_%s' %polarizations[0]][:].data
        line = self.ds['GCP_line_%s' %polarizations[0]][:].data
        lon = self.ds['GCP_longitude_%s' %polarizations[0]][:].data
        lat = self.ds['GCP_latitude_%s' %polarizations[0]][:].data
        lon = lon.reshape(np.unique(line[:].data).shape[0],
                np.unique(pixel[:].data).shape[0])
        lat = lat.reshape(np.unique(line[:].data).shape[0],
                np.unique(pixel[:].data).shape[0])

        """
        TODO: Repetition from mapper_sentinel1_l1.py... Consider generalising
        """
        sat_heading = initial_bearing(lon[:-1, :], lat[:-1, :], lon[1:, :], lat[1:, :])
        look_direction = scipy.ndimage.interpolation.zoom( np.mod(sat_heading + 90, 360),
                (np.shape(lon)[0] / (np.shape(lon)[0]-1.), 1))

        # Decompose, to avoid interpolation errors around 0 <-> 360
        look_direction_u = np.sin(np.deg2rad(look_direction))
        look_direction_v = np.cos(np.deg2rad(look_direction))
        look_u_VRT = VRT.from_array(look_direction_u)
        look_v_VRT = VRT.from_array(look_direction_v)
        lookVRT = VRT.from_lonlat(lon, lat)
        lookVRT.create_band([{'SourceFilename': look_u_VRT.filename,
                               'SourceBand': 1},
                              {'SourceFilename': look_v_VRT.filename,
                               'SourceBand': 1}],
                             {'PixelFunctionType': 'UVToDirectionTo'}
                             )

        # Blow up to full size
        lookVRT = lookVRT.get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 1)

        # Store VRTs so that they are accessible later
        self.band_vrts['look_u_VRT'] = look_u_VRT
        self.band_vrts['look_v_VRT'] = look_v_VRT
        self.band_vrts['lookVRT'] = lookVRT

        src = {
                'SourceFilename': self.band_vrts['lookVRT'].filename,
                'SourceBand': 1
            }
        dst = {
                'wkv': 'sensor_azimuth_angle',
                'name': 'look_direction'
            }
        self.create_band(src, dst)
        self.dataset.FlushCache()
        """ End repetition """

        self._remove_geotransform()
        self._remove_geolocation()
        self.dataset.SetProjection('')
        self.dataset.SetGCPs(self.get_gcps(), NSR().wkt)

        mditem = 'entry_title'
        if not self.dataset.GetMetadataItem(mditem):
            self.dataset.SetMetadataItem(mditem, filename)
        mditem = 'data_center'
        if not self.dataset.GetMetadataItem(mditem):
            self.dataset.SetMetadataItem(mditem, json.dumps(pti.get_gcmd_provider('NO/MET')))
        mditem = 'ISO_topic_category'
        if not self.dataset.GetMetadataItem(mditem):
            self.dataset.SetMetadataItem(mditem,
                pti.get_iso19115_topic_category('Imagery/Base Maps/Earth Cover')['iso_topic_category'])

        mm = pti.get_gcmd_instrument('sar')
        if self.ds.MISSION_ID=='S1A':
            ee = pti.get_gcmd_platform('sentinel-1a')
        else:
            ee = pti.get_gcmd_platform('sentinel-1b')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

        # Times in opendap.py are wrong for S1...
        self.dataset.SetMetadataItem('time_coverage_start',
                self.dataset.GetMetadataItem('ACQUISITION_START_TIME'))
        self.dataset.SetMetadataItem('time_coverage_end', 
                self.dataset.GetMetadataItem('ACQUISITION_STOP_TIME'))

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

