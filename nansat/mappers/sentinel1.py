import json
import numpy as np
try:
    import scipy
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

import gdal
from netCDF4 import Dataset

import pythesint as pti

from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.utils import initial_bearing
from nansat.exceptions import WrongMapperError, NansatReadError

class Sentinel1(VRT):
    """ Mapper class to access Sentinel-1 data with netCDF4 to work for both opendap streams and
    local files.
    """
    timeVarName = 'time'
    input_filename = ''

    def __init__(self, filename, flip_gcp_line=False):
        #if not 'S1A' in filename or not 'S1B' in filename:
        #    raise WrongMapperError('%s: Not Sentinel 1A or 1B' %filename)
        if not self.dataset.GetMetadataItem('SATELLITE_IDENTIFIER') or \
                not self.dataset.GetMetadataItem('SATELLITE_IDENTIFIER').lower()=='sentinel-1':
            raise WrongMapperError('%s: Not Sentinel 1A or 1B' %filename)
        if not IMPORT_SCIPY:
            raise NansatReadError('Sentinel-1 data cannot be read because scipy is not installed')


        self.input_filename = filename
        try:
            self.ds = Dataset(filename)
        except OSError:
            self.ds = Dataset(filename+'#fillmismatch')
            self.input_filename = filename+'#fillmismatch'
        try:
            lon = self.ds.variables['GCP_longitude_'+self.ds.polarisation[:2]]
        except (AttributeError, KeyError):
            raise WrongMapperError('%s: Not Sentinel 1A or 1B' %filename)

        self._remove_geotransform()
        self._remove_geolocation()
        self.dataset.SetProjection('')
        self.dataset.SetGCPs(self.get_gcps(flip_gcp_line=flip_gcp_line), NSR().wkt)
        self.add_incidence_angle_band()
        self.add_look_direction_band()
        self.set_gcmd_dif_keywords()

    def set_gcmd_dif_keywords(self):
        mditem = 'entry_title'
        if not self.dataset.GetMetadataItem(mditem):
            self.dataset.SetMetadataItem(mditem, self.input_filename)
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

        self.dataset.SetMetadataItem('time_coverage_start',
                self.dataset.GetMetadataItem('ACQUISITION_START_TIME'))
        self.dataset.SetMetadataItem('time_coverage_end', 
        self.dataset.GetMetadataItem('ACQUISITION_STOP_TIME'))

    def get_gcps(self, flip_gcp_line=False):
        """ Get Ground Control Points for the dataset. 

        Note that OPeNDAP streams and netCDF files are read differently by gdal. The OPeNDAP streams
        are read by specifying the get parameters to the OPeNDAP url. The get parameters specify the
        reference dimensions, e.g., x and y. Since these are specified, the raster data is correctly
        referenced to the GCPs. However, when gdal reads a raster band from netCDF, it reads it
        "blindly". This is risky, since the definition of origo may be different in gdal vs the
        original data (e.g., first line starts in upper left corner or in lower left corner). For
        Sentinel-1, the raster data is flipped in relation to the GCPs, so we need to flip the GCP
        line vector as well.

        """
        lon = self.ds.variables['GCP_longitude_'+self.ds.polarisation[:2]][:]
        lat = self.ds.variables['GCP_latitude_'+self.ds.polarisation[:2]][:]
        line = self.ds.variables['GCP_line_'+self.ds.polarisation[:2]][:]
        if flip_gcp_line:
            # Flip line vector
            line = self.ds.dimensions['y'].size - line
        pixel = self.ds.variables['GCP_pixel_'+self.ds.polarisation[:2]][:]

        gcps = []
        for i0 in range(0, self.ds.dimensions['gcp_index'].size):
            gcp = gdal.GCP(float(lon[i0]), float(lat[i0]), 0, float(pixel[i0]),
                    float(line[i0]))
            gcps.append(gcp)

        return gcps

    def add_incidence_angle_band(self):
        # Get GCP variables
        pixel = self.ds['GCP_pixel_'+self.ds.polarisation[:2]][:].data
        line = self.ds['GCP_line_'+self.ds.polarisation[:2]][:].data
        inci = self.ds['GCP_incidenceAngle_'+self.ds.polarisation[:2]][:].data
        inci = inci.reshape(np.unique(line[:].data).shape[0],
                np.unique(pixel[:].data).shape[0])

        # Add incidence angle band
        inciVRT = VRT.from_array(inci)
        inciVRT = inciVRT.get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 1)
        self.band_vrts['inciVRT'] = inciVRT
        src = {'SourceFilename': self.band_vrts['inciVRT'].filename,
               'SourceBand': 1}
        dst = {'wkv': 'angle_of_incidence',
               'name': 'incidence_angle'}
        self.create_band(src, dst)
        self.dataset.FlushCache()

    def get_full_size_GCPs(self):
        # Get GCP variables
        pixel = self.ds['GCP_pixel_' + self.ds.polarisation[:2]][:].data
        line = self.ds['GCP_line_' + self.ds.polarisation[:2]][:].data
        lon = self.ds['GCP_longitude_' + self.ds.polarisation[:2]][:].data
        lat = self.ds['GCP_latitude_' + self.ds.polarisation[:2]][:].data
        lon = lon.reshape(np.unique(line[:].data).shape[0],
                np.unique(pixel[:].data).shape[0])
        lat = lat.reshape(np.unique(line[:].data).shape[0],
                np.unique(pixel[:].data).shape[0])
        return lon, lat

    def add_look_direction_band(self):
        lon, lat = self.get_full_size_GCPs()
        """
        TODO: Also in mapper_sentinel1_l1.py... Use this code there..
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
