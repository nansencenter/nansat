# Name:         mapper_nora10_local_vpv.py
# Purpose:      Mapper for Nora10 NetCDF files
# Authors:      Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from datetime import datetime, timedelta
import numpy as np

from nansat.utils import gdal, ogr
from nansat.vrt import VRT

from nansat.exceptions import WrongMapperError

# Hardcoded for MET file sytem; may be overridden by other user with copy of data
baseFolder = '/vol/mis/dmf2/midlertidig/hindcast/hildeh/nora10/NetCDF/'
keywordBase = __name__[7:]


class Mapper(VRT):
    def __init__(self, filename, gdalDataset, gdalMetadata, logLevel=30,
                 **kwargs):
        if filename[0:len(keywordBase)] != keywordBase:
            raise WrongMapperError(__file__,
                                   "Not Nora10 data converted from felt to netCDF")

        requestedTime = datetime.strptime(filename[len(keywordBase)+1:],
                                          '%Y%m%d%H%M')
        # For correct rounding
        fileTime = requestedTime + timedelta(minutes=30)
        fileTime = fileTime - timedelta(minutes=fileTime.minute)

        nc_file = (baseFolder + 'windspeed_10m' +
                   fileTime.strftime('/%Y/%m/') + 'windspeed_' +
                   fileTime.strftime('%Y%m%d%H.nc'))
        nc_file_winddir = (baseFolder + 'winddir_10m' +
                           fileTime.strftime('/%Y/%m/') +
                           'winddir_' + fileTime.strftime('%Y%m%d%H.nc'))

        # Would prefer to use geotransform, but ob_tran
        # (General Oblique Transformation) is not supported by GDAL
        # Keeping lines below for potential future use:
        #proj4 = gdalDataset.GetMetadataItem('projection_rotated_ll#proj4')
        #proj4 = '+proj=ob_tran +o_proj=longlat +lon_0=-40 +o_lat_p=22 +a=6367470 +e=0'
        #rlatmin = -13.25; rlatmax  = 26.65; deltarlat = 0.1
        #rlonmin = 5.75; rlonmax  = 30.45; deltarlon = 0.1

        # Needed due to precence of time dimension in netCDF file
        gdal.SetConfigOption('GDAL_NETCDF_BOTTOMUP', 'No')

        # Read relevant arrays into memory
        g = gdal.Open('NETCDF:"' + nc_file + '":' + 'windspeed_10m')
        ws_10m = np.flipud(g.GetRasterBand(1).ReadAsArray())
        g = gdal.Open('NETCDF:"' + nc_file_winddir + '":' +
                        'wind_direction_10m')
        wd_10m = np.flipud(g.GetRasterBand(1).ReadAsArray())
        g = gdal.Open('NETCDF:"' + nc_file + '":' + 'latitude')
        lat = np.flipud(g.GetRasterBand(1).ReadAsArray())
        g = gdal.Open('NETCDF:"' + nc_file + '":' + 'longitude')
        lon = np.flipud(g.GetRasterBand(1).ReadAsArray())

        u10 = -ws_10m*np.sin(np.deg2rad(wd_10m))
        v10 = -ws_10m*np.cos(np.deg2rad(wd_10m))
        VRT_u10 = VRT(array=u10, lat=lat, lon=lon)
        VRT_v10 = VRT(array=v10, lat=lat, lon=lon)

        # Store band_vrts so that they are available after reprojection etc
        self.band_vrts = {'u_VRT': VRT_u10,
                        'v_VRT': VRT_v10}

        metaDict = []
        metaDict.append({'src': {'SourceFilename': VRT_u10.filename,
                                 'SourceBand': 1},
                         'dst': {'wkv': 'eastward_wind',
                                 'name': 'eastward_wind'}})
        metaDict.append({'src': {'SourceFilename': VRT_v10.filename,
                                 'SourceBand': 1},
                         'dst': {'wkv': 'northward_wind',
                                 'name': 'northward_wind'}})

        # Add pixel function with wind speed
        metaDict.append({
            'src': [{'SourceFilename': self.band_vrts['u_VRT'].filename,
                     'SourceBand': 1,
                     'DataType': 6},
                    {'SourceFilename': self.band_vrts['v_VRT'].filename,
                     'SourceBand': 1,
                     'DataType': 6}],
            'dst': {'wkv': 'wind_speed',
                    'name': 'windspeed',
                    'height': '10 m',
                    'PixelFunctionType': 'UVToMagnitude'}})

        # Add pixel function with wind direction
        metaDict.append({
            'src': [{'SourceFilename': self.band_vrts['u_VRT'].filename,
                     'SourceBand': 1,
                     'DataType': 6},
                    {'SourceFilename': self.band_vrts['v_VRT'].filename,
                     'SourceBand': 1,
                     'DataType': 6}],
            'dst': {'wkv': 'wind_from_direction',
                    'name': 'winddir',
                    'height': '10 m',
                    'PixelFunctionType': 'UVToDirectionFrom'}})

        # create empty VRT dataset with geolocation only
        self._init_from_lonlat(lon, lat)

        # add bands with metadata and corresponding values
        # to the empty VRT
        self.create_bands(metaDict)

        # Add time
        self.dataset.SetMetadataItem('time_coverage_start', fileTime.isoformat())
