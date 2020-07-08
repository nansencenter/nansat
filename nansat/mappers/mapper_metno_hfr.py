# Name:        mapper_metno_hfr
# Purpose:     Mapper for CODAR SeaSonde High-Frequency radar data provided by MET Norway
#              For data see: https://thredds.met.no/thredds/catalog/remotesensinghfradar/catalog.html
# Authors:     Artem Moiseev
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html

import os
from datetime import datetime, timedelta
import json
from netCDF4 import Dataset, num2date
import numpy as np
from osgeo import gdal, osr
import pythesint as pti
from scipy.interpolate import griddata
from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError


class Mapper(VRT):
    
    BAND_NAMES = ['direction', 'ersc', 'ertc', 'espc', 'etmp', 'maxv',
                  'minv', 'sprc', 'u', 'v', 'velo', 'vflg', 'xdst', 'ydst']
    SUPPORTED_LOCATIONS = ['RDLm_TORU', 'RDLm_FRUH', 'RDLm_BERL']

    def __init__(self, filename, gdal_dataset, gdal_metadata, GCP_COUNT=10, timestamp=None, **kwargs):
        filename_name = os.path.split(filename)[-1].split('.')[0]
        # Check if correct mapper
        correct_mapper = False 
        for location in self.SUPPORTED_LOCATIONS:
            # If it matches with one of locateions break the loop and flag True
            if filename_name.startswith(location):
                correct_mapper = True
                break
        if not correct_mapper:
            raise WrongMapperError

        # Import NetCDF4 dataset
        nc_dataset = Dataset(filename)
        # Define projection (depending on the HFR)
        if nc_dataset.getncattr('site') == 'TORU':
            proj4 = '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
            GRID_PX_SIZE = 1500 # Final raster px size in meters
        elif nc_dataset.getncattr('site') == 'FRUH':
            proj4 = '+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
            GRID_PX_SIZE = 5000 # Final raster px size in meters
        elif nc_dataset.getncattr('site') == 'BERL':
            proj4 = '+proj=utm +zone=35 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
            GRID_PX_SIZE = 5000 # Final raster px size in meters
        else:
            raise WrongMapperError

        srs = osr.SpatialReference()
        srs.ImportFromProj4(proj4)
        projection = srs.ExportToWkt()
        # Get x grid and y grid
        x_grd, y_grd = self.create_linear_grid(nc_dataset['x'][:], nc_dataset['y'][:], GRID_PX_SIZE)
        raster_x_size, raster_y_size = x_grd.shape
        # Define geotransform
        geotransform = (x_grd.min(), GRID_PX_SIZE, 0.0, y_grd.max(), 0.0, GRID_PX_SIZE * -1)
        # Define x and y size
        self._init_from_dataset_params(raster_x_size, raster_y_size, geotransform, projection)
        # If required timestamp was not specified then extract date from filename and use first time
        if timestamp is None:
            timestamp = self.date_from_filename(filename)
        # Comvert time info from the dataset to the datetime
        timestamps = num2date(nc_dataset['time'][:].data, nc_dataset['time'].units)
        # find band id for the required timestamp
        # Note add 1 because in gdal counting starts from 1 not from 0
        src_timestamp_id = np.where(timestamps == timestamp)[0][0] + 1
        # Iterate through all subdatasets and bands to the dataset
        for subdataset in gdal_dataset.GetSubDatasets():
            # Get name of subdataset
            subdataset_name = subdataset[0].split(':')[2]
            # Check if the subdataset in the accepted 3D vars list
            if subdataset_name not in self.BAND_NAMES:
                continue
            gdal_subdataset = gdal.Open(subdataset[0])
            # need to be float for the nan replasement
            band_data = gdal_subdataset.GetRasterBand(int(src_timestamp_id)).ReadAsArray().astype('float')
            # remove fill value (replace with nan)
            fill_value = int(gdal_subdataset.GetMetadata_Dict()['#'.join([subdataset_name, '_FillValue'])])
            band_data[band_data == fill_value] = np.nan
            # Interpolate data on the regular grid
            band_grid_data = self.band2grid((nc_dataset['x'][:], nc_dataset['y'][:]), 
                                             band_data, (x_grd, y_grd))
            # Create VRT ffor the regridded data
            band_vrt = VRT.from_array(band_grid_data)
            # Add VRT to the list of all dataset vrts
            self.band_vrts[subdataset_name + 'VRT'] = band_vrt
            # Add band to the dataset
            src = {'SourceFilename': self.band_vrts[subdataset_name + 'VRT'].filename, 
                   'SourceBand': 1}
            # Add band specific metadata
            dst = {'name': subdataset_name}
            for key in gdal_subdataset.GetMetadata_Dict().keys():
                if key.startswith(subdataset_name):
                    clean_metadata_name = key.split('#')[1]
                    dst[clean_metadata_name] = gdal_subdataset.GetMetadata_Dict()[key]
            # Create band
            self.create_band(src, dst)
            self.dataset.FlushCache()

        # Set GCMD metadata
        self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('SCR-HF')))
        self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform('CODAR SeaSonde')))
        self.dataset.SetMetadataItem('Data Center', json.dumps(pti.get_gcmd_provider('NO/MET')))
        self.dataset.SetMetadataItem('Entry Title', 'Near-Real Time Surface Ocean Radial Velocity')
        self.dataset.SetMetadataItem('gcmd_location',json.dumps(pti.get_gcmd_location('NORTH SEA')))
        # Set time coverage metadata
        self.dataset.SetMetadataItem('time_coverage_start', timestamp.isoformat())
        self.dataset.SetMetadataItem('time_coverage_end',
                                     (timestamp + timedelta(minutes=59, seconds=59)).isoformat())
        # Set NetCDF dataset metadata
        for key, value in gdal_dataset.GetMetadata_Dict().items():
            self.dataset.SetMetadataItem(key.split('#')[1], value)

    def date_from_filename(self, src):
        filename = os.path.splitext(os.path.basename(src))[0]
        year, month, day = filename.split('_')[-3:]
        src_timestamp = datetime(int(year), int(month), int(day), 0, 0, 0)
        return src_timestamp

    def create_linear_grid(self, x, y, px_size):
        x_grd, y_grd = np.meshgrid(np.arange(x.min(), x.max(), px_size),
                                   np.arange(y.max(), y.min(), px_size * -1))
        return x_grd, y_grd

    def band2grid(self, src_grd, var, dst_grd):
        # Points [(x, y), ... , ] from original file
        points = list(zip(src_grd[0].flatten(), src_grd[1].flatten()))
        repr_var = griddata(points, var.flatten(), (dst_grd[0], dst_grd[1]), method='nearest')
        return repr_var
