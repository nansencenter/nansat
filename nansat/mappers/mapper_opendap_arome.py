# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
import json
import pytz
import datetime

import numpy as np
import pythesint as pti

from netCDF4 import Dataset
from osgeo import gdal
from pyproj import CRS

from nansat.exceptions import WrongMapperError
from nansat.nsr import NSR
from nansat.vrt import VRT

from nansat.mappers.mapper_netcdf_cf import Mapper as MapperNCCF
#from nansat.mappers.mapper_arome import Mapper as MapperArome
from nansat.mappers.opendap import Opendap


class Mapper(MapperNCCF, Opendap):

    baseURLs = ['http://thredds.met.no/thredds/catalog/arome25/catalog.html',
                'https://thredds.met.no/thredds/dodsC/aromearcticarchive',
                'http://thredds.met.no/thredds/dodsC/aromearcticarchive',
                'https://thredds.met.no/thredds/dodsC/meps25epsarchive',
                'http://thredds.met.no/thredds/dodsC/meps25epsarchive']
    timeVarName = 'time'

    def __init__(self, filename, gdal_dataset, gdal_metadata, netcdf_dim=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        if netcdf_dim is None:
            netcdf_dim = {}
        
        self.input_filename = filename

        if gdal_dataset is None:
            gdal_dataset = gdal.Open(filename)

        try:
            ds = Dataset(filename)
        except OSError:
            raise WrongMapperError

        metadata = {}
        for attr in ds.ncattrs():
            metadata[attr] = self._fix_encoding(ds.getncattr(attr))

        if 'title' not in metadata.keys() or 'arome' not in metadata['title'].lower():
            raise WrongMapperError

        xsize = ds.dimensions['x'].size
        ysize = ds.dimensions['y'].size

        # Pick 10 meter height dimension only
        height_dim = 'height7'
        if height_dim not in ds.dimensions.keys():
            raise WrongMapperError
        if ds.dimensions[height_dim].size != 1:
            raise WrongMapperError
        if ds.variables[height_dim][0].data != 10:
            raise WrongMapperError

        varnames = []
        for var in ds.variables:
            var_dimensions = ds.variables[var].dimensions
            if var_dimensions == ('time', height_dim, 'y', 'x'):
                varnames.append(var)

        # Projection
        try:
            grid_mapping = ds.variables[ds.variables[varnames[0]].grid_mapping]
        except KeyError:
            raise WrongMapperError

        grid_mapping_dict = {}
        for index in grid_mapping.ncattrs():
            grid_mapping_dict[index] = grid_mapping.getncattr(index)
        crs = CRS.from_cf(grid_mapping_dict)
        nsr = NSR(crs.to_proj4())

        # Geotransform
        xx = ds.variables['x'][0:2]
        yy = ds.variables['y'][0:2]
        gtrans = xx[0], xx[1]-xx[0], 0, yy[0], 0, yy[1]-yy[0]

        self._init_from_dataset_params(xsize, ysize, gtrans, nsr.wkt)

        meta_dict = []
        if bands is None:
            bands = varnames
        for band in bands:
            if band not in ds.variables.keys():
                continue
            dimension_names, dim_sizes = self._get_dimension_info(band)
            self._pop_spatial_dimensions(dimension_names)
            index = self._get_index_of_dimensions(dimension_names, netcdf_dim, dim_sizes)
            fn = self._get_sub_filename(filename, band, dim_sizes, index)
            band_metadata = self.get_band_metadata_dict(fn, ds.variables[band])
            # Add time stamp to band metadata
            tt = datetime.datetime.fromisoformat(str(self.times()[index["time"]["index"]]))
            if tt.tzinfo is None:
                tt = pytz.utc.localize(tt)
            band_metadata["dst"]["time"] = tt.isoformat()
            meta_dict.append(band_metadata)

        self.create_bands(meta_dict)

        # Copy metadata
        for key in metadata.keys():
            self.dataset.SetMetadataItem(str(key), str(metadata[key]))

        mm = pti.get_gcmd_instrument('Computer')
        ee = pti.get_gcmd_platform('ecmwfifs')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    @staticmethod
    def get_date(filename):
        """Extract date and time parameters from filename and return
        it as a datetime object.
        """
        _, filename = os.path.split(filename)
        return datetime.datetime.strptime(filename.split('_')[-1], '%Y%m%dT%HZ.nc')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.epoch).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes
