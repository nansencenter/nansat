from pyproj import CRS
from netCDF4 import Dataset

from nansat.nsr import NSR

import requests

import xml.etree.ElementTree as ET

from dateutil.parser import parse

import json
import pythesint as pti

from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.exceptions import WrongMapperError

class Mapper(NetcdfCF):

    def __init__(self, filename, gdal_dataset, gdal_metadata, file_num=0, *args, **kwargs):

        if not filename.endswith(".ncml"):
            raise WrongMapperError

        ds = Dataset(filename)
        for attr in ds.ncattrs():
            gdal_metadata[attr] = ds.getncattr(attr)

        if not 'meps' in gdal_metadata['title'].lower():
            raise WrongMapperError

        self.input_filename = filename
    
        self._create_empty_from_projection_variable(gdal_dataset, gdal_metadata)

        import ipdb
        ipdb.set_trace()
        # Add bands with metadata and corresponding values to the empty VRT
        self.create_bands(self._band_list(gdal_dataset, metadata, *args, **kwargs))

        #self.dataset.SetMetadataItem('time_coverage_start', parse(
        #    gdal_metadata['NC_GLOBAL#min_time'], ignoretz=True, fuzzy=True).isoformat())
        #self.dataset.SetMetadataItem('time_coverage_end', parse(
        #    gdal_metadata['NC_GLOBAL#max_time'], ignoretz=True, fuzzy=True).isoformat()))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('models')

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))


    def _create_empty_from_projection_variable(self, gdal_dataset, gdal_metadata):

        # TODO: only open dataset once or close them..
        ds = Dataset(self.input_filename)

        shape=[]
        for key in ds.dimensions.keys():
            shape.append(ds.dimensions[key].size)

        for var in ds.variables.keys():
            if ds[var].shape==tuple(shape):
                break

        grid_mapping = ds.variables[ds.variables[var].grid_mapping]

        grid_mapping_dict = {}
        for index in grid_mapping.ncattrs():
            grid_mapping_dict[index] = grid_mapping.getncattr(index)
        crs = CRS.from_cf(grid_mapping_dict)

        xx = ds.variables["x"][0:2]
        yy = ds.variables["y"][0:2]
        geotr = (xx[0], xx[1]-xx[0], 0, yy[0], 0, yy[1]-yy[0])

        self._init_from_dataset_params(
                    x_size = ds.dimensions["x"].size,
                    y_size = ds.dimensions["y"].size,
                    geo_transform = geotr,
                    projection = NSR(crs.to_proj4()).wkt,
                    metadata = gdal_metadata)
