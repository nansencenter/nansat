import os
import json

import pythesint as pti

from osgeo import gdal
from pyproj import CRS
from netCDF4 import Dataset

from nansat.nsr import NSR
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

        url = self._get_odap_url(filename, file_num)
        gdal_dataset = gdal.Open(url)
        metadata = gdal_dataset.GetMetadata()
        import ipdb
        ipdb.set_trace()

        super(Mapper, self).__init__(url, gdal_dataset, gdal_metadata, *args, **kwargs)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('models')

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))


    def _get_odap_url(self, fn, file_num):
        url = os.path.split(fn)[0] + "/member_%02d" % int(os.path.basename(fn)[8:11]) + "/meps_" + os.path.basename(fn).split("_")[2] + "_%02d_" % file_num + os.path.basename(fn).split("_")[3][:-2]
        return url


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
