import json
import pythesint as pti

from osgeo import gdal
from pyproj import CRS
from netCDF4 import Dataset

from nansat.nsr import NSR
from nansat.exceptions import WrongMapperError
from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.mappers.opendap import Opendap


class Mapper(NetcdfCF, Opendap):

    def __init__(self, url, gdal_dataset, gdal_metadata, file_num=0, bands=None, *args, **kwargs):

        if not url.endswith(".nc"):
            raise WrongMapperError

        try:
            ds = Dataset(url)
        except OSError:
            raise WrongMapperError

        metadata = {}
        for attr in ds.ncattrs():
            content = ds.getncattr(attr)
            content.replace("æ", "ae")
            content.replace("ø", "oe")
            content.replace("å", "aa")
            metadata[attr] = content

        if 'title' not in metadata.keys() or 'meps' not in metadata['title'].lower():
            raise WrongMapperError

        self.input_filename = url

        xsize = ds.dimensions['x'].size
        ysize = ds.dimensions['y'].size

        # Pick 10 meter height dimension only
        height_dim = 'height6'
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
            index = self._get_index_of_dimensions(dimension_names, {}, dim_sizes)
            fn = self._get_sub_filename(url, band, dim_sizes, index)
            meta_dict.append(self.get_band_metadata_dict(fn, ds.variables[band]))

        self.create_bands(meta_dict)

        # Copy metadata
        for key in metadata.keys():
            self.dataset.SetMetadataItem(str(key), str(metadata[key]))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('models')

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
