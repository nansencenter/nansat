import json
import pythesint as pti

from osgeo import gdal
from netCDF4 import Dataset

from nansat.exceptions import WrongMapperError
from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF


class Mapper(NetcdfCF):

    def __init__(self, url, gdal_dataset, gdal_metadata, file_num=0, *args, **kwargs):

        if not url.endswith(".nc"):
            raise WrongMapperError

        ds = Dataset(url)
        for attr in ds.ncattrs():
            gdal_metadata[attr] = ds.getncattr(attr)

        if not 'meps' in gdal_metadata['title'].lower():
            raise WrongMapperError

        gdal_dataset = gdal.Open(url)
        metadata = gdal_dataset.GetMetadata()

        super(Mapper, self).__init__(url, gdal_dataset, gdal_metadata, *args, **kwargs)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('models')

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
