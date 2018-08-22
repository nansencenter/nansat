import json
import numpy as np
from datetime import datetime
import pythesint as pti

import gdal

from nansat.vrt import VRT
from nansat.geolocation import Geolocation
from nansat.nsr import NSR
from nansat.domain import Domain
#from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.mappers.scatterometers import Mapper as ScatterometryMapper
from nansat.exceptions import WrongMapperError

#class Mapper(NetcdfCF):
class Mapper(ScatterometryMapper):
    """ Nansat mapper for QuikScat """

    def __init__(self, filename, gdal_dataset, metadata, quartile=0, *args, **kwargs):

        if not 'quikscat' in metadata.get('NC_GLOBAL#source', '').lower():
            raise WrongMapperError

        super(Mapper, self).__init__(filename, gdal_dataset, metadata, *args, **kwargs)

        band_lat = self.dataset.GetRasterBand(self._latitude_band_number(gdal_dataset))
        # Check that it is actually latitudes
        if not band_lat.GetMetadata()['standard_name'] == 'latitude':
            raise ValueError('Cannot find latitude band')
        lat = band_lat.ReadAsArray()

        band_lon = self.dataset.GetRasterBand(self._longitude_band_number(gdal_dataset))
        # Check that it is actually longitudes
        if not band_lon.GetMetadata()['standard_name'] == 'longitude':
            raise ValueError('Cannot find longitude band')
        lon = band_lon.ReadAsArray()

        lon = ScatterometryMapper.shift_longitudes(lon)

        self.set_gcps(lon, lat, gdal_dataset)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('seawinds')
        ee = pti.get_gcmd_platform('quikscat')

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

