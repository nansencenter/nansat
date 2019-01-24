import re
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

        super(Mapper, self).__init__(filename, gdal_dataset, metadata, quartile=quartile, *args, **kwargs)

        lat = self.dataset.GetRasterBand(self._latitude_band_number(gdal_dataset)).ReadAsArray()
        lon = self.dataset.GetRasterBand(self._longitude_band_number(gdal_dataset)).ReadAsArray()
        lon = ScatterometryMapper.shift_longitudes(lon)
        self.set_gcps(lon, lat, gdal_dataset)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('seawinds')
        ee = pti.get_gcmd_platform('quikscat')
        provider = metadata['NC_GLOBAL#institution']
        if provider.lower()=='jpl':
            provider = 'NASA/JPL/QUIKSCAT'
        provider = pti.get_gcmd_provider(provider)

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
        self.dataset.SetMetadataItem('data_center', json.dumps(provider))
        self.dataset.SetMetadataItem('entry_title', metadata['NC_GLOBAL#title'])
        self.dataset.SetMetadataItem('ISO_topic_category',
                json.dumps(pti.get_iso19115_topic_category('Oceans')))
