import re
import json
import numpy as np
from datetime import datetime
import pythesint as pti

from netCDF4 import Dataset

import gdal

from nansat.vrt import VRT
from nansat.geolocation import Geolocation
from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.mappers.scatterometers import Mapper as ScatterometryMapper
from nansat.exceptions import WrongMapperError

class Mapper(ScatterometryMapper):
    """ Nansat mapper for ASCAT """

    def __init__(self, filename, gdal_dataset, metadata, quartile=0, *args, **kwargs):

        if not metadata.has_key('NC_GLOBAL#source') \
                or not 'ascat' in metadata['NC_GLOBAL#source'].lower():
            raise WrongMapperError

        super(Mapper, self).__init__(filename, gdal_dataset, metadata, *args, **kwargs)
 
        band_lat = self.dataset.GetRasterBand(self._latitude_band_number(gdal_dataset))
        # Check that it is actually latitudes
        if not band_lat.GetMetadata()['long_name'] == 'latitude':
            raise ValueError('Cannot find latitude band')
        lat = band_lat.ReadAsArray()

        band_lon = self.dataset.GetRasterBand(self._longitude_band_number(gdal_dataset))
        # Check that it is actually longitudes
        if not band_lon.GetMetadata()['long_name'] == 'longitude':
            raise ValueError('Cannot find longitude band')
        lon = band_lon.ReadAsArray()

        lon = ScatterometryMapper.shift_longitudes(lon)
        self.set_gcps(lon, lat, gdal_dataset)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        ii = pti.get_gcmd_instrument('ascat')
        pp = pti.get_gcmd_platform(metadata['NC_GLOBAL#source'].split(' ')[0])
        provider = pti.get_gcmd_provider(re.split('[^a-zA-Z]',
            metadata['NC_GLOBAL#institution'])[0])

        # TODO: Validate that the found instrument and platform are indeed what
        # we want....

        self.dataset.SetMetadataItem('instrument', json.dumps(ii))
        self.dataset.SetMetadataItem('platform', json.dumps(pp))
        self.dataset.SetMetadataItem('Data Center', json.dumps(provider))
        self.dataset.SetMetadataItem('Entry Title', metadata['NC_GLOBAL#title'])

    def times(self):
        """ Get times from time variable
        """
        ds = Dataset(self.input_filename)

        # Get datetime object of epoch and time_units string
        time_units = self._time_units(ds=ds)

        # Get all times - slight difference from NetCDF-CF mappers times method...
        times = ds.variables[self._timevarname(ds=ds)][:,0]

        # Create numpy array of np.datetime64 times (provide epoch to save time)
        tt = np.array([self._time_count_to_np_datetime64(tn,
            time_units=time_units) for tn in times])

        return tt
