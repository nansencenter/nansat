# ------------------------------------------------------------------------------
# Name:     mapper_sentinel1_iw.py
# Purpose:
#
# Author:       Artem Moiseev
#
# Created:  05.07.2018
# Copyright:    (c) NERSC
# License: GPL V3
# ------------------------------------------------------------------------------

import os
import re

from nansat import Domain
from nansat.exceptions import WrongMapperError
from nansat.mappers.mapper_generic import Mapper as GenericMapper


class Mapper(GenericMapper):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):
        Mapper.check_input(filename)
        GenericMapper.__init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs)
        dom = self.generate_domain()
        self.dataset.SetGCPs(dom.vrt.dataset.GetGCPs(), dom.vrt.dataset.GetProjection())

    @staticmethod
    def parse_filename(filename):
        _, filename = os.path.split(filename)
        filename_base, file_format = os.path.splitext(filename)
        return filename_base, file_format

    @staticmethod
    def check_input(filename):
        filename_base, file_format = Mapper.parse_filename(filename.strip())

        if file_format != '.nc':
            raise WrongMapperError('Only netCDF format is supported')
        elif not re.fullmatch(r'S1._IW_RVL.*', filename_base):
            raise WrongMapperError()

    def generate_domain(self):
        lon = self.dataset.GetRasterBand(1).ReadAsArray()
        lat = self.dataset.GetRasterBand(2).ReadAsArray()
        dom = Domain(lon=lon, lat=lat)
        return dom
