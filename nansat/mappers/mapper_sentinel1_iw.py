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

from nansat import Nansat, Domain
from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError


class Mapper(VRT):

    def __init__(self, filename, *args, **kwargs):
        Mapper.check_input(filename)
        n = Nansat(filename)
        dom = Mapper.generate_domain(n)
        n.vrt.dataset.SetGCPs(dom.vrt.dataset.GetGCPs(), dom.vrt.dataset.GetProjection())

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

    @staticmethod
    def generate_domain(n):
        lon = n['rvlLon']
        lat = n['rvlLat']
        dom = Domain(lon=lon, lat=lat)
        return dom
