# Name:        mapper_modisL1
# Purpose:     Mapping for MODIS-L1 data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse
import warnings

from nansat.utils import gdal, ogr
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT


class HDF4Mapper(VRT):

    def find_metadata(self, iMetadata, iKey, default=''):
        """ Find metadata which has similar key

        Parameters
        ----------
            iMetadata : dict
                input metadata, usually gdalMetadata
            iKey : str
                key to search for
            default : str
                default value

        """
        value = default
        for key in iMetadata:
            if iKey in key:
                value = iMetadata[key]
                break

        return value
