# Name:    nsr.py
# Purpose: Container of NSR class
# Authors:      Anton Korosov
# Created:      01.02.2014
# Copyright:    (c) NERSC 2011 - 2014
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
from __future__ import absolute_import, unicode_literals
import sys
import osr

from nansat.exceptions import NansatProjectionError


class NSR(osr.SpatialReference, object):
    """Nansat Spatial Reference. Overrides constructor of osr.SpatialReference.

    Parameters
    ----------
    srs : 0, PROJ4 or EPSG or WKT or osr.SpatialReference, NSR
        Specifies spatial reference system (SRS)
        PROJ4:
        string with proj4 options [http://trac.osgeo.org/proj/] e.g.:
        '+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs'
        '+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=0 +no_defs'
        EPSG:
        integer with EPSG number, [http://spatialreference.org/],
        e.g. 4326
        WKT:
        string with Well Know Text of SRS. E.g.:
        'GEOGCS["WGS 84",
            DATUM["WGS_1984",
                SPHEROID["WGS 84",6378137,298.257223563,
                    AUTHORITY["EPSG","7030"]],
                TOWGS84[0,0,0,0,0,0,0],
                AUTHORITY["EPSG","6326"]],
            PRIMEM["Greenwich",0,
                AUTHORITY["EPSG","8901"]],
            UNIT["degree",0.0174532925199433,
                AUTHORITY["EPSG","9108"]],
            AUTHORITY["EPSG","4326"]]'

    """

    def __init__(self, srs=0):
        """Create Spatial Reference System from input parameter"""
        if sys.version_info.major == 2:
            str_types = (str, unicode)
        else:
            str_types = str
        # create SRS
        osr.SpatialReference.__init__(self)

        # parse input parameters
        status = 1
        if srs is 0:
            # generate default WGS84 SRS
            status = self.ImportFromWkt(osr.SRS_WKT_WGS84)
        elif isinstance(srs, str_types):
            # parse as proj4 string
            status = self.ImportFromProj4(str(srs))
            if status > 0:
                # parse as WKT string
                status = self.ImportFromWkt(str(srs))
            if status > 0:
                raise NansatProjectionError('Proj4 or WKT (%s) is wrong' % srs)
        # TODO: catch long in python 3
        elif isinstance(srs, int):
            # parse as EPSG code
            status = self.ImportFromEPSG(srs)
            if status > 0:
                raise NansatProjectionError('EPSG %d is wrong' % srs)
        elif type(srs) in [osr.SpatialReference, NSR]:
            # parse from input Spatial Reference
            status = self.ImportFromWkt(srs.ExportToWkt())
            if status > 0:
                raise NansatProjectionError('NSR %s is wrong' % srs)

    @property
    def wkt(self):
        """Well Known Text representation of SRS"""
        return self.ExportToWkt()
