# Name:    domain.py
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
from __future__ import absolute_import
from nansat.tools import ProjectionError, osr


class NSR(osr.SpatialReference, object):
    '''Nansat Spatial Reference

    Overrides only constructor of osr.SpatialReference
    '''
    def __init__(self, srs=0):
        '''Create Spatial Reference System from input parameter

        Input
        -----
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

        See Also
        ---------
        [http://www.gdal.org/gdalwarp.html]
        [http://trac.osgeo.org/proj/]
        [http://spatialreference.org/]
        [http://www.gdal.org/ogr/osr_tutorial.html]

        '''
        # create SRS
        osr.SpatialReference.__init__(self)

        # parse input parameters
        status = 1
        if srs is 0:
            # generate default WGS84 SRS
            status = self.ImportFromWkt(osr.SRS_WKT_WGS84)
        elif type(srs) in [str, unicode]:
            # parse as proj4 string
            status = self.ImportFromProj4(str(srs))
            if status > 0:
                # parse as WKT string
                status = self.ImportFromWkt(str(srs))
            if status > 0:
                raise ProjectionError('Proj4 or WKT (%s) is wrong' % srs)
        elif type(srs) in [long, int]:
            # parse as EPSG code
            status = self.ImportFromEPSG(srs)
            if status > 0:
                raise ProjectionError('EPSG %d is wrong' % srs)
        elif type(srs) in [osr.SpatialReference, NSR]:
            # parse from input Spatial Reference
            status = self.ImportFromWkt(srs.ExportToWkt())
            if status > 0:
                raise ProjectionError('NSR %s is wrong' % srs)

        # set WKT
        self.wkt = self.ExportToWkt()
