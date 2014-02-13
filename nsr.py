# Name:    domain.py
# Purpose: Container of Domain class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

from nansat_tools import osr

class NSR(osr.SpatialReference, object):
    '''Nansat Spatial Reference
    
    Overrides only constructor of osr.SpatialReference
    '''
    def __init__(self, srs=0):
        '''Create Spatial Reference System from input parameter
        
        Input
        -----
        srs : 0, None, WKT, proj4, int, osr.SpatialReference, NSR
            0: default WGS 84 SRS is generated
            None: None is created
            string: parsed as proj4 or WKT
            integerer: parsed as EPSG code
            osr.SpatialReference or NSR: same as input
            
        '''
        # create SRS
        osr.SpatialReference.__init__(self)
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
        elif type(srs) in [long, int]:
            # parse as EPSG code
            status = self.ImportFromEPSG(srs)
        elif type(srs) in [osr.SpatialReference, NSR]:
            # parse from input Spatial Reference
            status = self.ImportFromWkt(srs.ExportToWkt())
