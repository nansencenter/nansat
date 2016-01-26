#-------------------------------------------------------------------------------
# Name:         mapper_gtopo30.py
# Purpose:      Mapping for the global 30 arc-second elevation
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	04.06.2015
# Last modified:08.06.2015 10:27
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#-------------------------------------------------------------------------------
import os.path

from nansat.vrt import VRT
from nansat.tools import WrongMapperError

class Mapper(VRT):
    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        '''
        Mapping for the global 30 arc-second elevation (see
        https://lta.cr.usgs.gov/GTOPO30).

        Parameters:
        -----------
        fileName : string
            Either the name of a gtopo30 DEM file, or <path>/gtopo30.vrt. The
            latter is an aggregation of the DEM-files available with gtopo30
            except the Antarctic one, which is in polarstereographic
            projection. You can create your own gtopo30.vrt file with gdal:
            > gdalbuildvrt gtopo30.vrt [E,W]*.DEM
        '''

        bn = os.path.basename(fileName)
        if not bn=='gtopo30.vrt' and not os.path.splitext(bn)[1]=='.DEM':
            raise WrongMapperError

        metaDict = [{'src': {'SourceFilename': fileName, 'SourceBand':  1},
                     'dst': {'wkv': 'height_above_reference_ellipsoid'}}]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
