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
import re
import os.path

from nansat.vrt import VRT
from nansat.tools import WrongMapperError

class Mapper(VRT):
    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        '''
        Mapping for the GTOPO30 (`<https://lta.cr.usgs.gov/GTOPO30>`_) and the GMTED2010
        (`<https://lta.cr.usgs.gov/GMTED2010>`_) global elevation models.

        Parameters: 
        ----------- 
        filename : string 
        
        Either the name of a GTOPO30 DEM file or GMTED2010 tif file, or <path>/<dem>.vrt. The latter
        is an aggregation of the DEM-files available from the given DEM. The GTOPO30 vrt does not
        contain the Antarctic, because this is in polarstereographic projection. 
        
        You can create your own gtopo30.vrt file with gdal:
            > gdalbuildvrt <dem>.vrt [E,W]*.DEM

        Remember to update this mapper by adding allowed filenames to the list of accepted filenames
        (accepted_names) if you create or apply new DEM datasets.

        '''

        bn = os.path.basename(filename)
        accepted_names = [
                'gmted2010_30.vrt',
                'gtopo30.vrt',
                '*.DEM',
                '*_gmted_mea*.tif',
            ]

        correct_mapper = False
        for accepted_name in accepted_names:
            m = re.search(accepted_name, filename)
            if m and m.group(0):
                correct_mapper = True
                break

        #if not bn=='gtopo30.vrt' and not os.path.splitext(bn)[1]=='.DEM':
        if not correct_mapper:
            raise WrongMapperError

        metaDict = [{'src': {'SourceFilename': filename, 'SourceBand':  1},
                     'dst': {'wkv': 'height_above_reference_ellipsoid'}}]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
