# Name:        mapper_mod44w
# Purpose:     Mapping for MOD44W watermask data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from nansat.vrt import *
import os.path


class Mapper(VRT):
    ''' VRT with mapping of WKV for MOD44W produc (MODIS watermask at 250 m)'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create VRT '''

        fileBaseName = os.path.basename(fileName)
        if not fileBaseName == 'MOD44W.vrt':
            raise AttributeError("MOD44W BAD MAPPER")

        metaDict = [{'src': {'SourceFilename': fileName, 'SourceBand':  1},
                     'dst': {'wkv': 'land_binary_mask'}}]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
