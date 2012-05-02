#-------------------------------------------------------------------------------
# Name:        nansat_mapper_merisL1
# Purpose:     Mapping for MOD44W watermask data
#
# Author:      antonk
#
# Created:     02.03.2012
# Copyright:   (c) NERSC 2012
# Licence:     GPL
#-------------------------------------------------------------------------------
from vrt import *
import os.path

class Mapper(VRT):
    ''' VRT with mapping of WKV for MOD44W produc (MODIS watermask at 250 m)'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        ''' Create VRT '''

        fileBaseName = os.path.basename(fileName)
        if not fileBaseName == 'MOD44W.vrt':
            raise AttributeError("MOD44W BAD MAPPER");
        
        metaDict = [{'source': fileName, 'sourceBand':  1, 'wkv': 'land_binary_mask', 'parameters': {'band_name': 'land_mask'}}];

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
