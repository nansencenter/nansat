# Name:         mapper_emodnet.py
# Purpose:      Mapper for bathymetry data from EMODNet
#               http://portal.emodnet-bathymetry.eu/
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
from dateutil.parser import parse

import numpy as np

from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.node import Node
from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError


class Mapper(VRT):
    def __init__(self, inputFileName, gdalDataset, gdalMetadata, logLevel=30,
                 **kwargs):
        # check if mapper fits
        if not gdalMetadata:
            raise WrongMapperError
        if not os.path.splitext(inputFileName)[1] == '.mnt':
            raise WrongMapperError
        try:
            mbNorthLatitude = float(gdalMetadata['NC_GLOBAL#mbNorthLatitude'])
            mbSouthLatitude = float(gdalMetadata['NC_GLOBAL#mbSouthLatitude'])
            mbEastLongitude = float(gdalMetadata['NC_GLOBAL#mbEastLongitude'])
            mbWestLongitude = float(gdalMetadata['NC_GLOBAL#mbWestLongitude'])
            mbProj4String = gdalMetadata['NC_GLOBAL#mbProj4String']
            Number_lines = int(gdalMetadata['NC_GLOBAL#Number_lines'])
            Number_columns = int(gdalMetadata['NC_GLOBAL#Number_columns'])
            Element_x_size = float(gdalMetadata['NC_GLOBAL#Element_x_size'])
            Element_y_size = float(gdalMetadata['NC_GLOBAL#Element_y_size'])
        except:
            raise WrongMapperError

        # find subdataset with DEPTH
        subDatasets = gdalDataset.GetSubDatasets()
        dSourceFile = None
        for subDataset in subDatasets:
            if subDataset[0].endswith('.mnt":DEPTH'):
                dSourceFile = subDataset[0]
        if dSourceFile is None:
            raise WrongMapperError
        dSubDataset = gdal.Open(dSourceFile)
        dMetadata = dSubDataset.GetMetadata()

        try:
            scale_factor = dMetadata['DEPTH#scale_factor']
            add_offset = dMetadata['DEPTH#add_offset']
        except:
            raise WrongMapperError

        geoTransform = [mbWestLongitude, Element_x_size, 0,
                        mbNorthLatitude, 0, -Element_y_size]

        # create empty VRT dataset with geolocation only
        self._init_from_dataset_params(Number_columns, Number_lines, geoTransform, NSR(mbProj4String).wkt, metadata=gdalMetadata)

        metaDict = [{'src': {'SourceFilename': dSourceFile,
                             'SourceBand': 1,
                             'ScaleRatio' : scale_factor,
                             'ScaleOffset' : add_offset},
                     'dst': {'wkv': 'depth'}}]

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)
