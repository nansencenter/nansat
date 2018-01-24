# Name:        mapper_pathfinder52
# Purpose:     Mapping for NOAA AVHRR PATHFINDER52 DATA
# Authors:      Anton Korosov, Dmitry Petrenko
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse

import numpy as np

import nansat.vrt as vrt
from nansat.nsr import NSR
from nansat.exceptions import WrongMapperError


class Mapper(vrt.VRT):
    ''' Mapper PATHFINDER (local files)

    TODO:
    * remote files
    '''

    def __init__(self, filename, gdalDataset, gdalMetadata, minQual=4,
                 **kwargs):
        ''' Create VRT '''

        if not 'AVHRR_Pathfinder-PFV5.2' in filename:
            raise WrongMapperError(filename)

        subDatasets = gdalDataset.GetSubDatasets()
        metaDict = []
        sstName = ''

        for subDataset in subDatasets:
            subDatasetName = subDataset[0].split(':')[2]

            if '//' in subDatasetName:
                h5Style = True
            else:
                h5Style = False

            if h5Style:
                subDatasetName = subDatasetName.replace('//', '')

            if subDatasetName == 'quality_level':
                qualName = subDataset[0]

            subGDALDataset = vrt.gdal.Open(subDataset[0])
            subGDALMetadata = subGDALDataset.GetRasterBand(1).GetMetadata()
            if h5Style:
                metaPrefix = subDatasetName + '_'
            else:
                metaPrefix = ''

            subWKV = subGDALMetadata.get(metaPrefix + 'standard_name', '')
            subScaleRatio = subGDALMetadata.get(metaPrefix + 'scale_factor',
                                                '1')
            subScaleOffset = subGDALMetadata.get(metaPrefix + 'add_offset',
                                                 '0')
            metaEntry = {'src': {'SourceFilename': subDataset[0],
                                 'sourceBand': 1,
                                 'ScaleRatio': subScaleRatio,
                                 'ScaleOffset': subScaleOffset},
                         'dst': {'wkv': subWKV}}

            # append band metadata to metaDict
            metaDict.append(metaEntry)

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(subGDALDataset)

        # add mask
        if qualName != '':
            qualDataset = vrt.gdal.Open(qualName)
            qualArray = qualDataset.ReadAsArray()
            qualArray[qualArray < minQual] = 1
            qualArray[qualArray >= minQual] = 128
            self.band_vrts = {'maskVRT': vrt.VRT(array=qualArray.astype('int8'))}
            metaDict.append({'src': {'SourceFilename': (self.
                                                        band_vrts['maskVRT'].
                                                        filename),
                                     'SourceBand': 1,
                                     'SourceType': 'SimpleSource',
                                     'DataType': 1},
                             'dst': {'name': 'mask'}})

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # append fixed projection and geotransform
        self.dataset.SetProjection(NSR().wkt)
        self.dataset.SetGeoTransform((-180, 0.0417, 0, 90, 0, -0.0417))

        # set TIMEstart_time
        if h5Style:
            startTimeKey = 'start_time'
        else:
            startTimeKey = 'NC_GLOBAL#start_time'
        self.dataset.SetMetadataItem('time_coverage_start', subGDALDataset.GetMetadataItem(startTimeKey))
