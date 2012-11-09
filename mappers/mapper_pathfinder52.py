#-------------------------------------------------------------------------------
# Name:        mapper_pathfinder52
# Purpose:     Mapping for MODIS-L1 data
#
# Author:      antonk
#
# Created:     13.12.2011
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from datetime import datetime, timedelta

import numpy as np

from vrt import VRT, gdal, parse

from nansat_tools import latlongSRS

class Mapper(VRT):
    ''' Mapper PATHFINDER (local files)
    
    TODO:
    * remote files
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create VRT '''
        
        assert 'AVHRR_Pathfinder-PFV5.2' in fileName, 'pathfinder52 BAD MAPPER'

        subDatasets = gdalDataset.GetSubDatasets()
        
        metaDict = []

        sstName = ''
        
        for subDataset in subDatasets:
            subDatasetName = subDataset[1].split(' ')[1]
            if '//' in subDatasetName:
                h5Style = True
            else:
                h5Style = False
            
            if h5Style:
                subDatasetName = subDatasetName.replace('//', '')

            if subDatasetName == 'sea_surface_temperature':
                sstName = subDataset[0]
                
            subGDALDataset = gdal.Open(subDataset[0])
            subGDALMetadata = subGDALDataset.GetRasterBand(1).GetMetadata()
            if h5Style:
                metaPrefix = subDatasetName + '_'
            else:
                metaPrefix = ''

            subWKV = subGDALMetadata.get(metaPrefix + 'standard_name', '')
            subScaleRatio = subGDALMetadata.get(metaPrefix + 'scale_factor', '1')
            subScaleOffset = subGDALMetadata.get(metaPrefix + 'add_offset', '0')
            metaEntry = {'src': {'SourceFilename': subDataset[0],
                                     'sourceBand':  1,
                                     'ScaleRatio': subScaleRatio,
                                     'ScaleOffset': subScaleOffset},
                         'dst': {'wkv': subWKV}}

            # append band metadata to metaDict
            metaDict.append(metaEntry)
            
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, subGDALDataset)
        
        # add mask
        if sstName != '':
            sstDataset = gdal.Open(sstName)
            sstArray = sstDataset.ReadAsArray()
            mask = np.zeros(sstArray.shape, 'uint8')
            mask[:] = 128 # all valid
            mask[sstArray < 0] = 1 # cloud pixels
            self.maskVRT = VRT(array=mask)
            metaDict.append(
                {'src': {'SourceFilename': self.maskVRT.fileName, 'SourceBand':  1},
                 'dst': {'name': 'mask'}})
            
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # append fixed projection and geotransform
        self.dataset.SetProjection(latlongSRS.ExportToWkt())
        self.dataset.SetGeoTransform((-180, 0.0417, 0, 90, 0, -0.0417))
        
        # set TIMEstart_time
        if h5Style:
            startTimeKey = 'start_time'
        else:
            startTimeKey = 'NC_GLOBAL#start_time'
        self._set_time(parse(subGDALDataset.GetMetadataItem(startTimeKey)))
