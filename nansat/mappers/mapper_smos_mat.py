# Name:        mapper_smos_mat
# Purpose:     Mapping for SMOS Matlab NERSC-internal format
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os

import numpy as np
from scipy.io import loadmat

from nansat.vrt import VRT
from nansat.tools import gdal, ogr, WrongMapperError


class Mapper(VRT):
    ''' MApper for Matlab files with SMOS data '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create SMOS VRT '''
        # check extension
        fName = os.path.split(fileName)[1]
        fExt = os.path.splitext(fileName)[1]
        if fExt == '.MAT' or fExt == '.mat' and 'OSUDP2' in fName:
            # load file
            matFile = loadmat(fileName)
        else:
            raise WrongMapperError

        # get geolocation
        geolocArray = matFile['geolocation'][0]
        srcProj4 = ('+proj=stere +lon_0=%f +lat_0=%f +datum=WGS84 +ellps=WGS84 +units=km +no_defs'
                    % (geolocArray[0], geolocArray[1]))
        srcProjection = osr.SpatialReference()
        srcProjection.ImportFromProj4(srcProj4)
        srcProjection = srcProjection.ExportToWkt()
        srcGeotransform = (geolocArray[2], geolocArray[4], 0,
                           geolocArray[3], 0, geolocArray[5])
        lon = matFile['longitude']
        #lat = matFile['latitude']
        srcRasterYSize, srcRasterXSize = lon.shape
        # create VRT from lat/lon
        # VRT.__init__(self, lon=lon, lat=lat)
        VRT.__init__(self,
                     srcGeoTransform=srcGeotransform,
                     srcProjection=srcProjection,
                     srcRasterXSize=srcRasterXSize,
                     srcRasterYSize=srcRasterYSize)

        # add the following variables
        varNames = ['SSS1', 'SSS2', 'SSS3', 'SST',
                    'Sigma_SSS1', 'Sigma_SSS2', 'Sigma_SSS3',
                    'Control_Flags_1', 'Control_Flags_2',
                    'Control_Flags_3', 'Control_Flags_4',
                    'Science_Flags_1', 'Science_Flags_2',
                    'Science_Flags_3', 'Science_Flags_4']
        metaDict = []
        for varName in varNames:
            var = matFile[varName]
            self.bandVRTs[varName] = VRT(array=var)
            metaDict.append({'src': {'SourceFilename':
                                     self.bandVRTs[varName].fileName,
                                     'sourceBand': 1},
                            'dst': {'name': varName}})

        # create mask
        cloudBits = [2, 3, 4, 5, 6]
        maxSigma = 3.0
        mask = np.zeros(lon.shape, 'uint16')
        mask[:] = 128
        mask[np.isnan(matFile['SSS1'])] = 0
        mask[matFile['Sigma_SSS1'] > maxSigma] = 1
        mask[matFile['Sigma_SSS2'] > maxSigma] = 1
        mask[matFile['Sigma_SSS3'] > maxSigma] = 1
        for cloudBit in cloudBits:
            for cfi in range(1, 5):
                bitMap = np.bitwise_and(matFile['Control_Flags_%d' % cfi],
                                        np.power(2, cloudBit))
                mask[bitMap > 0] = 1

        self.bandVRTs['mask'] = VRT(array=mask)
        metaDict.append({'src': {'SourceFilename':
                                 self.bandVRTs['mask'].fileName,
                                 'sourceBand': 1},
                        'dst': {'name': 'mask'}})

        self.logger.debug('metaDict: %s' % metaDict)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
