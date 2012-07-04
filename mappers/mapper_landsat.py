#-------------------------------------------------------------------------------
# Name:        mapper_landsat
# Purpose:     Mapping for LANDSAT.tar.gz
#
# Author:      antonk
#
# Created:     04.07.2012
# Copyright:   (c) NERSC 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import VRT
import tarfile

import gdal

class Mapper(VRT):
    ''' Mapper for LANDSAT3,4,5,6,7.tar.gz files'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        ''' Create LANDSAT VRT '''
        # try to open .tar or .tar.gz or .tgz file with tar
        tarFile = tarfile.open(fileName)
        
        tarNames = tarFile.getnames()
        print tarNames
        metaDict = []
        for tarName in tarNames:
            if tarName[0] == 'L' and (tarName[-4:] == '.TIF' or
                                      tarName[-4:] == '.tif'):
                print tarName
                metaDict.append({
                    'source': '/vsitar/%s/%s' % (fileName, tarName), 'sourceBand':  1, 'wkv': 'surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water'
                    })
        print metaDict
        tmpName = metaDict[0]['source']
        print tmpName
        gdalDatasetTmp = gdal.Open(tmpName)
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDatasetTmp, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
