# Name:         mapper_landsat
# Purpose:      Mapping for LANDSAT*.tar.gz
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
import tarfile
import warnings

from nansat.tools import WrongMapperError
from nansat.tools import gdal, np
from nansat.vrt import VRT

class Mapper(VRT):
    ''' Mapper for LANDSAT5,6,7.tar.gz files'''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                       resolution='low', **kwargs):
        ''' Create LANDSAT VRT from tar.gz files'''
        # try to open .tar or .tar.gz or .tgz file with tar
        try:
            tarFile = tarfile.open(fileName)
        except:
            raise WrongMapperError

        # collect names of bands and corresponding sizes
        # into bandsInfo dict and bandSizes list
        tarNames = tarFile.getnames()
        bandInfos = {}
        bandSizes = []
        for tarName in tarNames:
            # check if TIF files inside TAR qualify
            if   (tarName[0] in ['L', 'M'] and
                  os.path.splitext(tarName)[1] in ['.TIF', '.tif']):
                # let last part of file name be suffix
                bandSuffix = os.path.splitext(tarName)[0].split('_')[-1]
                # open TIF file from TAR using VSI
                sourceFilename = '/vsitar/%s/%s' % (fileName, tarName)
                gdalDatasetTmp = gdal.Open(sourceFilename)
                # keep name, GDALDataset and size
                bandInfos[bandSuffix] = [sourceFilename,
                                         gdalDatasetTmp,
                                         gdalDatasetTmp.RasterXSize]
                bandSizes.append(gdalDatasetTmp.RasterXSize)

        # if not TIF files found - not appropriate mapper
        if not bandSizes:
            raise WrongMapperError

        # get appropriate band size based on number of unique size and
        # required resoltuion
        if resolution == 'low':
            bandXSise = min(bandSizes)
        elif resolution in ['high', 'hi']:
            bandXSise = max(bandSizes)
        else:
            raise OptionError('Wrong resolution %s for file %s' % (resolution, fileName))

        # find bands with appropriate size and put to metaDict
        metaDict = []
        for bandInfo in bandInfos:
            if bandInfos[bandInfo][2] == bandXSise:
                metaDict.append({
                    'src': {'SourceFilename': bandInfos[bandInfo][0],
                            'SourceBand':  1},
                    'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                            'suffix': bandInfo}})
                gdalDatasetTmp = bandInfos[bandInfo][1]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDatasetTmp, **kwargs)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
