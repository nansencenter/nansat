# Name:        mapper_landsat
# Purpose:     Mapping for LANDSAT.tar.gz
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import tarfile
import warnings

from nansat.tools import WrongMapperError
from nansat.tools import gdal, ogr
from nansat.vrt import VRT
from nansat.node import Node


class Mapper(VRT):
    ''' Mapper for LANDSAT3,4,5,6,7,8.tar.gz files'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create LANDSAT VRT '''
        # try to open .tar or .tar.gz or .tgz file with tar
        try:
            tarFile = tarfile.open(fileName)
        except:
            raise WrongMapperError(__file__, "Wrong mapper")

        tarNames = tarFile.getnames()
        #print tarNames
        metaDict = []
        for tarName in tarNames:
            if ((tarName[0] == 'L' or tarName[0] == 'M') and
               (tarName[-4:] == '.TIF' or tarName[-4:] == '.tif')):
                #print tarName
                bandNo = tarName[-6:-4]
                metaDict.append({
                    'src': {'SourceFilename': '/vsitar/%s/%s' % (fileName,
                                                                 tarName),
                            'SourceBand':  1},
                    'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                            'suffix': bandNo}})
        #print metaDict
        sizeDiffBands = []
        for iFile in range(len(metaDict)):
            tmpName = metaDict[iFile]['src']['SourceFilename']
            gdalDatasetTmp = gdal.Open(tmpName)
            if iFile == 0:
                gdalDatasetTmp0 = gdalDatasetTmp
                xSize = gdalDatasetTmp.RasterXSize
                ySize = gdalDatasetTmp.RasterYSize
            elif (xSize != gdalDatasetTmp.RasterXSize or
                    ySize != gdalDatasetTmp.RasterYSize):
                sizeDiffBands.append(iFile)

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDatasetTmp0)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # 8th band of LANDSAT8 is a double size band.
        # Reduce the size to same as the 1st band.
        if len(sizeDiffBands) != 0:
            vrtXML = self.read_xml()
            node0 = Node.create(vrtXML)
            for iBand in sizeDiffBands:
                iBandNode = node0.nodeList('VRTRasterBand')[iBand]
                iNodeDstRect = iBandNode.node('DstRect')
                iNodeDstRect.replaceAttribute('xSize', str(xSize))
                iNodeDstRect.replaceAttribute('ySize', str(ySize))

            self.write_xml(node0.rawxml())
