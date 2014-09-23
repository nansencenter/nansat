# Name:        mapper_landsat
# Purpose:     Mapping for the highest resolution band of LANDSAT8.tar.gz
# Authors:      Anton Korosov, Asuka Yamakawa
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from vrt import VRT
import tarfile

# import standard and additional libraries
from nansat_tools import *

try:
    from osgeo import gdal
except ImportError:
    import gdal


class Mapper(VRT):
    ''' Mapper for high resolution band of LANDSAT8.tar.gz files'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create LANDSAT VRT '''
        # try to open .tar or .tar.gz or .tgz file with tar
        try:
            tarFile = tarfile.open(fileName)
        except:
            raise WrongMapperError

        tarNames = tarFile.getnames()
        metaDictAll = []
        for tarName in tarNames:
            if ((tarName[0] == 'L' or tarName[0] == 'M') and
               (tarName[-4:] == '.TIF' or tarName[-4:] == '.tif')):
                # crate metadataDict for all mappers
                bandNo = tarName[-6:-4]
                metaDictAll.append({
                    'src': {'SourceFilename': '/vsitar/%s/%s' % (fileName,
                                                                 tarName),
                            'SourceBand':  1},
                    'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                            'suffix': bandNo}})

        # copy metadataDict which has the highest resolution.
        for iFile in range(len(metaDictAll)):
            tmpName = metaDictAll[iFile]['src']['SourceFilename']
            gdalDatasetTmp = gdal.Open(tmpName)
            # set an initial size
            if iFile == 0:
                gdalDatasetTmp0 = gdalDatasetTmp
                xSize0 = gdalDatasetTmp.RasterXSize
                ySize0 = gdalDatasetTmp.RasterYSize
                xSize, ySize = xSize0, ySize0
                metaDict = [metaDictAll[0]]
                ratio = 1.0
            # if size of gdalDatasetTmp is larger than current size, replace
            if (xSize < gdalDatasetTmp.RasterXSize and
                ySize < gdalDatasetTmp.RasterYSize):
                    ratio = float(xSize0) / float(gdalDatasetTmp.RasterXSize)
                    xSize = gdalDatasetTmp.RasterXSize
                    ySize = gdalDatasetTmp.RasterYSize
                    metaDict = [metaDictAll[iFile]]
            # if size of gdalDatasetTmp is same as the current size, append metaDict
            elif (xSize == gdalDatasetTmp.RasterXSize and
                  ySize == gdalDatasetTmp.RasterYSize):
                    metaDict.append(metaDictAll[iFile])

        # modify geoTarnsform for the highest resplution
        geoTransform = list(gdalDatasetTmp.GetGeoTransform())
        geoTransform[1] = float(geoTransform[1]) * ratio
        geoTransform[5] = float(geoTransform[5]) * ratio

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDatasetTmp0)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # 8th band of LANDSAT8 is a double size band.
        # Reduce the size to same as the 1st band.
        vrtXML = self.read_xml()
        node0 = Node.create(vrtXML)
        node0.replaceAttribute('rasterXSize', str(xSize))
        node0.replaceAttribute('rasterYSize', str(ySize))
        self.write_xml(str(node0.rawxml()))

        # set new goeTransform
        if ratio != 1.0:
            self.dataset.SetGeoTransform(tuple(geoTransform))




