# Name:        mapper_landsat
# Purpose:     Mapping for LANDSAT.tar.gz
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import tarfile
import warnings
import ast

from nansat.tools import WrongMapperError
from nansat.tools import gdal, ogr, np
from nansat.vrt import VRT
from nansat.node import Node


class Mapper(VRT):
    ''' Mapper for LANDSAT3,4,5,6,7,8.tar.gz files'''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 resolution='standard', **kwargs):
        ''' Create LANDSAT VRT '''
        # try to open .tar or .tar.gz or .tgz file with tar
        try:
            tarFile = tarfile.open(fileName)
        except:
            raise WrongMapperError

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

        if not metaDict:
            raise WrongMapperError

        # Create a structured array with xSize, ySize, band numbers
        # and number of bands which has same size
        dtype=[('xSize', '>i4'), ('ySize', '>i4'),
               ('bands', '|S30'), ('NumOfBands', '>i4')]
        bandTypes = np.zeros((1,), dtype=dtype)

        for iBand in range(len(metaDict)):
            tmpName = metaDict[iBand]['src']['SourceFilename']
            gdalDatasetTmp = gdal.Open(tmpName)
            if iBand == 0:
                bandTypes = np.array([(gdalDatasetTmp.RasterXSize,
                                       gdalDatasetTmp.RasterYSize,
                                       str(iBand),
                                       1)],
                                     dtype=dtype )
            else:
                for i, iBandType in enumerate(bandTypes):
                    if (iBandType['xSize'] == gdalDatasetTmp.RasterXSize and
                        iBandType['ySize'] == gdalDatasetTmp.RasterYSize):
                        iBandType['bands'] += ',' +str(iBand)
                        iBandType['NumOfBands'] += 1
                        break
                    elif i == len(bandTypes) - 1:
                        bandTypes = np.append(
                                        bandTypes,
                                        np.array([(gdalDatasetTmp.RasterXSize,
                                                   gdalDatasetTmp.RasterYSize,
                                                   str(iBand),
                                                   1)],
                                                 dtype = dtype))

        sizeDiffBands = []
        if resolution == 'high':
            bandTypes = np.sort(bandTypes, order=['xSize', 'ySize'])
        else:
            bandTypes = np.sort(bandTypes, order=['NumOfBands'])
        # get created band size
        xSize = bandTypes[-1]['xSize']
        ySize = bandTypes[-1]['ySize']
        # get band numbers which should be removed
        for iBandType in np.delete(bandTypes, -1, 0):
            if type(ast.literal_eval(iBandType['bands'])) == list:
                sizeDiffBands.extend(ast.literal_eval(iBandType['bands']))
            else:
                sizeDiffBands.append(ast.literal_eval(iBandType['bands']))

        gdalDatasetTmp0 = gdal.Open(metaDict[0]['src']['SourceFilename'])

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDatasetTmp0)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        if len(sizeDiffBands) != 0:
            # return only high resolution bands.
            if resolution == 'high':
                xSizeDef = self.dataset.RasterXSize
                ySizeDef = self.dataset.RasterYSize

                vrtXML = self.read_xml()
                node0 = Node.create(vrtXML)
                node0.replaceAttribute('rasterXSize', str(xSize))
                node0.replaceAttribute('rasterYSize', str(ySize))
                # write contents
                self.write_xml(node0.rawxml())

                # modift GeoTarnsform for the new rasterSize
                geoTransform = list(self.dataset.GetGeoTransform())
                geoTransform[1] = geoTransform[1] / (float(xSize) / xSizeDef)
                geoTransform[5] = geoTransform[5] / (float(ySize) / ySizeDef)
                self.dataset.SetGeoTransform(tuple(geoTransform))
                # delete low resolution bands
                self.delete_bands(list(np.array(sizeDiffBands[0]) + 1))

            # 8th band of LANDSAT8 is a double size band.
            # Reduce the size to same as the 1st band.
            else:
                vrtXML = self.read_xml()
                node0 = Node.create(vrtXML)
                for iBand in sizeDiffBands:
                    iBandNode = node0.nodeList('VRTRasterBand')[iBand]
                    iNodeDstRect = iBandNode.node('DstRect')
                    iNodeDstRect.replaceAttribute('xSize', str(xSize))
                    iNodeDstRect.replaceAttribute('ySize', str(ySize))
                # write contents
                self.write_xml(node0.rawxml())

