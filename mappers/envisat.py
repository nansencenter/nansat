from dateutil.parser import parse
from struct import unpack
from vrt import VRT

class Envisat():
    '''Methods shared between Envisat mappers'''
    def _set_envisat_time(self, gdalMetadata):
        ''' Get time from metadata, set time to VRT'''
        # set time
        productTime = gdalMetadata["SPH_FIRST_LINE_TIME"]
        self._set_time(parse(productTime))

    def read_scaling_gads(self, fileName, indeces):
        ''' Read Scaling Factor GADS to get scalings of MERIS L1/L2'''

        maxGADS = max(indeces) + 1

        # open file, find offset
        f = file(fileName, 'rt')
        headerLines = f.readlines(100)
        gadsDSNameString = 'DS_NAME="Scaling Factor GADS         "\n'
        if gadsDSNameString in headerLines:
            i1 = headerLines.index(gadsDSNameString)
            iGadsOffset = i1 + 3
            scalingGADSOffset = int(headerLines[iGadsOffset].\
                               replace('DS_OFFSET=', '').replace('<bytes>', ''))
        f.close()

        # fseek to gads, read all into a list
        f = file(fileName, 'rb')
        f.seek(scalingGADSOffset, 0)
        allGADSValues = []
        for i in range(maxGADS):
            fbString = f.read(4)
            fbVal = unpack('>f', fbString)[0]
            allGADSValues.append(fbVal)
        f.close()

        #get only values required for the mapper
        return [allGADSValues[i] for i in indeces]

    def create_geoDataset(self, fileName, gadsDSName, parameters, dataType):
        ''' Create a dataset with a geolocation band.

        Create a dataset with VRTRawRasterBand for geolocation.
        The size of the band is small.
        the column is 11 and the row depends on the data.

        Parameters
        ----------
            fileName: string.
                      fileName of the underlying data.
            gadsDSName : string.
                         dataset name.
            parameters : dictionary.
                         key is a name of the dataset.
                         values is an integer that shows the offset.
            dataType : string. (python dataType)
                        "uint16", "int16", "uint32","int32",
                        "float32", "float64" or "complex64"

        Returns
        -------
            geoDataset : dataset with a small geolocation band

        '''
        # open file, find offset
        f = file(fileName, 'rt')
        headerLines = f.readlines(150)
        offsetDict = {}
        bandNames = []

        # create a dictionary which has offsets
        if gadsDSName in headerLines:
            gridOffset = headerLines.index(gadsDSName)
            for iDic in parameters:
                keys = iDic.keys()
                subDict = None
                if "substream" in keys:
                    subDict = iDic["substream"]
                    keys.remove("substream")
                for ikey in keys:
                    offset = gridOffset + iDic[ikey]
                    offsetDict[ikey] = int(headerLines[offset].replace(ikey+'=', '').replace('<bytes>', ''))
                    if subDict is not None:
                        for j, jkey in enumerate (subDict.keys()):
                            bandNames.append(jkey)
                            offsetDict[bandNames[j]] = offsetDict[ikey] + int(subDict[bandNames[j]])

        # convert python dataType to gdal dataType index
        dataType = { "uint16": 2, "int16": 3 , "uint32": 4,
                     "int32": 5 , "float32": 6, "float64": 7,
                     "complex64": 11}.get(dataType, 6)
        # get size of the dataType in byte
        pixOffset = {2:2, 3:2, 4:4, 5:4, 6:4, 7:8, 11:8}.get(dataType, 4)

        # Create a small empty VRT for incident angle
        geoDataset = VRT(srcRasterXSize=11, srcRasterYSize=offsetDict["NUM_DSR"])

        metaDict = []
        for iBandName in bandNames:
            # Make a dictionary that is parameters for creating a band
            parameters = { "ImageOffset" : offsetDict[iBandName],
                           "PixelOffset" : pixOffset,
                           "LineOffset" : offsetDict["DSR_SIZE"],
                           "ByteOrder" : "MSB", "dataType": dataType,
                           "band_name": iBandName}
            metaDict.append({'source': fileName, 'sourceBand': 0, 'wkv': '', 'parameters': parameters })
        # Add VRTRawRasterband
        geoDataset._create_bands(metaDict)
        return geoDataset

#m = MERIS();
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_1PNPDK20110817_110451_000004053105_00339_49491_7010.N1', range(7, 22))
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(7, 20) + [20, 21, 22, 20])
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(33, 46) + [46, 47, 48, 46])
