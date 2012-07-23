from dateutil.parser import parse
from struct import unpack

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
        #geoDataset = VRT(srcRasterXSize=11, srcRasterYSize=offsetDict["NUM_DSR"])

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
    '''
    Reference:
    Meris :    http://earth.eo.esa.int/pcs/envisat/meris/documentation/meris_3rd_reproc/Vol11_Meris_6a.pdf
                --> p52
    ASAR :      http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def get_parameters(self, fileName, fileType, data_key):
        if fileType == "MER_":
            gadsDSName = 'DS_NAME="Tie points ADS              "\n'
            ADSR_list = {
             "Dim" : 71,
             "latitude"                  : {"offset" : 13         , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
             "longitude"                 : {"offset" : 13+284*1   , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
             "DME altitude"              : {"offset" : 13+284*2   , "datatype" : "int32" , "unit" : "m"},
             "DME roughness"             : {"offset" : 13+284*3   , "datatype" : "uint32", "unit" : "m"},
             "DME latitude corrections"  : {"offset" : 13+284*4   , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
             "DME longitude corrections" : {"offset" : 13+284*5   , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
             "sun zenith angles"         : {"offset" : 13+284*6   , "datatype" : "uint32", "unit" : "(10)^-6 deg"},
             "sun azimuth angles"        : {"offset" : 13+284*7   , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
             "viewing zenith angles"     : {"offset" : 13+284*8   , "datatype" : "uint32", "unit" : "(10)^-6 deg"},
             "viewing azimuth angles"    : {"offset" : 13+284*9   , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
             "zonal winds"               : {"offset" : 13+284*10+142*0 , "datatype" : "int16" , "unit" : "m*s-1"},
             "meridional winds"          : {"offset" : 13+284*10+142*1 , "datatype" : "int16" , "unit" : "m*s-1"},
             "mean sea level pressure"   : {"offset" : 13+284*10+142*2 , "datatype" : "uint16", "unit" : "hPa"},
             "total ozone"               : {"offset" : 13+284*10+142*3 , "datatype" : "uint16", "unit" : "DU"},
             "relative humidity"         : {"offset" : 13+284*10+142*4 , "datatype" : "uint16", "unit" : "%"}
            }
        elif fileType == "ASA_":
            gadsDSName = 'DS_NAME="GEOLOCATION GRID ADS        "\n'
            ADSR_list = {
             "Dim" : 11,
             "first_line_samp_numbers"      : {"offset" : 25               , "datatype" : "int32"  , "unit" : ""},
             "first_samp_numbers"           : {"offset" : 25+11*4*0        , "datatype" : "float32", "unit" : "ns"},
             "first_slant_range_times"      : {"offset" : 25+11*4*1        , "datatype" : "float32", "unit" : "deg"},
             "first_line_incidenceAngle"    : {"offset" : 25+11*4*2        , "datatype" : "float32", "unit" : "deg"},
             "first_line_lats"              : {"offset" : 25+11*4*3        , "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
             "first_line_longs"             : {"offset" : 25+11*4*4        , "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
             "last_line_samp_numbers"       : {"offset" : 25+11*5+34+11*4*0, "datatype" : "int32"  , "unit" : ""},
             "last_line_slant_range_times"  : {"offset" : 25+11*5+34+11*4*1, "datatype" : "float32", "unit" : "ns"},
             "last_line_incidenceAngle"     : {"offset" : 25+11*5+34+11*4*2, "datatype" : "float32", "unit" : "deg"},
             "last_line_lats"               : {"offset" : 25+11*5+34+11*4*3, "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
             "last_line_longs"              : {"offset" : 25+11*5+34+11*4*4, "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
            }

        dim = ADSR_list["Dim"]
        dataDict = {}
        for key in data_key:
            if key in ADSR_list:
                dataDict[key] = ADSR_list[key]

        f = file(fileName, 'rt')
        headerLines = f.readlines(150)
        offsetDict = {}
        bandNames = []

        parameters = {"DS_OFFSET": 3, "NUM_DSR" : 5, "DSR_SIZE" : 6}
        # create a dictionary which has offsets
        if gadsDSName in headerLines:
            gridOffset = headerLines.index(gadsDSName)
            offsetDict["NUM_DSR"]  = int(headerLines[gridOffset + parameters["NUM_DSR"]].replace('NUM_DSR=', '').replace('<bytes>', ''))
            offsetDict["DSR_SIZE"] = int(headerLines[gridOffset + parameters["DSR_SIZE"]].replace('DSR_SIZE=', '').replace('<bytes>', ''))
            offsetDict["DS_OFFSET"] = int(headerLines[gridOffset + parameters["DS_OFFSET"]].replace('DS_OFFSET=', '').replace('<bytes>', ''))
            for ikey in dataDict:
                offsetDict[ikey] = int(offsetDict["DS_OFFSET"]) + dataDict[ikey]["offset"]

        metaDict = []
        for ikey in data_key:
            # convert python dataType to gdal dataType index
            dataType = { "uint16": 2, "int16": 3 , "uint32": 4,
                         "int32": 5 , "float32": 6, "float64": 7,
                         "complex64": 11}.get(dataDict[ikey]["datatype"], 6)
            # get size of the dataType in byte
            pixOffset = {2:2, 3:2, 4:4, 5:4, 6:4, 7:8, 11:8}.get(dataType, 4)

            # Make a dictionary that is parameters for creating a band
            parameters = { "ImageOffset" : offsetDict[ikey],
                           "PixelOffset" : pixOffset,
                           "LineOffset" : offsetDict["DSR_SIZE"],
                           "ByteOrder" : "MSB", "dataType": dataType,
                           "band_name": ikey,
                           "unit": dataDict[ikey]["unit"]}
            metaDict.append({'source': fileName, 'sourceBand': 0, 'wkv': ikey, 'parameters': parameters, "SourceType":"RawRasterBand" })

        return dim,  offsetDict["NUM_DSR"], metaDict

#m = MERIS();
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_1PNPDK20110817_110451_000004053105_00339_49491_7010.N1', range(7, 22))
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(7, 20) + [20, 21, 22, 20])
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(33, 46) + [46, 47, 48, 46])
