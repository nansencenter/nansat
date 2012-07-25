from dateutil.parser import parse
from struct import unpack

class Envisat():
    '''Methods shared between Envisat mappers'''
    def _set_envisat_time(self, gdalMetadata):
        ''' Get time from metadata, set time to VRT'''
        # set time
        productTime = gdalMetadata["SPH_FIRST_LINE_TIME"]
        self._set_time(parse(productTime))

    def find_offset(self, fileName, gadsDSName, stepDict, dataDict={}):
        '''
        Return a band wchich include parameters (offsets and steps)
        for creating Geolocation Array metadata

        Parameters
        ----------
            fileName: string
                       fileName of the underlying data
            gadsDSName : string
            stepDict : dictionary
                key is a name of ADS.
                value is offset form the beginning of gadsDSName
            dataDict : dictionary
                key is a name in "DS_OFFSET".
                value is offset from the beginning of "DS_OFFSET".

        Returns
        -------
            offsetDict : dictionary
                keys are keys in stepDict and dataDict.
                values are offset from the beginning of the file
        '''
        # open file, find offset
        f = file(fileName, 'rt')
        headerLines = f.readlines(150)
        offsetDict = {}

        # create a dictionary which has offsets
        if gadsDSName in headerLines:
            gridOffset = headerLines.index(gadsDSName)
            for iKey in stepDict:
                offsetDict[iKey]  = int(headerLines[gridOffset + stepDict[iKey]].replace(iKey+"=", '').replace('<bytes>', ''))
            if dataDict != {}:
                for jkey in dataDict:
                    offsetDict[jkey] = int(offsetDict["DS_OFFSET"]) + dataDict[jkey]["offset"]
        f.close()
        return offsetDict


    def read_allList(self, fileName, offsetDict, keyName, calcsize, numOfReadData):
        '''
        Read data based on indices of data

        Parameters
        ----------
            fileName : string
            offsetDict : dictionary
                key is title(name) of ADS, values is index of the key
            keyName : string
                key in ADSR_list
            calcsize: data type format
            numOfReadData: int
                number of reading data
                step = {'LINES_PER_TIE_PT':-4, 'SAMPLES_PER_TIE_PT':-3}
        Returns
        -------
            allGADSValues : list
                includes values which are read from the file.
                the number of elements is numOdReadData.
        '''
        # fseek to gads, read all into a list
        f = file(fileName, 'rb')
        f.seek(offsetDict[keyName], 0)
        allGADSValues = []
        for i in range(numOfReadData):
            fbString = f.read(4)
            fbVal = unpack(calcsize, fbString)[0]
            allGADSValues.append(fbVal)
        f.close()
        return allGADSValues

    def read_scaling_gads(self, fileName, indeces):
        ''' Read Scaling Factor GADS to get scalings of MERIS L1/L2

        Parameters
        ----------
            fileName : string
            indeces : list

        Returns
        -------
            list
        '''
        maxGADS = max(indeces) + 1
        # open file, find offset
        gadsDSName = 'DS_NAME="Scaling Factor GADS         "\n'
        step = {"DS_OFFSET" : 3}
        offsetDict = self.find_offset(fileName, gadsDSName, step)
        allGADSValues = self.read_allList(fileName, offsetDict, "DS_OFFSET", '>f', maxGADS)

        #get only values required for the mapper
        return [allGADSValues[i] for i in indeces]

    def get_ADSRlist(self, fileType):
        '''
        Parameters
        ----------
            fileType : string ("MER_" or "ASA_")

        Returns
        -------
            gadsDSName : string
            ADSR_list : list
                includes "key name for geolocation", "offset", "datatype" and "unit"

        See also
        ---------
            Meris : http://earth.eo.esa.int/pcs/envisat/meris/documentation/meris_3rd_reproc/Vol11_Meris_6a.pdf (--> p52)
            ASAR :  http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
        '''

        if fileType == "MER_":
            gadsDSName = 'DS_NAME="Tie points ADS              "\n'
            ADSR_list = {
             "Dim" : 71,
             "temp"                      : {"offset" : 0          , "datatype" : "int32" , "unit" : "(10)^-6 deg"},
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
             "num_lines"                    : {"offset" : 13                 , "datatype" : "int"    , "unit" : ""},
             "first_samp_numbers"           : {"offset" : 25+11*4*0          , "datatype" : "float32", "unit" : ""},
             "first_slant_range_times"      : {"offset" : 25+11*4*1          , "datatype" : "float32", "unit" : "ns"},
             "first_line_incidenceAngle"    : {"offset" : 25+11*4*2          , "datatype" : "float32", "unit" : "deg"},
             "first_line_lats"              : {"offset" : 25+11*4*3          , "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
             "first_line_longs"             : {"offset" : 25+11*4*4          , "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
             "last_line_samp_numbers"       : {"offset" : 25+11*4*5+34+11*4*0, "datatype" : "int32"  , "unit" : ""},
             "last_line_slant_range_times"  : {"offset" : 25+11*4*5+34+11*4*1, "datatype" : "float32", "unit" : "ns"},
             "last_line_incidenceAngle"     : {"offset" : 25+11*4*5+34+11*4*2, "datatype" : "float32", "unit" : "deg"},
             "last_line_lats"               : {"offset" : 25+11*4*5+34+11*4*3, "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
             "last_line_longs"              : {"offset" : 25+11*4*5+34+11*4*4, "datatype" : "int32"  , "unit" : "(10)^-6 deg"},
            }

        return gadsDSName, ADSR_list

    def get_parameters(self, fileName, fileType, data_key):
        '''
        Preparation for _creats_bands_().
        Return band size and dictionary which includes
        'source','sourceBand','wkv','parameters' and 'SourceType' as the keys.

        Parameters
        ----------
            fileName: string
                       fileName of the underlying data
            fileType : string
                       "MER_" or "ASA_"
            data_key : list
                       element should be one/some of keys in ADSR_list

        Returns
        -------
            dim, offsetDict["NUM_DSR"] : int
                        XSize and YSize of the band
            metaDict : dictionary
                        parameters for _creats_bands_()

        '''
        gadsDSName, ADSR_list = self.get_ADSRlist(fileType)

        dim = ADSR_list["Dim"]
        dataDict = {}
        for key in data_key:
            if key in ADSR_list:
                dataDict[key] = ADSR_list[key]
        step = {'NUM_DSR':5, 'DSR_SIZE':6, 'DS_OFFSET':3}

        offsetDict = self.find_offset(fileName, gadsDSName, step, dataDict)

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
            metaDict.append({'source': fileName, 'sourceBand': 0, 'wkv': ikey,
                             'parameters': parameters, "SourceType":"RawRasterBand" })

        return dim,  offsetDict["NUM_DSR"], metaDict

    def get_GeoArrayParameters(self, fileName, fileType, data_key=[]):
        '''
        Return a band wchich include parameters (offsets and steps)
        for creating Geolocation Array metadata

        Parameters
        ----------
            fileName: string
                       fileName of the underlying data
            fileType : string
                       "MER_" or "ASA_"
            data_key : list
                       elements should be latitude and longitude key names in ADSR_list

        Returns
        -------
            geolocParameter : list
                [pixelOffset, lineOffset, pixelStep, lineStep]

        '''
        if fileType == "MER_":
            gadsDSName = 'DS_NAME="Quality ADS                 "\n'
            step = {'LINES_PER_TIE_PT':-4, 'SAMPLES_PER_TIE_PT':-3}
            offsetDict = self.find_offset(fileName, gadsDSName, step)

            geolocParameter = [0, 0, offsetDict["SAMPLES_PER_TIE_PT"],
                               offsetDict["LINES_PER_TIE_PT"]]

        elif fileType == "ASA_":
            gadsDSName, ADSR_list = self.get_ADSRlist(fileType)

            dataDict = {}
            for ikey in data_key:
                if ikey in ADSR_list:
                    dataDict[ikey] = ADSR_list[ikey]
            step = {'DS_OFFSET':3}
            offsetDict = self.find_offset(fileName, gadsDSName, step, dataDict)
            # fseek to gads, read all into a list
            allGADSValues = self.read_allList(fileName, offsetDict, "num_lines", '>i', 14)

            geolocParameter = [allGADSValues[3]-1, allGADSValues[0]-1,
                               allGADSValues[4]-allGADSValues[3],
                               allGADSValues[1]]

        return geolocParameter

#m = MERIS();
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_1PNPDK20110817_110451_000004053105_00339_49491_7010.N1', range(7, 22))
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(7, 20) + [20, 21, 22, 20])
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(33, 46) + [46, 47, 48, 46])
