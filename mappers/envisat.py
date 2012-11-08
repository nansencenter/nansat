from dateutil.parser import parse
from struct import unpack
from vrt import VRT, Geolocation
import gdal

class Envisat():
    '''Methods shared between Envisat mappers'''
    def _set_envisat_time(self, gdalMetadata):
        ''' Get time from metadata, set time to VRT'''
        # set time
        productTime = gdalMetadata["SPH_FIRST_LINE_TIME"]
        self._set_time(parse(productTime))

    def read_text(self, fileName, gadsDSName, textOffset, subValOffset={}):
        ''' Return values or offsets of keys in textOffset and subValOffset .

        Find a location of gadsDSName.
        Adjust the location with textOffset and read the text at the location.
        Convert the text to integer and set it into valueDict.
        If subValOffset is given, adjust the location by adding subValOffset
        and set it into valueDict.
        Return the valueDict which has names of variables and
        values or offset of the variables.

        Parameters
        ----------
            fileName: string
                       fileName of the underlying data
            gadsDSName : string
            textOffset : dictionary
                key is a name of ADS.
                value is offset form the location of gadsDSName
            subValOffset : dictionary
                key is a name in "DS_OFFSET".
                value is offset from the location of "DS_OFFSET".

        Returns
        -------
            valueDict : dictionary
                keys are keys in textOffset and subValOffset.
                values are values of the keys or
                  offsets of the keys from the beginning of the file.
        '''
        # open file and read
        f = file(fileName, 'rt')
        headerLines = f.readlines(150)
        valueDict = {}

        # create a dictionary which has offsets
        if gadsDSName in headerLines:
            # get location of gadsDSName
            gridOffset = headerLines.index(gadsDSName)
            # Adjust the location of the varaibles by adding textOffset.
            # Read a text at the location and convert the text into integer.
            for iKey in textOffset:
                valueDict[iKey]  = int(headerLines[gridOffset + textOffset[iKey]].replace(iKey+"=", '').replace('<bytes>', ''))
            # if subValOffset is given, the offset given by the above step is adjusted
            if subValOffset != {}:
                for jkey in subValOffset:
                    valueDict[jkey] = int(valueDict["DS_OFFSET"]) + subValOffset[jkey]["offset"]
        f.close()
        return valueDict


    def read_all_list(self, fileName, offsetDict, keyName, fmt, numOfReadData):
        '''
        Read binary data and return values corresponding to keyNames

        Parameters
        ----------
            fileName : string
            offsetDict : dictionary
                keys are names of variables, values are offsets from the beginning of the file.
            keyName : string
                name of required variable (key in ADSR_list)
            fmt: data type format
            numOfReadData: int
                a number of reading data

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
            fbVal = unpack(fmt, fbString)[0]
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
        textOffset = {"DS_OFFSET" : 3}
        offsetDict = self.read_text(fileName, gadsDSName, textOffset)
        allGADSValues = self.read_all_list(fileName, offsetDict, "DS_OFFSET", '>f', maxGADS)

        #get only values required for the mapper
        return [allGADSValues[i] for i in indeces]

    def get_ADSRlist(self, fileType):
        '''Return a dictionary which incldes name of geolocation and their offset, dataType and unit.

        Parameters
        ----------
            fileType : string ("MER_" or "ASA_")

        Returns
        -------
            gadsDSName : string
            ADSR_list : list
                includes "key name for geolocation", "offset", "dataType" and "unit"

        See also
        ---------
            Meris : http://earth.eo.esa.int/pcs/envisat/meris/documentation/meris_3rd_reproc/Vol11_Meris_6a.pdf (--> p52)
            ASAR :  http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
        '''

        if fileType == "MER_":
            gadsDSName = 'DS_NAME="Tie points ADS              "\n'
            ADSR_list = {
             "Dim" : 71,
             "latitude"                  : {"offset" : 13         , "dataType" : gdal.GDT_Int32 , "unit" : "(10)^-6 deg"},
             "longitude"                 : {"offset" : 13+284*1   , "dataType" : gdal.GDT_Int32 , "unit" : "(10)^-6 deg"},
             "DME altitude"              : {"offset" : 13+284*2   , "dataType" : gdal.GDT_Int32 , "unit" : "m"},
             "DME roughness"             : {"offset" : 13+284*3   , "dataType" : gdal.GDT_UInt32, "unit" : "m"},
             "DME latitude corrections"  : {"offset" : 13+284*4   , "dataType" : gdal.GDT_Int32 , "unit" : "(10)^-6 deg"},
             "DME longitude corrections" : {"offset" : 13+284*5   , "dataType" : gdal.GDT_Int32 , "unit" : "(10)^-6 deg"},
             "sun zenith angles"         : {"offset" : 13+284*6   , "dataType" : gdal.GDT_UInt32, "unit" : "(10)^-6 deg"},
             "sun azimuth angles"        : {"offset" : 13+284*7   , "dataType" : gdal.GDT_Int32 , "unit" : "(10)^-6 deg"},
             "viewing zenith angles"     : {"offset" : 13+284*8   , "dataType" : gdal.GDT_UInt32, "unit" : "(10)^-6 deg"},
             "viewing azimuth angles"    : {"offset" : 13+284*9   , "dataType" : gdal.GDT_Int32 , "unit" : "(10)^-6 deg"},
             "zonal winds"               : {"offset" : 13+284*10+142*0 , "dataType" : gdal.GDT_Int16 , "unit" : "m*s-1"},
             "meridional winds"          : {"offset" : 13+284*10+142*1 , "dataType" : gdal.GDT_Int16 , "unit" : "m*s-1"},
             "mean sea level pressure"   : {"offset" : 13+284*10+142*2 , "dataType" : gdal.GDT_UInt16, "unit" : "hPa"},
             "total ozone"               : {"offset" : 13+284*10+142*3 , "dataType" : gdal.GDT_UInt16, "unit" : "DU"},
             "relative humidity"         : {"offset" : 13+284*10+142*4 , "dataType" : gdal.GDT_UInt16, "unit" : "%"}
            }
        elif fileType == "ASA_":
            gadsDSName = 'DS_NAME="GEOLOCATION GRID ADS        "\n'
            ADSR_list = {
             "Dim" : 11,
             "num_lines"                    : {"offset" : 13                 , "dataType" : gdal.GDT_Int16  , "unit" : ""},
             "first_samp_numbers"           : {"offset" : 25+11*4*0          , "dataType" : gdal.GDT_Float32, "unit" : ""},
             "first_slant_range_times"      : {"offset" : 25+11*4*1          , "dataType" : gdal.GDT_Float32, "unit" : "ns"},
             "first_line_incidenceAngle"    : {"offset" : 25+11*4*2          , "dataType" : gdal.GDT_Float32, "unit" : "deg"},
             "first_line_lats"              : {"offset" : 25+11*4*3          , "dataType" : gdal.GDT_Int32  , "unit" : "(10)^-6 deg"},
             "first_line_longs"             : {"offset" : 25+11*4*4          , "dataType" : gdal.GDT_Int32  , "unit" : "(10)^-6 deg"},
             "last_line_samp_numbers"       : {"offset" : 25+11*4*5+34+11*4*0, "dataType" : gdal.GDT_Int32  , "unit" : ""},
             "last_line_slant_range_times"  : {"offset" : 25+11*4*5+34+11*4*1, "dataType" : gdal.GDT_Float32, "unit" : "ns"},
             "last_line_incidenceAngle"     : {"offset" : 25+11*4*5+34+11*4*2, "dataType" : gdal.GDT_Float32, "unit" : "deg"},
             "last_line_lats"               : {"offset" : 25+11*4*5+34+11*4*3, "dataType" : gdal.GDT_Int32  , "unit" : "(10)^-6 deg"},
             "last_line_longs"              : {"offset" : 25+11*4*5+34+11*4*4, "dataType" : gdal.GDT_Int32  , "unit" : "(10)^-6 deg"},
            }

        return gadsDSName, ADSR_list

    def get_rawband_parameters(self, fileName, fileType, data_key):
        '''Return all parameters to create a VRTRawRasterBand.

        Get ADSR_list of required variables given by data_key.
        Make parameters to create RawRasterBands.
        Create dictionary(metaDict) which has format for_create_bands().
        Return band size to create VRT and dictionary to add bands.

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
            metaDict : list
                        elements are dictionaries that are parameters for _creats_bands_()

        '''
        # Get gadsDSName and ADSR_list corresoinding to the given fileType
        gadsDSName, ADSR_list = self.get_ADSRlist(fileType)

        # get dimension of the data
        dim = ADSR_list["Dim"]
        adsrDict = {}

        # pick up required dictionaries from ADSR_List
        for key in data_key:
            if key in ADSR_list:
                adsrDict[key] = ADSR_list[key]

        # create a dictionary which represent offset from location of gadsDSName
        textOffset = {'NUM_DSR':5, 'DSR_SIZE':6, 'DS_OFFSET':3}

        # Get offsets of required variables from the beginning of the file
        offsetDict = self.read_text(fileName, gadsDSName, textOffset, adsrDict)

        metaDict = []
        # prepare parameters to create bands
        for ikey in data_key:
            # get size of the dataType in bytes
            pixOffset = {gdal.GDT_Byte:     1,
                         gdal.GDT_Int16:    2,
                         gdal.GDT_UInt16:   2,
                         gdal.GDT_Int32:    4,
                         gdal.GDT_UInt32:   4,
                         gdal.GDT_Float32:  4,
                         gdal.GDT_Float64:  8,
                         gdal.GDT_CInt16:   2,
                         gdal.GDT_CInt32:   4,
                         gdal.GDT_CFloat32: 4,
                         gdal.GDT_CFloat64: 8}[adsrDict[ikey]["dataType"]]
                         
            # Append dictionary with parameters for this band
            metaDict.append({'src': {'SourceFilename': fileName,
                                'SourceBand': 0,
                                "SourceType": "RawRasterBand",
                                "ImageOffset" : offsetDict[ikey],
                                "PixelOffset" : pixOffset,
                                "LineOffset" : offsetDict["DSR_SIZE"],
                                "ByteOrder" : "MSB"},
                             'dst': {"dataType": adsrDict[ikey]["dataType"],
                                "name": ikey,
                                'wkv': ikey,
                                "unit": adsrDict[ikey]["unit"]}})
        return dim, offsetDict["NUM_DSR"], metaDict

    def get_geoarray_parameters(self, fileName, fileType, stepSize, data_key=[]):
        ''' Return parameters for Geolocation Domain Metadata

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
        if stepSize != 0:
            geolocParameter = [0, 0, stepSize, stepSize]

        # if MERIS
        elif fileType == "MER_":
            gadsDSName = 'DS_NAME="Quality ADS                 "\n'
            textOffset = {'LINES_PER_TIE_PT':-4, 'SAMPLES_PER_TIE_PT':-3}
            # get values of 'LINES_PER_TIE_PT' and 'SAMPLES_PER_TIE_PT'
            valueDict = self.read_text(fileName, gadsDSName, textOffset)
            # create parameters for Geolocation Domain Metadata. (offset = 0)
            geolocParameter = [0, 0, valueDict["SAMPLES_PER_TIE_PT"],
                               valueDict["LINES_PER_TIE_PT"]]

        # if ASAR
        elif fileType == "ASA_":
            gadsDSName, ADSR_list = self.get_ADSRlist(fileType)
            dataDict = {data_key[0] : ADSR_list[data_key[0]]}
            textOffset = {'DS_OFFSET':3}
            # get data format
            fmt = {"int":">i", gdal.GDT_Int32:">i", gdal.GDT_UInt32 : ">I",
                   gdal.GDT_Float32:">f"}.get(ADSR_list[data_key[0]]["dataType"], ">f")
            # read text data and get the offset of required key
            offsetDict = self.read_text(fileName, gadsDSName, textOffset, dataDict)
            # Read binary data from offset
            allGADSValues = self.read_all_list(fileName, offsetDict, data_key[0], fmt, 14)
            # create parameters for Geolocation Domain Metadata
            geolocParameter = [allGADSValues[3]-1, allGADSValues[0]-1,
                               allGADSValues[4]-allGADSValues[3],
                               allGADSValues[1]]

        return geolocParameter

    def create_VRT_with_rawbands(self, fileName, fileType, data_key):
        ''' Create VRT with some small bands

        Get parameters for createing VRT.
        Create a empty VRT and add bands.

        Parameters
        ----------
            fileName: string
                       fileName of the underlying data
            fileType : string
                       "MER_" or "ASA_"
            data_key : list
                       elements should be one/some of keys in ADSR_list

        Returns:
        --------
            VRT : includes some VRTRawRasterBands

        '''
        # Get parameters for VRTRawRasterBand
        XSize, YSize, metaDict = self.get_rawband_parameters(fileName, fileType, data_key)
        print 'metaDict', metaDict
        # Create dataset with small band
        vrt = VRT(srcRasterXSize=XSize, srcRasterYSize=YSize)
        # Add VRTRawRasterBand
        vrt._create_bands(metaDict)
        return vrt

    def add_geoarray_dataset(self, fileName, fileType, XSize, YSize, latlonName, srs, parameters=[], stepSize = 2):
        ''' Add geolocation domain metadata to the dataset

        Create VRT which has lat and lon VRTRawRasterBands.
        Get parameter for geolocation domain metadata (steps and offsets).
        Create latitude and longitude VRT which has original units. (/1000000.0)
        Add the bands in the latitude and longitude VRT
        as X_DATASET and Y_DATASET in geolocation domain metadata.

        Return a band wchich include parameters (offsets and steps)
        for creating Geolocation Array metadata

        Parameters
        ----------
            fileName: string
                       fileName of the underlying data
            fileType : string
                       "MER_" or "ASA_"
            latlonName : dictionary
                        keys are "latitude" and "longitude" and
                        the values are keys which correspond to "longitude" and "latitude" in ADSR_list
            srs :  string. SRS
            parameters : list, optional
                       elements keys in ADSR_list
            stepSize : int > 0
                        step size for geolocation array. default is 2.
                        If stepSize=0, geolocation array is not enlarged.
                        If stepSize > 0, geolocation array is enlarged to (XSize/stepSize, YSize/stepSize).
                        If stepSize=1, geolocation size is same as the band. it takes time to calculate the reprojection.

        Modifies:
        ---------
            Add Geolocation Array metadata

        '''
        # Create dataset with VRTRawRasterbands
        geoVRT = self.create_VRT_with_rawbands(fileName, fileType, [latlonName["latitude"], latlonName["longitude"]])

        # Get geolocParameter which is required for adding geolocation array metadata
        geolocParameter = self.get_geoarray_parameters(fileName, fileType, stepSize, parameters)

        # Create lat and lon VRT with original units
        lonVRT = VRT(array=geoVRT.dataset.GetRasterBand(2).ReadAsArray()/1000000.0)
        latVRT = VRT(array=geoVRT.dataset.GetRasterBand(1).ReadAsArray()/1000000.0)

        if stepSize != 0:
            latVRT = latVRT.resized(int(XSize/stepSize), int(YSize/stepSize))
            lonVRT = lonVRT.resized(int(XSize/stepSize), int(YSize/stepSize))

        # Add geolocation domain metadata to the dataset
        self.add_geolocation(Geolocation(xVRT=lonVRT,
                              yVRT=latVRT,
                              xBand=1, yBand=1,
                              srs=srs,
                              lineOffset=geolocParameter[1],
                              lineStep=int(geolocParameter[3]),
                              pixelOffset=geolocParameter[0],
                              pixelStep=int(geolocParameter[2])))


#m = MERIS()
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_1PNPDK20110817_110451_000004053105_00339_49491_7010.N1', range(7, 22))
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(7, 20) + [20, 21, 22, 20])
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(33, 46) + [46, 47, 48, 46])
