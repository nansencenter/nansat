# Name:         mapper_asar_doppler.py
# Purpose:      Mapper for ASAR data to read the Doppler centroid grid in ASAR WSM data
# Authors:      Peng Yu, Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from dateutil.parser import parse
import struct
from nansat.vrt import VRT, GeolocationArray
import gdal
import numpy as np
from scipy import interpolate

class Doppler():
    '''This class is similar to the ENVISAT mapper
    We modify the envisat mapper to read doppler information from ADS (additional data sets)
    which are DOP CENTROID_GRID_ADS and GEOLOCATION_GRID_ADS (ASAR)
    You can see the method and data direction in envisat.py 
    '''
    allADSParams = {
        'ASA_': {
            'name': 'DS_NAME="GEOLOCATION GRID ADS        "\n',
            'width' : 11,
            'list': {
                "num_lines"                   : {"offset": 13                 , "dataType": gdal.GDT_Int16  , "units": ""},
                "first_line_samp_numbers"     : {"offset": 25+11*4*0          , "dataType": gdal.GDT_Float32, "units": ""},
                "first_line_slant_range_times": {"offset": 25+11*4*1          , "dataType": gdal.GDT_Float32, "units": "ns"},
                "first_line_incidence_angle"  : {"offset": 25+11*4*2          , "dataType": gdal.GDT_Float32, "units": "deg"},
                "first_line_lats"             : {"offset": 25+11*4*3          , "dataType": gdal.GDT_Int32  , "units": "(10)^-6 deg"},
                "first_line_longs"            : {"offset": 25+11*4*4          , "dataType": gdal.GDT_Int32  , "units": "(10)^-6 deg"},
                "last_line_samp_numbers"      : {"offset": 25+11*4*5+34+11*4*0, "dataType": gdal.GDT_Int32  , "units": ""},
                "last_line_slant_range_times" : {"offset": 25+11*4*5+34+11*4*1, "dataType": gdal.GDT_Float32, "units": "ns"},
                "last_line_incidence_angle"   : {"offset": 25+11*4*5+34+11*4*2, "dataType": gdal.GDT_Float32, "units": "deg"},
                "last_line_lats"              : {"offset": 25+11*4*5+34+11*4*3, "dataType": gdal.GDT_Int32  , "units": "(10)^-6 deg"},
                "last_line_longs"             : {"offset": 25+11*4*5+34+11*4*4, "dataType": gdal.GDT_Int32  , "units": "(10)^-6 deg"}
            }
        },
	'zdt': {
            'name'  : 'DS_NAME="DOP CENTROID GRID ADS       "\n',
	    'width' : 1,
	    'list': {
                "first_zero_doppler_time": {"offset": 0, "dataType": gdal.GDT_Int32, "units": "MJD"},
		"last_zero_doppler_time" : {"offset": 13+100*4*2, "dataType": gdal.GDT_Int32, "units": "MJD"}
	    }
        },
	'DC': {
            'name': 'DS_NAME="DOP CENTROID GRID ADS       "\n',
            'width' : 100,
            'list': {
                "slant_range_time" : {"offset": 13 , "dataType": gdal.GDT_Float32 , "units": "ns"},
		"dop_coef" : {"offset": 13+100*4*1 , "dataType": gdal.GDT_Float32 , "units": "Hz"}
            }
        }
    }

    # map: GDAL TYPES ==> struct format strings
    structFmt = {gdal.GDT_Int16: ">h",
                 gdal.GDT_UInt16: ">H",
                 gdal.GDT_Int32: ">i",
                 gdal.GDT_UInt32: ">I",
                 gdal.GDT_Float32: ">f"}
    # names of grids with longitude/latitude in ASAR and MERIS ADS
    lonlatNames = {'ASA_': ['first_line_longs', 'first_line_lats']}

    def __init__(self, fileName):
        '''Select set of params'''
        self.iFileName = fileName
        self.lonlatNames = self.lonlatNames['ASA_']

    def _set_envisat_time(self, gdalMetadata):
        ''' Get time from metadata, set time to VRT'''
        # set time
        productTime = gdalMetadata["SPH_FIRST_LINE_TIME"]
        self._set_time(parse(productTime))

    def add_zdtimeADSParams(self):
        #read First zero doppler time offset
        self.dsOffsetDict = {'zdt':
                self.read_offset_from_header(self.allADSParams['zdt']['name'])}
		   
    def add_dopplerADSParams(self):
        #read doppler offset
        self.dsOffsetDict = {'DC':
                self.read_offset_from_header(self.allADSParams['DC']['name'])}
		
    def add_asarADSParams(self):
        #read asar offset
        self.dsOffsetDict = {'ASA_': self.read_offset_from_header(self.allADSParams['ASA_']['name'])}
 
    def read_offset_from_header(self, gadsDSName):
        # number of lines after line with 'DS_NAME'
        textOffset = {'DS_OFFSET': 3, 'DS_SIZE': 4,
                      'NUM_DSR': 5, 'DSR_SIZE': 6}

        # open file and read 150 header lines
        f = file(self.iFileName, 'rt')
        headerLines = f.readlines(150)
        offsetDict = {}
        # create a dictionary with offset, size, number of records,
        # size of records
        if gadsDSName in headerLines:
            # get location of gadsDSName
            gridOffset = headerLines.index(gadsDSName)
            # Adjust the location of the varaibles by adding textOffset.
            # Read a text at the location and convert the text into integer.
            for iKey in textOffset:
                offsetDict[iKey] = int(headerLines[gridOffset +
                                       textOffset[iKey]].replace(iKey+"=", '').
                                       replace('<bytes>', ''))
        f.close()
        return offsetDict

    def read_binary_line(self, offset, fmtString, length):
        # fseek, read all into a list
        f = file(self.iFileName, 'rb')
        f.seek(offset, 0)
        binaryValues = []
        for i in range(length):
            fString = f.read(struct.calcsize(fmtString))
            fVal = struct.unpack(fmtString, fString)[0]
            binaryValues.append(fVal)
        f.close()
        
        return binaryValues

    def resize_dop_array(self, arrayHeight, arrayWidth, array):
        """resize the x*y array to doppler x*y array"""
        #set doppler array size
        dopHeight = self.read_offset_from_header('DS_NAME="DOP CENTROID GRID ADS       "\n')["NUM_DSR"]
        dopWidth = self.allADSParams['DC']['width']
        dopInc = np.zeros((arrayHeight, arrayWidth))
        x = np.arange(arrayHeight).reshape(1,arrayHeight)
        y = np.arange(arrayWidth).reshape(1,arrayWidth)
        xx = np.linspace(x.min(),x.max(),dopHeight)
        yy = np.linspace(y.min(),y.max(),dopWidth)
        resizeArray = interpolate.RectBivariateSpline(x,y,array, kx=3,ky=3)
        newArray = resizeArray(xx,yy)
        dopArray = newArray
        return dopArray

    def create_VRT_from_ADS(self, adsName, keyName=None):
        # Get parameters of arrays in ADS
        adsWidth = self.allADSParams[keyName]['width']
        adsParams = self.allADSParams[keyName]['list'][adsName]
        dsOffsetDict = self.dsOffsetDict[keyName]
        adsHeight = dsOffsetDict["NUM_DSR"]
        fmtString = self.structFmt[adsParams['dataType']]
        #set the size of doppler array
        dopWidth = self.allADSParams['DC']['width']
        dopHeight = self.read_offset_from_header('DS_NAME="DOP CENTROID GRID ADS       "\n')["NUM_DSR"]
        # create an array whose elements are fetched from ADS
        array = np.array([])

        for i in range(adsHeight):
            lineOffset = (dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] +
                          dsOffsetDict["DSR_SIZE"] * i)
            binaryLine = self.read_binary_line(lineOffset, fmtString, adsWidth)
            array = np.append(array, binaryLine)
        array = array.reshape(adsHeight, adsWidth)

        # read 'last_line_...'
        if keyName == 'ASA_':
            adsName = adsName.replace('first_line', 'last_line')
            adsParams = self.allADSParams[keyName]['list'][adsName]
            lineOffset = (dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] +
                          dsOffsetDict["DSR_SIZE"] * i)
            binaryLine = self.read_binary_line(lineOffset, fmtString, adsWidth)
            array = np.append(array, binaryLine)
            adsHeight += 1       
            array = array.reshape(adsHeight, adsWidth)
            #array = self.resize_dop_array(adsHeight, adsWidth, array)

        # adjust the scale
        if '(10)^-6' in adsParams['units']:
            array /= 1000000.0
            adsParams['units'] = adsParams['units'].replace('(10)^-6 ', '')

        # create VRT from the array
        adsVrt = VRT(array=array)
        # add "name" and "units" to band metadata
        bandMetadata = {"name": adsName, "units": adsParams['units']}
        adsVrt.dataset.GetRasterBand(1).SetMetadata(bandMetadata)
        return adsVrt

    def get_ads_vrts(self, gdalDataset, adsNames, keyName=None,
                     step=1, **kwargs):
        dopHeight = self.read_offset_from_header('DS_NAME="DOP CENTROID GRID ADS       "\n')["NUM_DSR"]
        dopWidth = self.allADSParams['DC']['width']
        # list with VRT with arrays of lon/lat
        adsVRTs = []
        for adsName in adsNames:
            # create VRT with array from ADS
            adsVRTs.append(self.create_VRT_from_ADS(adsName, keyName))
            # resize the VRT to match <step>
            adsVRTs[-1] = adsVRTs[-1].get_resized_vrt(dopWidth/step, dopHeight/step, **kwargs)
        return adsVRTs

    def set_zero_doppler_time(self, gdalDataset):
        if gdalDataset.GetMetadataItem('MPH_SOFTWARE_VER')[5:9] > 4.06:
            fzdtArray = self.get_zero_doppler_time('first_zero_doppler_time')
            lzdtArray = self.get_zero_doppler_time('last_zero_doppler_time')
            #check that the difference between first and last zdt is correct
            #Time increment from first to last
            delta = lzdtArray - fzdtArray 
            #New centred zero Doppler times
            zdtArray = fzdtArray + delta/2
        else:
            fzdtArray = self.get_zero_doppler_time('first_zero_doppler_time')
            zdtArray = fzdtArray  + float(6)*0.192/86400
            lzdtArray  = zdtArray + float(6)*0.192/86400
        #array from dopHeight*1 to  dopHeight*100
        #arrayWidth = np.ones((1,100))
        #zdtArray  = zdtArray*arrayWidth
        # create VRT from the array
        adsVrt = VRT(array=zdtArray)
        # add "name" and "units" to band metadata
        bandMetadata = {"name": 'zero doppler time', "units": 'MJD'}
        adsVrt.dataset.GetRasterBand(1).SetMetadata(bandMetadata)
        adsVRTs = []
        adsVRTs.append(adsVrt)
        return adsVRTs

    def read_time_value(self, offset, fmtString, length):
        # fseek, read the number value
        f = file(self.iFileName, 'rb')
        f.seek(offset, 0)
        for i in range(length):
            fString = f.read(struct.calcsize(fmtString))
            fVal = struct.unpack(fmtString, fString)[0]
        f.close()
        return fVal
        
    def get_zero_doppler_time(self, adsName):
        adsParams = self.allADSParams['zdt']['list'][adsName]
        dsOffsetDict = self.dsOffsetDict['Zrdt']
        adsHeight = dsOffsetDict["NUM_DSR"]
        # create day, second, micronsecond array
        zdtArray = np.array([],dtype='>f')
        for j in range(adsHeight):
            dayOffset = (dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] +
                          dsOffsetDict["DSR_SIZE"] * j)
            dayValue = self.read_time_value(dayOffset,'>i', 1)

            secondOffset = (dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] + 4 +
                          dsOffsetDict["DSR_SIZE"] * j)
            secondValue = self.read_time_value(secondOffset,'>I', 1)

            micsecondOffset = (dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] + 8 +
                          dsOffsetDict["DSR_SIZE"] * j)
            micsecondValue = self.read_time_value(micsecondOffset,'>I', 1)
            # create fzd time line with second as the time unit
            #zdtValue = float(dayValue)*24*3600 + float(secondValue)+ float(micsecondValue)/1000000
            #dayValue is ignored for higher second precison
            #one can caculate the days using _set_envisat_time or dayValue plus 01,01,2000
            zdtValue = float(secondValue)+ float(micsecondValue)/1000000
            zdtLine = [zdtValue]*100
            zdtArray = np.append(zdtArray, zdtLine)
        # create fzd time array with second as the time unit
        zdtArray = zdtArray.reshape(adsHeight, 100)
        return zdtArray


from nansat.domain import Domain

class Mapper(VRT, Doppler):
    ''' VRT with mapping of WKV for ASAR Level 0

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR   and
            https://earth.esa.int/handbooks/asar/CNTR6-6-25.htm#eph.asar.asardf.asarrec.ASAR_Dop_Cen_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 full_incAng=True, geolocation=False,step=1, **kwargs):
        '''
        Parameters
        -----------
        fileName : string
        gdalDataset : gdal dataset
        gdalMetadata : gdal metadata
        full_incAng : bool (default is True)
            if True, add full size incedence angle
        geolocation : bool (default is False)
            if True, add gdal geolocation
        step: int (used in class(Doppler))
            step of pixel and line in GeolocationArrays. lat/lon grids are
            generated at that step
        '''
        product = gdalMetadata.get("MPH_PRODUCT")
        if product[0:6] != "ASA_WS":
            raise AttributeError("BAD ASA_WSM MAPPER")

        Doppler.__init__(self, fileName)

        # get polarization string (remove '/', since NetCDF doesnt support that in metadata)
        polarization = gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", "")

       # Create VRTdataset with small VRTRawRasterbands
        self.add_zdtimeADSParams()
        self.subVRTs = {'adsVRTs' : self.set_zero_doppler_time(gdalDataset)}

        self.add_dopplerADSParams()
        for iADSName in ["slant_range_time", "dop_coef"]:
            self.subVRTs['adsVRTs'].append(self.get_ads_vrts(
                                            gdalDataset,
                                            [iADSName],
                                            keyName = 'DC',
                                            step=step, **kwargs)[0])
    
        self.add_asarADSParams()
        self.subVRTs['adsVRTs'].append(self.get_ads_vrts(gdalDataset,
                                        ["first_line_incidence_angle"],
                                        keyName = 'ASA_',
                                        step=step, **kwargs)[0])

        # create empty VRT dataset with geolocation only and resize the RawCounts array
        VRT.__init__(self, gdalDataset)
        dopHeight = self.read_offset_from_header('DS_NAME="DOP CENTROID GRID ADS       "\n')["NUM_DSR"]
        dopWidth = self.allADSParams['DC']['width']
        rawBand = gdalDataset.GetRasterBand(1)
        rawArray = rawBand.ReadAsArray()
        doprawArray = self.resize_dop_array(gdalDataset.RasterYSize, gdalDataset.RasterXSize, rawArray)
        RawadsVrt = VRT(array=doprawArray)
        # add "name" and "units" to band metadata
        bandMetadata = {"name": 'RawCounts'}
        RawadsVrt.dataset.GetRasterBand(1).SetMetadata(bandMetadata)
        RawadsVrt= RawadsVrt.get_resized_vrt(dopWidth/step, dopHeight/step, **kwargs)
        self.subVRTs['adsVRTs'].append(RawadsVrt)

        # get calibration constant
        gotCalibration = True
        try:
            calibrationConst = float(gdalDataset.GetMetadataItem(
                "MAIN_PROCESSING_PARAMS_ADS_CALIBRATION_FACTORS.1.EXT_CAL_FACT", "records"))
        except:
            try:
                # Apparently some ASAR files have calibration constant stored in another place
                calibrationConst = float(gdalDataset.GetMetadataItem(
                    "MAIN_PROCESSING_PARAMS_ADS_1_CALIBRATION_FACTORS.1.EXT_CAL_FACT", "records"))
            except:
                self.logger.warning('Cannot get calibrationConst')
                gotCalibration = False

        # add empty dictionary
        metaDict = []

        if full_incAng:
            for adsVRT in self.subVRTs['adsVRTs']:
                metaDict.append({'src': {'SourceFilename': adsVRT.fileName,
                                         'SourceBand': 1},
                                 'dst': {'name': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name').replace('last_line_', ''),
                                         'short_name': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name').replace('last_line_', ''),
                                         'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}})
        if gotCalibration:
                # add dicrtionary for sigma0
                sphPass = gdalMetadata['SPH_PASS']
                metaDict.append({'src': [{'SourceBand': 1, 'ScaleRatio': np.sqrt(1.0 /calibrationConst), 
                                                             'SourceFilename': self.subVRTs['adsVRTs'][-1].fileName},
                                                            {'SourceBand': 1, 'SourceFilename': self.subVRTs['adsVRTs'][-1].fileName }],
                                                 'dst': {'short_name': 'sigma0',
                                                            'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                                                            'PixelFunctionType': 'RawcountsIncidenceToSigma0',
                                                            'polarization': polarization,
                                                            'suffix': polarization,
                                                            'pass': sphPass,
                                                            'dataType': 6}})
        #add GeolocationArray
        xyVRTs = self.get_ads_vrts(gdalDataset, self.lonlatNames, 'ASA_',
                                       step)
        GeoloArray = GeolocationArray(xVRT=xyVRTs[0],
                                        yVRT=xyVRTs[1],
                                        xBand=1, yBand=1,
                                        srs=gdalDataset.GetGCPProjection(),
                                        lineOffset=0, 
                                        lineStep=1,
                                        pixelOffset=0,
                                        pixelStep=1)

        # create empty VRT dataset using GeolocationArray only
        VRT.__init__(self,
                     srcRasterXSize=dopWidth,
                     srcRasterYSize=dopHeight,
                     geolocationArray=GeoloArray,
                     srcProjection=GeoloArray.d['SRS'])

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays
        # get VRTs with lon and lat
        if geolocation:
            xyVRTs = self.get_ads_vrts(gdalDataset, self.lonlatNames, 'ASA_',
                                       step)

            # Add geolocation domain metadata to the dataset
            self.add_geolocationArray(GeolocationArray(xVRT=xyVRTs[0],
                                      yVRT=xyVRTs[1],
                                      xBand=1, yBand=1,
                                      srs=gdalDataset.GetGCPProjection(),
                                      lineOffset=0,
                                      lineStep=1,
                                      pixelOffset=0,
                                      pixelStep=step))

        # Add SAR look direction to metadata domain
        self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT') # ASAR is always right-looking
        self.dataset.SetMetadataItem('ORBIT_DIRECTION', gdalMetadata['SPH_PASS'].upper())

        # "SAR_center_look_direction" below is obsolete, and may soon be deleted
        #
        # Note that this is the look direction in the center of the domain. For
        # longer domains, especially at high latitudes, the azimuth direction
        # may vary a lot over the domain, and using the center angle will be a
        # coarse approximation.
        self.dataset.SetMetadataItem('SAR_center_look_direction',
                                     str(np.mod(Domain(ds=gdalDataset).
                                         upwards_azimuth_direction() + 90,
                                                360)))

        ###################################################################
        # Add sigma0_VV
        ###################################################################
        if 'VV' not in polarization and 'HH' in polarization:
            src = [{'SourceBand': 1, 'ScaleRatio': np.sqrt(1.0 /calibrationConst),
                   'SourceFilename': self.subVRTs['adsVRTs'][-1].fileName},
                   {'SourceBand': 1, 
                   'SourceFilename': self.subVRTs['adsVRTs'][-1].fileName }]
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sigma0HHToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self._create_band(src, dst)
            self.dataset.FlushCache()
