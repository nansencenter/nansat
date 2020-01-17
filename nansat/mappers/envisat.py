# Name:         envisat.py
# Purpose:      Contains Envisat class definition
# Authors:      Asuka Yamakava, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse
import struct

import numpy as np
try:
    import scipy.ndimage
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

from nansat.vrt import VRT
from nansat.geolocation import Geolocation
from nansat.utils import gdal, ogr

from nansat.exceptions import WrongMapperError, NansatReadError

class Envisat(object):
    """Methods/data shared between Envisat mappers

    This class is needed to read awkward N1 format of ENVISAT
    Mostly it support reading of variables from ADS (additional data sets)
    which are TIE_POINTS_ADS (MERIS) or GEOLOCATION_GRID_ADS (ASAR)
    """

    # specific name of geolocation and offsets,
    # dataTypes and units of each dataset
    allADSParams = {
        'MER_': {'name': 'DS_NAME="Tie points ADS              "\n',
                 'width': 71,
                 'list': {"latitude":
                          {
                          "offset": 13,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "longitude":
                          {
                          "offset": 13+284*1,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "DME altitude":
                          {
                          "offset": 13+284*2,
                          "dataType": gdal.GDT_Int32,
                          "units": "m"
                          },
                          "DME roughness":
                          {
                          "offset": 13+284*3,
                          "dataType": gdal.GDT_UInt32,
                          "units": "m"
                          },
                          "DME latitude corrections":
                          {
                          "offset": 13+284*4,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "DME longitude corrections":
                          {
                          "offset": 13+284*5,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "sun zenith angles":
                          {
                          "offset": 13+284*6,
                          "dataType": gdal.GDT_UInt32,
                          "units": "(10)^-6 deg"
                          },
                          "sun azimuth angles":
                          {
                          "offset": 13+284*7,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "viewing zenith angles":
                          {
                          "offset": 13+284*8,
                          "dataType": gdal.GDT_UInt32,
                          "units": "(10)^-6 deg"
                          },
                          "viewing azimuth angles":
                          {
                          "offset": 13+284*9,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "zonal winds":
                          {
                          "offset": 13+284*10+142*0,
                          "dataType": gdal.GDT_Int16,
                          "units": "m*s-1"
                          },
                          "meridional winds":
                          {
                          "offset": 13+284*10+142*1,
                          "dataType": gdal.GDT_Int16,
                          "units": "m*s-1"
                          },
                          "mean sea level pressure":
                          {
                          "offset": 13+284*10+142*2,
                          "dataType": gdal.GDT_UInt16,
                          "units": "hPa"
                          },
                          "total ozone":
                          {
                          "offset": 13+284*10+142*3,
                          "dataType": gdal.GDT_UInt16,
                          "units": "DU"
                          },
                          "relative humidity":
                          {
                          "offset": 13+284*10+142*4,
                          "dataType": gdal.GDT_UInt16,
                          "units": "%"
                          }
                          }},
        'ASA_': {'name': 'DS_NAME="GEOLOCATION GRID ADS        "\n',
                 'width': 11,
                 'list': {"num_lines":
                          {
                          "offset": 13,
                          "dataType": gdal.GDT_Int16,
                          "units": ""
                          },
                          "first_line_samp_numbers":
                          {
                          "offset": 25+11*4*0,
                          "dataType": gdal.GDT_Float32,
                          "units": ""
                          },
                          "first_line_slant_range_times":
                          {
                          "offset": 25+11*4*1,
                          "dataType": gdal.GDT_Float32,
                          "units": "ns"
                          },
                          "first_line_incidence_angle":
                          {
                          "offset": 25+11*4*2,
                          "dataType": gdal.GDT_Float32,
                          "units": "deg"
                          },
                          "first_line_lats":
                          {
                          "offset": 25+11*4*3,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "first_line_longs":
                          {
                          "offset": 25+11*4*4,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "last_line_samp_numbers":
                          {
                          "offset": 25+11*4*5+34+11*4*0,
                          "dataType": gdal.GDT_Int32,
                          "units": ""
                          },
                          "last_line_slant_range_times":
                          {
                          "offset": 25+11*4*5+34+11*4*1,
                          "dataType": gdal.GDT_Float32,
                          "units": "ns"
                          },
                          "last_line_incidence_angle":
                          {
                          "offset": 25+11*4*5+34+11*4*2,
                          "dataType": gdal.GDT_Float32,
                          "units": "deg"
                          },
                          "last_line_lats":
                          {
                          "offset": 25+11*4*5+34+11*4*3,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          },
                          "last_line_longs":
                          {
                          "offset": 25+11*4*5+34+11*4*4,
                          "dataType": gdal.GDT_Int32,
                          "units": "(10)^-6 deg"
                          }
                          }
                 }}

    # map: GDAL TYPES ==> struct format strings
    structFmt = {gdal.GDT_Int16: ">h",
                 gdal.GDT_UInt16: ">H",
                 gdal.GDT_Int32: ">i",
                 gdal.GDT_UInt32: ">I",
                 gdal.GDT_Float32: ">f"}
    # names of grids with longitude/latitude in ASAR and MERIS ADS
    lonlatNames = {'ASA_': ['first_line_longs', 'first_line_lats'],
                   'MER_': ['longitude', 'latitude']}

    def setup_ads_parameters(self, filename, gdalMetadata):
        """Select set of params and read offset of ADS"""
        if not gdalMetadata or not ('MPH_PRODUCT' in gdalMetadata.keys()):
            raise WrongMapperError

        self.product = gdalMetadata["MPH_PRODUCT"]
        self.iFileName = filename
        self.prodType = gdalMetadata["MPH_PRODUCT"][0:4]
        self.allADSParams = self.allADSParams[self.prodType]
        self.dsOffsetDict = self.read_offset_from_header(
                                                    self.allADSParams['name'])
        self.lonlatNames = self.lonlatNames[self.prodType]

    def _set_envisat_time(self, gdalMetadata):
        """ Get time from metadata, set time to VRT"""
        # set valid time
        self.dataset.SetMetadataItem('time_coverage_start', parse(gdalMetadata["SPH_FIRST_LINE_TIME"]).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', parse(gdalMetadata["SPH_LAST_LINE_TIME"]).isoformat())


    def read_offset_from_header(self, gadsDSName):
        """ Read offset of ADS from text header.

        Find a location of gadsDSName.
        Adjust the location with textOffset and read the text at the location.
        Convert the text to integer and set it into offsetDict.

        Returns
        -------
            offsetDict : dictionary
                offset of DS, size of DS, number of records, size of record
        """
        # number of lines after line with 'DS_NAME'
        textOffset = {'DS_OFFSET': 3, 'DS_SIZE': 4,
                      'NUM_DSR': 5, 'DSR_SIZE': 6}

        # open file and read 150 header lines
        with open(self.iFileName, "rb") as f:
            headerLines = [f.readline().decode('utf-8') for i in range(150)]

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
        return offsetDict

    def read_binary_line(self, offset, fmtString, length):
        """Read line with binary data at given offset

        Open file
        Read values of given size, at given offset, given times
        Convert to list of values of given format and return

        Parameters
        ----------
            offset: int, start of reading
            fmtString: str, data type format
            lentgh: number of values to read

        Returns
        -------
            binaryValues : list
                values which are read from the file.
                the number of elements is length
        """
        # fseek, read all into a list
        f = open(self.iFileName, 'rb')
        f.seek(offset, 0)
        binaryValues = []
        for i in range(length):
            fString = f.read(struct.calcsize(fmtString))
            fVal = struct.unpack(fmtString, fString)[0]
            binaryValues.append(fVal)
        f.close()

        return binaryValues

    def read_scaling_gads(self, indeces):
        """ Read Scaling Factor GADS to get scalings of MERIS L1/L2

        Parameters
        ----------
            filename : string
            indeces : list

        Returns
        -------
            list
        """
        maxGADS = max(indeces) + 1
        dsOffsetDict = self.read_offset_from_header(
            'DS_NAME="Scaling Factor GADS         "\n')
        allGADSValues = self.read_binary_line(dsOffsetDict["DS_OFFSET"],
                                              '>f', maxGADS)
        #get only values required for the mapper
        return [allGADSValues[i] for i in indeces]

    def get_array_from_ADS(self, adsName):
        """ Create VRT with a band from Envisat ADS metadata

        Read offsets of the <adsName> ADS.
        Read 2D matrix of binary values from ADS from file.
        Read last line ADS (in case of ASAR).
        Zoom array with ADS data to <zoomSize>. Zooming is needed to create
        smooth matrices. Array is zoomed to small size because it is stred in
        memory. Later the VRT with zoomed array is VRT.get_resized_vrt()
        in order to match the size of the Nansat onject.
        Create VRT from the ADS array.

        Parameters
        ----------
            adsName : str
                name of variable from ADS to read. should match allADSParams

        Returns
        -------
            adsVrt : VRT
                vrt with a band created from ADS array

        """
        # Get parameters of arrays in ADS
        adsWidth = self.allADSParams['width']
        adsParams = self.allADSParams['list'][adsName]

        # get data type format string and size
        fmtString = self.structFmt[adsParams['dataType']]

        # create an array whose elements are fetched from ADS
        array = np.array([])

        # read sequence of 1D arrays from ADS
        adsHeight = self.dsOffsetDict["NUM_DSR"]
        for i in range(adsHeight):
            lineOffset = (self.dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] +
                          self.dsOffsetDict["DSR_SIZE"] * i)
            binaryLine = self.read_binary_line(lineOffset, fmtString, adsWidth)
            array = np.append(array, binaryLine)

        # read 'last_line_...'
        if self.prodType == 'ASA_':
            adsName = adsName.replace('first_line', 'last_line')
            adsParams = self.allADSParams['list'][adsName]
            lineOffset = (self.dsOffsetDict['DS_OFFSET'] +
                          adsParams['offset'] +
                          self.dsOffsetDict["DSR_SIZE"] * i)
            binaryLine = self.read_binary_line(lineOffset, fmtString, adsWidth)
            array = np.append(array, binaryLine)
            adsHeight += 1

        # adjust the scale
        if '(10)^-6' in adsParams['units']:
            array /= 1000000.0
            # Commenting out line below, otherwise subsequent calls
            # for lon and lat results in modified unit,
            # and hence necessary scaling is not performed
            # adsParams['units'] = adsParams['units'].replace('(10)^-6 ', '')

        # reshape the array into 2D matrix
        array = array.reshape(adsHeight, adsWidth)
        return array

    def create_VRT_from_ADS(self, adsName, zoomSize=500):
        """ Create VRT with a band from Envisat ADS metadata

        Read offsets of the <adsName> ADS.
        Read 2D matrix of binary values from ADS from file.
        Zoom array with ADS data to <zoomSize>. Zooming is needed to create smooth matrices. Array
        is zoomed to small size because it is stred in memory. Later the VRT with zoomed array is
        VRT.get_resized_vrt() in order to match the size of the Nansat object.

        Create VRT from the ADS array.

        Parameters
        ----------
            adsName : str
                name of variable from ADS to read. should match allADSParams
            zoomSize :  int, optional, 500
                size, to which original matrix from ADSR is zoomed using
                scipy.zoom

        Returns
        -------
            adsVrt : VRT, vrt with a band created from ADS array

        """
        adsHeight = self.dsOffsetDict["NUM_DSR"]
        adsParams = self.allADSParams['list'][adsName]
        array = self.get_array_from_ADS(adsName)

        if not IMPORT_SCIPY:
            raise NansatReadError('ENVISAT data cannot be read because scipy is not installed...')

        # zoom the array
        array = scipy.ndimage.interpolation.zoom(array,
                                                 zoomSize / float(adsHeight),
                                                 order=1)

        # create VRT from the array
        adsVrt = VRT.from_array(array=array)
        # add "name" and "units" to band metadata
        bandMetadata = {"name": adsName, "units": adsParams['units']}
        adsVrt.dataset.GetRasterBand(1).SetMetadata(bandMetadata)

        return adsVrt

    def get_ads_vrts(self, gdalDataset, adsNames, zoomSize=500,
                     step=1, **kwargs):
        """Create list with VRTs with zoomed and resized ADS arrays

        For given names of variables (which should match self.allADSParams):
            Get VRT with zoomed ADS array
            Get resized VRT

        Parameters
        ----------
            gdalDataset: GDAL Dataset
                input dataset
            adsNames: list with strings
                names of varaiables from self.allADSParams['list']
            zoomSize: int, 500
                size to which the ADS array will be zoomed by scipy.zoom
            step: int, 1
                step, at which data will be given

        Returns
        --------
            adsVRTs: list with VRT
                list with resized VRT with zoomed arrays

        """
        XSize = gdalDataset.RasterXSize
        YSize = gdalDataset.RasterYSize
        # list with VRT with arrays of lon/lat
        adsVRTs = []
        for adsName in adsNames:
            # create VRT with array from ADS
            adsVRTs.append(self.create_VRT_from_ADS(adsName, zoomSize))
            # resize the VRT to match <step>
            adsVRTs[-1] = adsVRTs[-1].get_resized_vrt(XSize/step, YSize/step)
        return adsVRTs

    def add_geolocation_from_ads(self, gdalDataset, zoomSize=500, step=1):
        """
        Add geolocation domain metadata to the dataset

        Get VRTs with zoomed arrays of lon and lat
        Create geolocation object and add to the metadata

        Parameters
        ----------
        gdalDataset: GDAL Dataset
            input dataset
        zoomSize: int, optional, 500
            size, to which the ADS array will be zoomed using scipy
            array of this size will be stored in memory
        step: int
            step of pixel and line in GeolocationArrays. lat/lon grids are
            generated at that step

        Modifies
        --------
        Adds Geolocation Array metadata
        """
        # get VRTs with lon and lat
        xyVRTs = self.get_ads_vrts(gdalDataset, self.lonlatNames, zoomSize,
                                   step)

        # Add geolocation domain metadata to the dataset
        self._add_geolocation(Geolocation(x_vrt=xyVRTs[0],
                                  y_vrt=xyVRTs[1],
                                  x_band=1, y_band=1,
                                  srs=gdalDataset.GetGCPProjection(),
                                  line_offset=0,
                                  line_step=step,
                                  pixel_offset=0,
                                  pixel_step=step))

