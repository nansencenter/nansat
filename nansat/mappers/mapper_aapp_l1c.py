# Name:         mapper_aapp_L1C.py
# Purpose:      Nansat mapper for AAPP Level 1C output
# Author:       Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

# Description of file format:
# http://research.metoffice.gov.uk/research/interproj/nwpsaf/aapp/NWPSAF-MF-UD-003_Formats.pdf (page 120-)
import sys
import struct
import datetime
import warnings

from nansat.exceptions import WrongMapperError
from nansat.geolocation import Geolocation
from nansat.vrt import VRT

dataFormats = {1: 'LAC', 2: 'GAC', 3: 'HRPT'}
dataSetQualityIndicatorOffset = 114
recordLength = 29808
headerLength = recordLength
imageOffset = headerLength + 1092


class Mapper(VRT):
    ''' VRT with mapping of WKV for AVHRR L1C output from AAPP '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):

        ########################################
        # Read metadata from binary file
        ########################################
        try:
            fp = open(filename, 'rb')
        except IOError:
            raise WrongMapperError
        fp.seek(24)
        try:
            satID = int(struct.unpack('<l', fp.read(4))[0])
        except:
            raise WrongMapperError

        ##################
        # Read time
        ##################
        fp.seek(44)
        year = int(struct.unpack('<l', fp.read(4))[0])
        dayofyear = int(struct.unpack('<l', fp.read(4))[0])
        millisecondsOfDay = int(struct.unpack('<l', fp.read(4))[0])
        try:
            time = (datetime.datetime(year, 1, 1) +
                    datetime.timedelta(dayofyear - 1,
                                       milliseconds=millisecondsOfDay))
        except:
            raise WrongMapperError

        fp.seek(72)
        numScanLines = int(struct.unpack('<l', fp.read(4))[0])
        missingScanLines = int(struct.unpack('<l', fp.read(4))[0])
        numCalibratedScanLines = int(struct.unpack('<l', fp.read(4))[0])
        if missingScanLines != 0:
            print('WARNING: Missing scanlines: ' + str(missingScanLines))

        fp.seek(88)
        dataFormatNum = int(struct.unpack('<l', fp.read(4))[0])
        dataFormat = dataFormats[dataFormatNum]

        # Determine if we have channel 3A (daytime) or channel 3B (nighttime)
        def int2bitstring(s):
            return str(s) if s <= 1 else int2bitstring(s >> 1) + str(s & 1)
        fp.seek(headerLength + 20)
        scanlinebitFirstline = int(struct.unpack('<L', fp.read(4))[0])
        fp.seek(headerLength + recordLength*(numCalibratedScanLines-2) + 20)
        scanlinebitLastline = int(struct.unpack('<L', fp.read(4))[0])

        if int2bitstring(scanlinebitFirstline)[-1] == '0':
            startsWith3A = True
        else:
            startsWith3A = False
        if int2bitstring(scanlinebitLastline)[-1] == '0':
            endsWith3A = True
        else:
            endsWith3A = False

        if startsWith3A != endsWith3A:
            print('############################################')
            print('WARNING: channel 3 switches ')
            print('between daytime and nighttime (3A <-> 3B)')
            print('###########################################')

        ###########################
        # Make Geolocation Arrays
        ###########################
        srcRasterYSize = numCalibratedScanLines

        # Making VRT with raw (unscaled) lon and lat
        # (smaller bands than full dataset)
        self.band_vrts = {'RawGeolocVRT': VRT(srcRasterXSize=51,
                                            srcRasterYSize=srcRasterYSize)}
        RawGeolocMetaDict = []
        for lonlatNo in range(1, 3):
            RawGeolocMetaDict.append(
                {'src': {'SourceFilename': filename,
                         'SourceBand': 0,
                         'SourceType': "RawRasterBand",
                         'DataType': gdal.GDT_Int32,
                         'ImageOffset': (headerLength + 676 +
                                         (lonlatNo - 1) * 4),
                         'PixelOffset': 8,
                         'LineOffset': recordLength,
                         'ByteOrder': 'LSB'},
                 'dst': {}})

        self.band_vrts['RawGeolocVRT'].create_bands(RawGeolocMetaDict)

        # Make derived GeolocVRT with scaled lon and lat
        self.band_vrts['GeolocVRT'] = VRT(srcRasterXSize=51,
                                        srcRasterYSize=srcRasterYSize)
        GeolocMetaDict = []
        for lonlatNo in range(1, 3):
            GeolocMetaDict.append(
                {'src': {'SourceFilename': (self.band_vrts['RawGeolocVRT'].
                                            filename),
                         'SourceBand': lonlatNo,
                         'ScaleRatio': 0.0001,
                         'ScaleOffset': 0,
                         'DataType': gdal.GDT_Int32},
                 'dst': {}})

        self.band_vrts['GeolocVRT'].create_bands(GeolocMetaDict)

        GeolocObject = Geolocation(x_vrt=self.band_vrts['GeolocVRT'],
                                        y_vrt=self.band_vrts['GeolocVRT'],
                                        x_band=2, y_band=1,  # x = lon, y = lat
                                        line_offset=0, pixel_offset=25,
                                        line_step=1, pixel_step=40)

        #######################
        # Initialize dataset
        #######################
        # create empty VRT dataset with geolocation only
        # (from Geolocation Array)
        self._init_from_dataset_params(2048, numCalibratedScanLines,
                                        (0,1,0,numCalibratedScanLines,0,-1),
                                        GeolocObject.d['SRS'])
        self._add_geolocation(GeolocObject)


        ##################
        # Create bands
        ##################
        self.band_vrts['RawBandsVRT'] = VRT(
            srcRasterXSize=2048,
            srcRasterYSize=numCalibratedScanLines)
        RawMetaDict = []
        metaDict = []

        centralWavelengths = [0.63, 0.86, np.NaN, 10.8, 12.0]
        if startsWith3A:
            centralWavelengths[2] = 1.6
            firstIRband = 4
        else:
            centralWavelengths[2] = 3.7
            firstIRband = 3

        for bandNo in range(1, 6):
            RawMetaDict.append(
                {'src': {'SourceFilename': filename,
                         'SourceBand': 0,
                         'SourceType': "RawRasterBand",
                         'dataType': gdal.GDT_UInt16,
                         'ImageOffset': imageOffset + (bandNo - 1) * 2,
                         'PixelOffset': 10,
                         'LineOffset': recordLength,
                         'ByteOrder': 'LSB'},
                 'dst': {'dataType': gdal.GDT_UInt16}})

            if bandNo < firstIRband:
                wkv = 'albedo'
                minmax = '0 60'
            else:
                wkv = 'brightness_temperature'
                minmax = '290 210'

            metaDict.append(
                {'src': {'SourceFilename': (self.band_vrts['RawBandsVRT'].
                                            filename),
                         'SourceBand': bandNo,
                         'ScaleRatio': 0.01,
                         'ScaleOffset': 0,
                         'DataType': gdal.GDT_UInt16},
                 'dst': {'originalFilename': filename,
                         'dataType': gdal.GDT_Float32,
                         'wkv': wkv,
                         'colormap': 'gray',
                         'wavelength': centralWavelengths[bandNo-1],
                         'minmax': minmax}})

        # Add temperature difference between ch3 and ch 4 as pixelfunction
        if not startsWith3A:  # Only if ch3 is IR (nighttime)
            metaDict.append(
                {'src': [{'SourceFilename': (self.band_vrts['RawBandsVRT'].
                                             filename),
                          'ScaleRatio': 0.01,
                          'ScaleOffset': 0,
                          'SourceBand': 4},
                         {'SourceFilename': (self.band_vrts['RawBandsVRT'].
                                             filename),
                          'ScaleRatio': 0.01,
                          'ScaleOffset': 0,
                          'SourceBand': 3}],
                 'dst': {'PixelFunctionType': 'diff',
                         'originalFilename': filename,
                         'dataType': gdal.GDT_Float32,
                         'name': 'ch4-ch3',
                         'short_name': 'ch4-ch3',
                         'long_name': 'AVHRR ch4 - ch3 temperature difference',
                         'colormap': 'gray',
                         'units': 'kelvin',
                         'minmax': '-3 3'}})

        self.self.band_vrts['RawBandsVRT'].create_bands(RawMetaDict)
        self.create_bands(metaDict)

        globalMetadata = {}
        globalMetadata['satID'] = str(satID)
        globalMetadata['daytime'] = str(int(startsWith3A))
        self.dataset.SetMetadata(globalMetadata)

        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start', time.isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', time.isoformat())

        return
