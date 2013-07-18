# Name:         mapper_asar.py
# Purpose:      Mapper for Envisat/ASAR data
# Authors:      Asuka Yamakava, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from vrt import VRT
from envisat import Envisat
import numpy as np

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, caliblation='one-line'):

        product = gdalMetadata.get("MPH_PRODUCT")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        # init ADS parameters
        Envisat.__init__(self, fileName, product[0:4])

        # get polarization string (remove '/', since NetCDF doesnt support that in metadata)
        polarization = gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", "")

        # Create VRTdataset with small VRTRawRasterbands
        self.adsVRTs = self.get_ads_vrts(gdalDataset,
                                         ["first_line_incidenceAngle",
                                          "first_line_incidenceAngle"],
                                          lineBand=[True, False])

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # get calibration constant
        gotCalibration = True
        try:
            calibrationConst = float(gdalDataset.GetMetadataItem("MAIN_PROCESSING_PARAMS_ADS_1_CALIBRATION_FACTORS.1.EXT_CAL_FACT", "records"))
        except:
            self.logger.warning('Cannot get calibrationConst')
            gotCalibration = False

        # add dictionary for raw counts
        metaDict = [{'src': {'SourceFilename': fileName, 'SourceBand': 1},
                     'dst': {'name': 'RawCounts_%s' % polarization}}]

        if gotCalibration:
            # create one pixel band with calibrationConst value and add it to dictionary
            array = np.ones((2,2)) * calibrationConst
            self.calibrationVRT = VRT(array=array)
            metaDict.append({'src': {'SourceFilename': self.calibrationVRT.fileName,
                                     'SourceBand': 1},
                             'dst': {'name': 'calibrationConst'}
                            })

            # add dictionary for IncidenceAngle (and other ADS variables)
            for adsVRT in self.adsVRTs:
                metaDict.append({'src': {'SourceFilename': adsVRT.fileName,
                                         'SourceBand': 1},
                                 'dst': {'name':  adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name'),
                                         'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}
                                })

            # add dicrtionary for sigma0, ice and water
            names = ['sigma0', 'ice', 'water']
            wkt = ['surface_backwards_scattering_coefficient_of_radar_wave',
                    'ice','water']
            sphPass = [gdalMetadata['SPH_PASS'], '', '']

            if caliblation == 'one-line':
                adsBandNo = 0
                sourceFileNames = [self.adsVRTs[adsBandNo].fileName,
                                   self.calibrationVRT.fileName,
                                   fileName]
                pixelFunctionTypes = ['RawcountsIncidenceToSigma0_FromLine',
                                      'Sigma0NormalizedIce_FromLine']
                if polarization=='HH':
                    pixelFunctionTypes.append('Sigma0HHNormalizedWater_FromLine')
                elif polarization=='VV':
                    pixelFunctionTypes.append('Sigma0VVNormalizedWater_FromLine')
            else:
                adsBandNo = 1
                sourceFileNames = [fileName,
                                   self.calibrationVRT.fileName,
                                   self.adsVRTs[adsBandNo].fileName,]
                pixelFunctionTypes = ['RawcountsIncidenceToSigma0',
                                      'Sigma0NormalizedIce']
                if polarization=='HH':
                    pixelFunctionTypes.append('Sigma0HHNormalizedWater')
                elif polarization=='VV':
                    pixelFunctionTypes.append('Sigma0VVNormalizedWater')

            for iPixFunc in range(len(pixelFunctionTypes)):
                metaDict.append({'src': [{'SourceFilename': sourceFileNames[0],
                                          'SourceBand': 1},
                                         {'SourceFilename': sourceFileNames[1],
                                          'SourceBand': 1},
                                         {'SourceFilename': sourceFileNames[2],
                                          'SourceBand': 1}],
                                 'dst': {'name':names[iPixFunc],
                                         'wkv': wkt[iPixFunc],
                                         'PixelFunctionType': pixelFunctionTypes[iPixFunc],
                                         'polarization': polarization,
                                         'suffix': polarization,
                                         'pass': sphPass[iPixFunc],
                                         'dataType': 6}})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays
        #self.add_geolocation_from_ads(gdalDataset, step=1)
