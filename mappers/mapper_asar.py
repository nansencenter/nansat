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

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        '''
        Parameters (**kwargs)
        ---------------------
        ASA_full_incAng : bool (default False)
            if True, use full-size incidence angle band.
            if False, use one-line incidence angle band.
        '''
        product = gdalMetadata.get("MPH_PRODUCT")
        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        kwDict = {'geolocation' : False}
        # choose kwargs for envisat and asar and change keyname
        for key in kwargs:
            if key.startswith('envisat') or key.startswith('asar'):
                keyName = key.replace('envisat_', '').replace('asar_', '')
                kwDict[keyName] = kwargs[key]
            else:
                kwDict[key] = kwargs[key]

        Envisat.__init__(self, fileName, product[0:4], **kwDict)
        # get polarization string (remove '/', since NetCDF doesnt support that in metadata)
        polarization = gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", "")

        # Create VRTdataset with small VRTRawRasterbands
        self.adsVRTs = self.get_ads_vrts(gdalDataset,
                                         ["first_line_incidenceAngle"])

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, **kwDict)

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
            # add dicrtionary for sigma0, ice and water
            names = ['sigma0', 'ice', 'water']
            wkt = ['surface_backwards_scattering_coefficient_of_radar_wave',
                    'ice', 'water']
            sphPass = [gdalMetadata['SPH_PASS'], '', '']

            sourceFileNames = [fileName,
                               self.adsVRTs[0].fileName]

            pixelFunctionTypes = ['RawcountsIncidenceToSigma0',
                                  'Sigma0NormalizedIce']
            if polarization=='HH':
                pixelFunctionTypes.append('Sigma0HHNormalizedWater')
            elif polarization=='VV':
                pixelFunctionTypes.append('Sigma0VVNormalizedWater')

            # add pixelfunction bands to metaDict
            for iPixFunc in range(len(pixelFunctionTypes)):
                srcFiles = []
                for iFileName in sourceFileNames:
                    sourceFile = {'SourceFilename': iFileName,
                                  'SourceBand': 1}
                    # if ASA_full_incAng, set 'ScaleRatio' into source file dict
                    if iFileName == fileName:
                        sourceFile['ScaleRatio'] = np.sqrt(1.0/calibrationConst)
                    srcFiles.append(sourceFile)

                metaDict.append({'src': srcFiles,
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
        if self.d['geolocation']:
            self.add_geolocation_from_ads(gdalDataset)
