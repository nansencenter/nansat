#-------------------------------------------------------------------------------
# Name:        nansat_mapper_asar
# Purpose:     Mapping for Envisat ASAR-L1 data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import VRT
from envisat import Envisat
import numpy as np

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        product = gdalMetadata.get("MPH_PRODUCT")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        # init ADS parameters
        Envisat.__init__(self, fileName, product[0:4])

        # get polarization string (remove '/', since NetCDF doesnt support that in metadata)
        polarization = gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", "")

        # Create VRTdataset with small VRTRawRasterbands
        self.adsVRTs = self.get_ads_vrts(gdalDataset, ["first_line_incidenceAngle"])

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
            # add dictionary for IncidenceAngle (and other ADS variables)
            for adsVRT in self.adsVRTs:
                metaDict.append({'src': {'SourceFilename': adsVRT.fileName,
                                         'SourceBand': 1},
                                 'dst': {'name':  adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name'),
                                         'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}
                                })
    
            # add dictionary for sigma0
            metaDict.append({'src': [{'SourceFilename': fileName,
                                      'SourceBand': 1,
                                      'ScaleRatio': np.sqrt(1.0/calibrationConst)},
                                     {'SourceFilename': self.adsVRTs[0].fileName,
                                      'SourceBand': 1}],
    
                            'dst': { 'wkv':  'surface_backwards_scattering_coefficient_of_radar_wave',
                                     'PixelFunctionType': 'RawcountsIncidenceToSigma0',
                                     'polarization': polarization,
                                     'suffix': polarization,
                                     'pass': gdalMetadata['SPH_PASS'],
                                     'dataType': 6}})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays
        #self.add_geolocation_from_ads(gdalDataset, step=1)
