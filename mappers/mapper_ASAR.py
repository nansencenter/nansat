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
from vrt import VRT, Geolocation
from envisat import Envisat
import numpy as np
import os
import gdal

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        product = gdalMetadata.get("MPH_PRODUCT", "Not_ASAR")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        # Create VRTdataset with small VRTRawRasterbands
        incAngleDataset = self.create_VRT_from_ADSRarray(fileName, "incidenceAngle")

        # Enlarge the band to the underlying data band size
        self.incAngleDataset = incAngleDataset.resized(gdalDataset.RasterXSize, gdalDataset.RasterYSize)

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # get calibration constant
        calibrationConst = float(gdalDataset.GetMetadataItem("MAIN_PROCESSING_PARAMS_ADS_1_CALIBRATION_FACTORS.1.EXT_CAL_FACT", "records"))

        # compute calibrated sigma0 (netcdf does not accept "/")
        metaDict = [{'src': {'SourceFilename': fileName, 'SourceBand': 1},
                     'dst': {'name': 'RawCounts_%s' % gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", "")}}]

        # add dictionary for IncidenceAngle
        for iBand in range(self.incAngleDataset.dataset.RasterCount):
            bandMetadata = self.incAngleDataset.dataset.GetRasterBand(iBand+1).GetMetadata()
            metaDict.append({'src': {'SourceFilename': self.incAngleDataset.fileName, 'SourceBand': iBand+1},
                             'dst': bandMetadata})
        # add dictionary for real sigma0
        metaDict.append({'src': [{'SourceFilename': fileName, 'SourceBand': 1,
                                    'ScaleRatio': np.sqrt(1.0/calibrationConst)},
                                 {'SourceFilename': self.incAngleDataset.fileName, 'SourceBand': 1}],
                         'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                                 'PixelFunctionType': 'RawcountsIncidenceToSigma0',
                                 'polarisation': gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", ""),
                                 'name': 'sigma0_%s' % gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", ""),
                                 'pass': gdalMetadata['SPH_PASS'],
                                 'dataType': gdal.GDT_Float32}})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        ''' Set GeolocationArray '''
        latlonName = {"latitude":"lats","longitude":"longs"}
        self.add_geoarray_dataset(fileName, product[0:4], gdalDataset.RasterXSize, gdalDataset.RasterYSize, latlonName, gdalDataset.GetGCPProjection(), ["num_lines"])


