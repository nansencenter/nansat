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

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1 '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        product = gdalMetadata.get("MPH_PRODUCT", "Not_ASAR")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER");

        metaDict = [{'source': fileName,
                     'sourceBand': 1,
                     'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                     'parameters':{
                         'band_name': 'sigma0',
                         'polarisation': gdalMetadata['SPH_MDS1_TX_RX_POLAR'],
                         'pass': gdalMetadata['SPH_PASS'],
                         'warning': 'fake sigma0, not yet calibrated'}}];


        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)
