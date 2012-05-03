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
from vrt import *

class Mapper(VRT):
    ''' VRT with mapping of WKV for ASAR Level 1 '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        
        product = gdalMetadata.get("MPH_PRODUCT", "Not_ASAR")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER");

        metaDict = [{'source': fileName,
                     'sourceBand': 1,
                     'wkv': 'normalized_radar_cross_section',
                     'parameters':{
                         'band_name': 'sigma0',
                         'polarisation': gdalMetadata['SPH_MDS1_TX_RX_POLAR'],
                         'pass': gdalMetadata['SPH_PASS'],
                         'warning': 'fake sigma0, not yet calibrated'}}];

            
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        productTime = gdalMetadata["SPH_FIRST_LINE_TIME"]
        self._set_time(dateutil.parser.parse(productTime))