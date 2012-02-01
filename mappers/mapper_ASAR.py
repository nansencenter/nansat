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

    def __init__(self, rawVRTFileName, fileName, dataset, metadata, vrtBandList):
        ''' Create ASAR VRT '''
        VRT.__init__(self, dataset, metadata, rawVRTFileName);
        
        product = metadata.get("MPH_PRODUCT", "Not_ASAR")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER");

        metaDict = [{'source': fileName, 'sourceBand': 1, 'wkv': 'sigma0', \
                     'parameters':{'polarisation': metadata['SPH_MDS1_TX_RX_POLAR'], \
                    'pass': metadata['SPH_PASS'],'warning': 'fake sigma0, not yet calibrated'}}];

        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);
            
        self._createVRT(metaDict, vrtBandList);
