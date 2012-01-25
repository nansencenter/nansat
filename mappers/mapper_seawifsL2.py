#-------------------------------------------------------------------------------
# Name:        mapper_modisL1
# Purpose:     Mapping for SeaWifs-L2 data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

class Mapper(VRT):
    ''' VRT with mapping of WKV for SeaWifs Level2        '''

    def __init__(self, ds, fileName, metadata, vrtBandList, rawVRTName):
        ''' Create SeaWifs_L2 VRT '''
        VRT.__init__(self, metadata, rawVRTName);
        
        if metadata['Sensor Name'] != 'SeaWiFS':
            raise AttributeError("SeaWiFS BAD MAPPER");
        
        subDsString = 'HDF4_SDS:SEAWIFS_L2:"%s":%s'
        
        # Trying to read band 21 = chlor_a The parameter height is added just to check that it does not crash due to no parameters
        metaDict = [\
        {'source': subDsString % (fileName, '21'), 'sourceBand': 1, 'wkv': 'chlor_a'}\
        ];
                
        #set number of bands
        if vrtBandList == None:
            vrtBandList = range(1,len(metaDict)+1);

        #open subdataset
        subDs = gdal.Open(ds.GetSubDatasets()[0][0]);    

        self.createVRT_and_add_bands(subDs, metaDict, vrtBandList);

        return
