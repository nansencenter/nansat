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

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        ''' Create SeaWifs_L2 VRT '''
        if gdalMetadata['Sensor Name'] != 'SeaWiFS':
            raise AttributeError("SeaWiFS BAD MAPPER");
        
        subDsString = 'HDF4_SDS:SEAWIFS_L2:"%s":%s'
        
        # Trying to read band 21 = chlor_a The parameter height is added just to check that it does not crash due to no parameters
        metaDict = [
        {'source': subDsString % (fileName, '21'), 'sourceBand': 1, 'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'parameters': {'band_name': 'chlor_a'}}
        ]
                
        #open subdataset
        gdalSubDataset = gdal.Open(gdalDataset.GetSubDatasets()[0][0]);    

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalSubDataset, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        return