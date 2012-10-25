#-------------------------------------------------------------------------------
# Name:        mapper_kmssL1
# Purpose:     Mapping for KMSS-L1 data
#
# Author:      evgenym
#
# Created:     01.08.2012
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#
#%assign values to L_500, L_655, L_800
#L_555 = img(:, :, 1);
#L_655 = img(:, :, 2);
#L_800 = img(:, :, 3);

#
#from vrt import *
#import os.path


from datetime import datetime
from numpy import mod

try:
    from osgeo import gdal
except ImportError:
    import gdal

from vrt import *
from domain import Domain


class Mapper(VRT):
    ''' VRT with mapping of WKV for KMSS TOA tiff data'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=10):
        ''' Create VRT '''
       # print 'KMSS_MAPPER_Started'
       # product = gdalMetadata.get("SATELLITE_IDENTIFIER", "Not_KMSS_tiff")
       # print "gdalMetadata_out"
       # print product
        #if it is not KMSS, return
       # if product!= 'KMSS':
       #     raise AttributeError("KMSS BAD MAPPER");
        product = gdalDataset.GetDriver().LongName
        print 'KMSS ' + product
        if product!= 'GeoTIFF':
            raise AttributeError("Not_KMSS_tiff");
      
        metaDict = [
        {'source': fileName, 'sourceBand':  1, 'wkv': 'toa_outgoing_spectral_radiance', 'parameters': {'wavelength': '555'}},
        {'source': fileName, 'sourceBand':  2, 'wkv': 'toa_outgoing_spectral_radiance', 'parameters': {'wavelength': '655'}},
        {'source': fileName, 'sourceBand':  3, 'wkv': 'toa_outgoing_spectral_radiance', 'parameters': {'wavelength': '800'}},
        ];
        # from https://gsics.nesdis.noaa.gov/wiki/Development/StandardVariableNames
            
        # add 'band_name' to 'parameters'
        #for bandDict in metaDict:
        #    if bandDict['parameters'].has_key('wavelength'):
        #        bandDict['parameters']['band_name'] = 'radiance_' + bandDict['parameters']['wavelength']    
        
        # set scale factor to the band metadata (only radiances)
        #for i, bandDict in enumerate(metaDict[:-1]):
         #   bandDict['parameters']['scale'] = str(scales[i])
        
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._add_all_bands(metaDict)


       
