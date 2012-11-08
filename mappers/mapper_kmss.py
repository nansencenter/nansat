#-------------------------------------------------------------------------------
# Name:        mapper_kmssL1
# Purpose:     Mapping for KMSS-L1 data
#
# Author:      evgenym(me)
#
# Created:     01.08.2012
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

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
        product = gdalDataset.GetDriver().LongName

        raise AttributeError("Not_KMSS_tiff");
        if product!= 'GeoTIFF':
            raise AttributeError("Not_KMSS_tiff");

        metaDict = [
        {'src': {'SourceFilename': fileName, 'SourceBand':  1},
         'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '555'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  2},
         'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '655'}},
        {'src': {'SourceFilename': fileName, 'SourceBand':  3},
         'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '800'}},
        ]
	# from https://gsics.nesdis.noaa.gov/wiki/Development/StandardVariableNames

        # add DataType into 'src' and name into 'dst'
        for bandDict in metaDict:
            if 'DataType' not in bandDict['src']:
                bandDict['src']['DataType'] = 2
            if bandDict['dst'].has_key('wavelength'):
                bandDict['dst']['name'] = 'toa_radiance_' + bandDict['dst']['wavelength']

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset);

         # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
