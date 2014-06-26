# Name:        mapper_kmssL1
# Purpose:     Mapping for KMSS-L1 data
# Author:      Evgeny Morozov
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html
from datetime import datetime

from osgeo import gdal

from nansat.vrt import VRT
from nansat.domain import Domain


class Mapper(VRT):
    ''' VRT with mapping of WKV for KMSS TOA tiff data'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create VRT '''
        product = gdalDataset.GetDriver().LongName
        if cmp(os.path.split(fileName)[1][0:4], '101_') != 0:
            if cmp(os.path.split(fileName)[1][0:4], '102_') != 0:
                raise AttributeError("NTSOMZ GeoTIFF KMSS filename usually starts with '101' or '102'")

        if product != 'GeoTIFF':
            raise AttributeError("Not_GeoTIFF")
        if cmp(fileName[-3:], 'tif') != 0:
            raise AttributeError("for NTSOMZ GeoTIFF extension must be tif")
        if gdalDataset.RasterCount != 3:
            raise AttributeError("Not_NTSOMZ_KMSS_geotiff! does not have 3 bands!")

        metaDict = [{'src': {'SourceFilename': fileName, 'SourceBand': 1},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '555'}},
                    {'src': {'SourceFilename': fileName, 'SourceBand': 2},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '655'}},
                    {'src': {'SourceFilename': fileName, 'SourceBand': 3},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '800'}}
                    ]
        # from https://gsics.nesdis.noaa.gov/wiki/Development/StandardVariableNames

        # add DataType into 'src' and name into 'dst'
        for bandDict in metaDict:
            if 'DataType' not in bandDict['src']:
                bandDict['src']['DataType'] = 2
            if 'wavelength' in bandDict['dst']:
                bandDict['dst']['name'] = ('toa_radiance_' +
                                           bandDict['dst']['wavelength'])

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

         # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
