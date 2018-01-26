# Name:        mapper_kmssL1
# Purpose:     Mapping for KMSS-L1 data
# Author:      Evgeny Morozov
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html
import os
from datetime import datetime

from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT
from nansat.domain import Domain


class Mapper(VRT):
    ''' VRT with mapping of WKV for KMSS TOA tiff data'''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' Create VRT '''
        if (os.path.split(filename)[1][0:4] != '101_' or
                os.path.split(filename)[1][0:4] != '102_'):
                raise WrongMapperError

        try:
            product = gdalDataset.GetDriver().LongName
        except:
            raise WrongMapperError

        if (product != 'GeoTIFF' or filename[-3:] != 'tif' or
                gdalDataset.RasterCount != 3):
            raise WrongMapperError

        metaDict = [{'src': {'SourceFilename': filename, 'SourceBand': 1},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '555'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 2},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '655'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 3},
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
        self._init_from_gdal_dataset(gdalDataset)

         # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)
