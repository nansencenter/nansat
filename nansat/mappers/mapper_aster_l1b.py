# Name:        mapper_modisL1
# Purpose:     Mapping for MODIS-L1 data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse
import warnings
import json

import pythesint as pti

from nansat.utils import gdal, ogr
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT
from hdf4_mapper import HDF4Mapper


class Mapper(HDF4Mapper):
    ''' VRT with mapping of WKV for MODIS Level 1 (QKM, HKM, 1KM) '''

    def __init__(self, filename, gdalDataset, gdalMetadata, emrange='VNIR', **kwargs):
        ''' Create MODIS_L1 VRT '''
        # check mapper
        try:
            INSTRUMENTSHORTNAME = gdalMetadata['INSTRUMENTSHORTNAME']
        except:
            raise WrongMapperError
        if INSTRUMENTSHORTNAME != 'ASTER':
            raise WrongMapperError
        try:
            SHORTNAME = gdalMetadata['SHORTNAME']
        except:
            raise WrongMapperError
        if SHORTNAME != 'ASTL1B':
            raise WrongMapperError

        # set up metadict for data with various resolution
        subDSString = 'HDF4_EOS:EOS_SWATH:"%s":%s:%s'
        metaDictVNIR = [
        {'src': {'SourceFilename': subDSString % (filename, 'VNIR_Swath', 'ImageData1' )}, 'dst': {'wavelength': '560'}},
        {'src': {'SourceFilename': subDSString % (filename, 'VNIR_Swath', 'ImageData2' )}, 'dst': {'wavelength': '660'}},
        {'src': {'SourceFilename': subDSString % (filename, 'VNIR_Swath', 'ImageData3N')}, 'dst': {'wavelength': '820'}},
        {'src': {'SourceFilename': subDSString % (filename, 'VNIR_Swath', 'ImageData3B')}, 'dst': {'wavelength': '820'}},
        ]

        metaDictSWIR = [
        {'src': {'SourceFilename': subDSString % (filename, 'SWIR_Swath', 'ImageData4')}, 'dst': {'wavelength': '1650'}},
        {'src': {'SourceFilename': subDSString % (filename, 'SWIR_Swath', 'ImageData5')}, 'dst': {'wavelength': '2165'}},
        {'src': {'SourceFilename': subDSString % (filename, 'SWIR_Swath', 'ImageData6')}, 'dst': {'wavelength': '2205'}},
        {'src': {'SourceFilename': subDSString % (filename, 'SWIR_Swath', 'ImageData7')}, 'dst': {'wavelength': '2260'}},
        {'src': {'SourceFilename': subDSString % (filename, 'SWIR_Swath', 'ImageData8')}, 'dst': {'wavelength': '2330'}},
        {'src': {'SourceFilename': subDSString % (filename, 'SWIR_Swath', 'ImageData9')}, 'dst': {'wavelength': '2395'}},
        ]

        metaDictTIR = [
        {'src': {'SourceFilename': subDSString % (filename, 'TIR_Swath',  'ImageData10')}, 'dst': {'wavelength': '8300'}},
        {'src': {'SourceFilename': subDSString % (filename, 'TIR_Swath',  'ImageData11')}, 'dst': {'wavelength': '8650'}},
        {'src': {'SourceFilename': subDSString % (filename, 'TIR_Swath',  'ImageData12')}, 'dst': {'wavelength': '9100'}},
        {'src': {'SourceFilename': subDSString % (filename, 'TIR_Swath',  'ImageData13')}, 'dst': {'wavelength': '10600'}},
        {'src': {'SourceFilename': subDSString % (filename, 'TIR_Swath',  'ImageData14')}, 'dst': {'wavelength': '11300'}},
        ]

        # select appropriate metaDict based on <emrange> parameter
        metaDict = {'VNIR': metaDictVNIR,
                    'SWIR': metaDictSWIR,
                    'TIR': metaDictTIR,
                    }[emrange]

        # get 1st EOS subdataset and parse to VRT.__init__()
        # for retrieving geo-metadata
        try:
            gdalSubDataset0 = gdal.Open(metaDict[0]['src']['SourceFilename'])
        except (AttributeError, IndexError):
            raise WrongMapperError

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalSubDataset0, metadata=gdalMetadata)

        # add source band, wkv and suffix
        for metaEntry in metaDict:
            metaEntry['src']['SourceBand'] = 1
            metaEntry['dst']['wkv'] = 'toa_outgoing_spectral_radiance'
            metaEntry['dst']['suffix'] = metaEntry['dst']['wavelength']

            if 'ImageData3N' in metaEntry['src']['SourceFilename']:
                metaEntry['dst']['suffix'] += 'N'

            if 'ImageData3B' in metaEntry['src']['SourceFilename']:
                metaEntry['dst']['suffix'] += 'B'

        # add scale and offset
        for metaEntry in metaDict:
            bandNo = metaEntry['src']['SourceFilename'].strip().split(':')[-1].replace('ImageData', '')
            metaEntry['src']['ScaleRatio'] = float(gdalMetadata['INCL' + bandNo])
            metaEntry['src']['ScaleOffset'] = float(gdalMetadata['OFFSET' + bandNo])

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # set time
        datetimeString = self.find_metadata(gdalMetadata, "SETTINGTIMEOFPOINTING")
        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start',
                                     parse(datetimeString+'+00').isoformat())
        self.dataset.SetMetadataItem('time_coverage_end',
                                     parse(datetimeString+'+00').isoformat())

        mm = pti.get_gcmd_instrument('ASTER')
        ee = pti.get_gcmd_platform('TERRA')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

        self._remove_geolocation()
