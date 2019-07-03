# Name:        nansat_mapper_merisL1
# Purpose:     Mapping for Meris-L1 data
# Authors:     Anton Korosov
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse
import json

from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError
from nansat.mappers.envisat import Envisat

import pythesint as pti

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for MERIS Level 1 (FR or RR) '''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 geolocation=False, zoomSize=500, step=1, **kwargs):

        ''' Create MER1 VRT

        Parameters
        -----------
        filename : string

        gdalDataset : gdal dataset

        gdalMetadata : gdal metadata

        geolocation : bool (default is False)
            if True, add gdal geolocation

        zoomSize: int (used in envisat.py)
            size, to which the ADS array will be zoomed using scipy
            array of this size will be stored in memory

        step: int (used in envisat.py)
            step of pixel and line in GeolocationArrays. lat/lon grids are
            generated at that step

        '''

        self.setup_ads_parameters(filename, gdalMetadata)

        if (self.product[0:9] != "MER_FRS_1" and
                self.product[0:9] != "MER_RR__1"):
            raise WrongMapperError

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalDataset)

        metaDict = [{'src': {'SourceFilename': filename, 'SourceBand': 1},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '413'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 2},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '443'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 3},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '490'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 4},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '510'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 5},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '560'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 6},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '620'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 7},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '665'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 8},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '681'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 9},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '709'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 10},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '753'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 11},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '761'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 12},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '778'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 13},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '864'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 14},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '849'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 15},
                     'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                             'wavelength': '900'}},
                    {'src': {'SourceFilename': filename, 'SourceBand': 16,
                             'DataType': 1},
                     'dst': {'wkv': 'quality_flags', 'suffix': 'l1'}}
                    ]

        # add DataType into 'src' and suffix into 'dst'
        for bandDict in metaDict:
            if 'DataType' not in bandDict['src']:
                # default for meris L1 DataType UInt16
                bandDict['src']['DataType'] = 2
            if 'wavelength' in bandDict['dst']:
                bandDict['dst']['suffix'] = bandDict['dst']['wavelength']

        # get scaling GADS from header
        scales = self.read_scaling_gads(range(7, 22))
        # set scale factor to the band metadata (only radiances)
        for i, bandDict in enumerate(metaDict[:-1]):
            bandDict['src']['ScaleRatio'] = str(scales[i])

        # get list with resized VRTs from ADS
        self.band_vrts = {'adsVRTs': self.get_ads_vrts(gdalDataset,
                                                     ['sun zenith angles',
                                                      'sun azimuth angles',
                                                      'zonal winds',
                                                      'meridional winds'],
                                                     zoomSize=zoomSize,
                                                     step=step)}
        # add bands from the ADS VRTs
        for adsVRT in self.band_vrts['adsVRTs']:
            metaDict.append({'src': {'SourceFilename': adsVRT.filename,
                                     'SourceBand': 1},
                             'dst': {'name': (adsVRT.dataset.GetRasterBand(1).
                                              GetMetadataItem('name')),
                                     'units': (adsVRT.dataset.GetRasterBand(1).
                                               GetMetadataItem('units'))}
                             })

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('meris')
        ee = pti.get_gcmd_platform('envisat')

        # TODO: Validate that the found instrument and platform are indeed what we
        # want....

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

        # add geolocation arrays
        if geolocation:
            self.add_geolocation_from_ads(gdalDataset,
                                          zoomSize=zoomSize, step=step)
