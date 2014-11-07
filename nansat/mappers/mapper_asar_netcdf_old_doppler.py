#------------------------------------------------------------------------------
# Name:		mapper_asar_netcdf_old_doppler.py
# Purpose:      To read previously processed ASAR WS Doppler data saved as
#               netcdf files
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	09.10.2014
# Last modified:13.10.2014 17:09
# Copyright:    (c) NERSC
# License:
#------------------------------------------------------------------------------
import warnings

import os
import glob
import numpy as np
import scipy
from dateutil.parser import parse

from nansat.vrt import VRT
from nansat.tools import gdal, WrongMapperError, initial_bearing
from nansat.nsr import NSR
from nansat.node import Node


class Mapper(VRT):
    '''
        Create VRT with mapping of ASAR wide swath Doppled data
    '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):

        # Check this is ASAR old doppler netcdf
        if not len(filename.split('.')) == 3:
            raise WrongMapperError
        if not (filename.split('.')[2] == 'nc' and
                filename.split('.')[1] == 'doppler' and
                filename.split('.')[0][0:3] == 'ASA'):
            raise WrongMapperError

        # Remove 'NC_GLOBAL#' and 'GDAL_' and
        # 'NANSAT_' from keys in gdalDataset
        tmpGdalMetadata = {}
        for key in gdalMetadata.keys():
            newKey = key.replace('NC_GLOBAL#', '')
            tmpGdalMetadata[newKey] = gdalMetadata[key]
        gdalMetadata = tmpGdalMetadata

        # Get file names from subdatasets
        subDatasets = gdalDataset.GetSubDatasets()
        filenames = [f[0] for f in subDatasets]

        lon = [(gdal.Open(filenames.pop(ii)).ReadAsArray()
                for ii, fn in enumerate(filenames)
                if 'lon' in fn)][0]
        lat = [(gdal.Open(filenames.pop(ii)).ReadAsArray()
                for ii, fn in enumerate(filenames)
                if 'lat' in fn)][0]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, lon=lon, lat=lat)

        # Add list of calibration files to global metadata
        #self.dataset.SetMetadataItem(
        #        'Orbit based range bias calibration files',
        #        [filenames.pop(ii) for ii,fn in enumerate(filenames) if
        #            'calibration_file_orbit' in fn][0])
        remove_calfile_info = [(filenames.pop(ii)
                                for ii, fn in enumerate(filenames)
                                if 'calibration_file_orbit' in fn)][0]

        name2wkv_dict = {'azimuth': 'platform_azimuth_angle',
                         'incidence_angles': 'angle_of_incidence',
                         'sat2target_elevation': 'sensor_view_angle',
                         'slant_range_time': '',
                         'dop_coef_observed': '',
                         'dop_coef_predicted': '',
                         'valid': '',
                         'raw_counts': '',
                         'azivar_raw_counts': '',
                         'azibias': '',
                         'range_bias_orbit': '',
                         'range_bias_std_orbit': '',
                         'valid_orbit': '',
                         'range_bias_scene': '',
                         'range_bias_std_scene': '',
                         'valid_scene': '',
                         }
        metaDict = []
        bandNo = {}
        for i, filename in enumerate(filenames):
            subDataset = gdal.Open(filename)
            subBand = subDataset.GetRasterBand(1)
            # generate src metadata
            src = {'SourceFilename': filename, 'SourceBand': 1}
            src['DataType'] = subBand.DataType

            bandMetadata = subBand.GetMetadata_Dict()
            bandNo[bandMetadata.get('NETCDF_VARNAME')] = i + 1
            # generate dst metadata
            dst = {'wkv': name2wkv_dict[bandMetadata.get('NETCDF_VARNAME')],
                   'name': bandMetadata.get('NETCDF_VARNAME'),
                   }
            metaDict.append({'src': src, 'dst': dst})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        metaDict = []
        dco = (self.dataset.GetRasterBand(bandNo['dop_coef_observed']).
               ReadAsArray())
        dcp = (self.dataset.GetRasterBand(bandNo['dop_coef_predicted']).
               ReadAsArray())
        azibias = self.dataset.GetRasterBand(bandNo['azibias']).ReadAsArray()
        range_bias = (self.dataset.GetRasterBand(bandNo['range_bias_scene']).
                      ReadAsArray())
        fdg = dco - dcp - azibias - range_bias
        fdg[range_bias > 10000] = np.nan
        fdg[azibias > 10000] = np.nan
        fdgVRT = VRT(array=fdg, lat=lat, lon=lon)

        mask = np.ones(np.shape(dco))
        mask[range_bias > 10000] = 0.
        mask[azibias > 10000] = 0.
        maskVRT = VRT(array=mask, lat=lat, lon=lon)

        self.subVRTs['fdgVRT'] = fdgVRT
        metaDict.append({'src': {'SourceFilename': (
                                 self.subVRTs['fdgVRT'].fileName),
                                 'SourceBand': 1
                                 },
                         'dst': {'name': 'fdg',
                                 'long_name': (
                                 'Line of sight geophysical Doppler shift'),
                                 'units': 'Hz',
                                 }
                         })
        self.subVRTs['maskVRT'] = maskVRT
        metaDict.append({'src': {'SourceFilename': (
                                 self.subVRTs['maskVRT'].fileName),
                                 'SourceBand': 1
                                 },
                        'dst': {'name': 'mask',
                                'long_name': 'Mask for use in plotting',
                                }
                         })

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        self.dataset.SetMetadataItem('mapper', 'asar_netcdf_old_doppler')
