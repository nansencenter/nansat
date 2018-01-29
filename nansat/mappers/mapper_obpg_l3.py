# Name:        mapper_obpg_l3
# Purpose:     Mapping for L3 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from __future__ import unicode_literals, absolute_import, division, print_function
import datetime
import os.path
import glob

import numpy as np

from nansat.tools import gdal, ogr
from nansat.vrt import VRT
from nansat.nsr import NSR

from nansat.exceptions import WrongMapperError

class Mapper(VRT):
    ''' Mapper for Level-3 Standard Mapped Image from
    http://oceancolor.gsfc.nasa.gov'''

    # detect wkv from metadata 'Parameter'
    param2wkv = {'Chlorophyll a concentration': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                 'Diffuse attenuation coefficient': 'volume_attenuation_coefficient_of_downwelling_'
                                                    'radiative_flux_in_sea_water',
                 'Remote sensing reflectance': 'surface_ratio_of_upwelling_radiance_emerging_from_'
                                               'sea_water_to_downwelling_radiative_flux_in_air',
                 'CDOM Index': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_'
                               'to_dissolved_organic_matter',
                 'Sea Surface Salinity': 'sea_surface_salinity',
                 'Sea Surface Temperature': 'sea_surface_temperature',
                 'Instantaneous Photosynthetically Available Radiation': 'instantaneous_photosynthetically_available_radiation',
                 'Particle backscatter at 443 nm': 'volume_backscattering_coefficient_of_radiative_flux_in_sea_water_due_to_suspended_particles',
                 'Chlorophyll a concentration, Garver-Siegel-Maritorena Model': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                 'Photosynthetically Available Radiation': 'downwelling_photosynthetic_photon_radiance_in_sea_water',
                 'Instantaneous Photosynthetically Available Radiation': 'instantaneous_downwelling_photosynthetic_photon_radiance_in_sea_water',
                 }

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' OBPG L3 VRT '''

        try:
            assert 'Level-3 Standard Mapped Image' in gdalMetadata['Title']
        except:
            raise WrongMapperError

        # get list of similar (same date) files in the directory
        iDir, iFile = os.path.split(filename)
        iFileName, iFileExt = os.path.splitext(iFile)
        simFilesMask = os.path.join(iDir, iFileName)
        simFiles = glob.glob(simFilesMask + iFileExt[0:6] + '*')
        #print 'simFilesMask, simFiles', simFilesMask, simFiles

        metaDict = []
        for simFile in simFiles:
            #print 'simFile', simFile
            # open file, get metadata and get parameter name
            simSupDataset = gdal.Open(simFile)
            if simSupDataset is None:
                # skip this similar file
                #print 'No dataset: %s not a supported SMI file' % simFile
                continue
            # get subdatasets from the similar file
            simSubDatasets = simSupDataset.GetSubDatasets()
            if len(simSubDatasets) > 0:
                for simSubDataset in simSubDatasets:
                    #print 'simSubDataset', simSubDataset
                    if 'l3m_data' in simSubDataset[1]:
                        # get SourceFilename from subdataset
                        tmpSourceFilename = simSubDataset[0]
                        break
            else:
                # get SourceFilename from dataset
                tmpSourceFilename = simFile

            # open subdataset with GDAL
            #print 'tmpSourceFilename', tmpSourceFilename
            tmpGdalDataset = gdal.Open(tmpSourceFilename)

            try:
                # get metadata, get 'Parameter'
                tmpGdalMetadata = tmpGdalDataset.GetMetadata()
                simParameter = tmpGdalMetadata['Parameter']
            except:
                print('No parameter: %s not a supported SMI file')
                continue
            else:
                # set params of the similar file
                simSourceFilename = tmpSourceFilename
                simGdalDataset = tmpGdalDataset
                simGdalMetadata = tmpGdalMetadata

            # get WKV from the similar file
            #print 'simParameter', simParameter
            for param in self.param2wkv:
                #print 'param', param
                if param in simParameter:
                    simWKV = self.param2wkv[param]
                    break

            # generate entry to metaDict
            metaEntry = {'src': {'SourceFilename': simSourceFilename,
                                 'SourceBand':  1,
                                 'ScaleRatio': float(simGdalMetadata['Slope']),
                                 'ScaleOffset': float(simGdalMetadata['Intercept'])},
                         'dst': {'wkv': simWKV}}

            # add wavelength and BandName
            if ' at ' in simParameter and ' nm' in simParameter:
                simWavelength = simParameter.split(' at ')[1].split(' nm')[0]
                metaEntry['dst']['suffix'] = simWavelength
                metaEntry['dst']['wavelength'] = simWavelength

            # add band with Rrsw
            metaEntry2 = None
            if simWKV == 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air':
                metaEntry2 = {'src': [metaEntry['src']]}
                metaEntry2['dst'] = {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_water',
                                     'suffix': simWavelength,
                                     'wavelength': simWavelength,
                                     'PixelFunctionType': 'NormReflectanceToRemSensReflectance',
                                     }

            # append entry to metaDict
            metaDict.append(metaEntry)
            if metaEntry2 is not None:
                metaDict.append(metaEntry2)

        #get array with data and make 'mask'
        a = simGdalDataset.ReadAsArray()
        mask = np.zeros(a.shape, 'uint8') + 64
        mask[a < -32000] = 1
        self.band_vrts = {'mask': VRT(array=mask)}

        metaDict.append(
            {'src': {'SourceFilename': self.band_vrts['mask'].filename,
                     'SourceBand':  1},
             'dst': {'name': 'mask'}})

        # create empty VRT dataset with geolocation only
        # print 'simGdalMetadata', simGdalMetadata
        latitudeStep = float(simGdalMetadata.
                             get('Latitude Step',
                                 simGdalMetadata.get('Latitude_Step', 1)))
        longitudeStep = float(simGdalMetadata.
                              get('Longitude Step',
                                  simGdalMetadata.get('Longitude_Step', 1)))
        numberOfColumns = int(simGdalMetadata.
                              get('Number of Columns',
                                  simGdalMetadata.get('Number_of_Columns', 1)))
        numberOfLines = int(simGdalMetadata.
                            get('Number of Lines',
                                simGdalMetadata.get('Number_of_Lines', 1)))
        #longitudeStep = float(simGdalMetadata['Longitude Step'])
        # x_size, y_size, geo_transform, projection, gcps=None, gcp_projection='', **kwargs
        self._init_from_dataset_params(numberOfColumns, numberOfLines,
                                (-180.0, longitudeStep, 0.0, 90.0, 0.0, -longitudeStep),
                                NSR().wkt)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # Add valid time
        startYear = int(simGdalMetadata.get('Start Year',
                                            simGdalMetadata.
                                            get('Start_Year', 1)))
        startDay = int(simGdalMetadata.get('Start Day',
                                           simGdalMetadata.
                                           get('Start)Day', 1)))
        self.dataset.SetMetadataItem('time_coverage_start',
            (datetime.datetime(startYear, 1, 1) +
                       datetime.timedelta(startDay)).isoformat())
