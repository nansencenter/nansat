# Name:        mapper_masr2_l3
# Purpose:     Mapping for L3 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import datetime
import os.path
import glob
import gdal
import numpy as np
from nansat.vrt import VRT, GeolocationArray, NSR, parse


class Mapper(VRT):
    ''' Mapper for Level-3 AMSR2 data from https://gcom-w1.jaxa.jp'''

    # detect wkv from metadata 'Parameter'
    param2wkv = {'Chlorophyll a concentration': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                 'Diffuse attenuation coefficient': 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water',
                 'Remote sensing reflectance': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                 'CDOM Index': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter',
                 'Sea Surface Salinity': 'sea_surface_salinity',
                 'Sea Surface Temperature': 'sea_surface_temperature',
                 'Instantaneous Photosynthetically Available Radiation': 'instantaneous_photosynthetically_available_radiation',
                 'Particle backscatter at 443 nm': 'volume_backscattering_coefficient_of_radiative_flux_in_sea_water_due_to_suspended_particles',
                 'Chlorophyll a concentration, Garver-Siegel-Maritorena Model': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                 'Photosynthetically Available Radiation': 'photosynthetically_available_radiation',
                 }
    freqs = [6, 7, 10, 18, 23, 36, 89]

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' OBPG L3 VRT '''

        # test the product
        assert gdalMetadata['PlatformShortName'] == 'GCOM-W1'
        assert gdalMetadata['SensorShortName'] == 'AMSR2'
        assert gdalMetadata['ProductName'] == 'AMSR2-L3'
        

        # get list of similar (same date, A/D orbit) files in the directory
        iDir, iFile = os.path.split(fileName)
        iFileMask = iFile[:30] + '%02d' + iFile[32:]
        simFiles = []
        for freq in self.freqs:
            simFile = os.path.join(iDir, iFileMask % freq)
            #print simFile
            if os.path.exists(simFile):
                simFiles.append(simFile)

        metaDict = []
        for freq in self.freqs:
            simFile = os.path.join(iDir, iFileMask % freq)
            if simFile not in simFiles:
                continue
            #print 'simFile', simFile
            # open file, get metadata and get parameter name
            simSupDataset = gdal.Open(simFile)
            if simSupDataset is None:
                # skip this similar file
                #print 'No dataset: %s not a supported SMI file' % simFile
                continue
            # get subdatasets from the similar file
            simSubDatasets = simSupDataset.GetSubDatasets()
            for simSubDataset in simSubDatasets:
                #print 'simSubDataset', simSubDataset
                if 'Brightness_Temperature' in simSubDataset[0]:
                    # get SourceFilename from subdataset
                    metaEntry = {'src': {'SourceFilename': simSubDataset[0],
                                         'SourceBand'  : 1,
                                         'ScaleRatio'  : 0.0099999998,
                                         'ScaleOffset' : 0},
                                 'dst': {'wkv'          : 'brightness_temperature',
                                         'frequency'    : '%02d' % freq,
                                         'polarisation' : simSubDataset[0][-2:-1],
                                         'suffix'       : '%02d%s' % (freq, simSubDataset[0][-2:-1])}}
                    metaDict.append(metaEntry)

        # initiate VRT for the NSIDC 10 km grid
        VRT.__init__(self,
                     srcGeoTransform=(-3850000, 10000, 0.0,
                                      5850000, 0.0, -10000),
                     srcProjection=NSR(3411).wkt,
                     srcRasterXSize=760,
                     srcRasterYSize=1120)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # Add valid time
        self._set_time(parse(gdalMetadata['ObservationStartDateTime']))
