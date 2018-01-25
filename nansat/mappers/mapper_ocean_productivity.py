# Name:        mapper_ocean_prodcutivity
# Purpose:     Mapping for MODIS and SeaWiFS Level-3
#              from Ocean Productivity website (Oregon State University)
# Authors:      Dmitry Petrenko, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import datetime
import os.path
import glob

import numpy as np

from nansat.tools import gdal, ogr
from nansat.vrt import VRT

from nansat.exceptions import WrongMapperError

class Mapper(VRT):
    ''' Mapper for Ocean Productivity website
    http://www.science.oregonstate.edu/ocean.productivity/'''
    # detect wkv from metadata 'Parameter'
    param2wkv = {'chl': 'mass_concentration_of_chlorophyll_a_in_sea_water',
                 'sst': 'sea_surface_temperature',
                 'par': 'downwelling_photosynthetic_photon_radiance_in_sea_water',
                 'ipar': 'instantaneous_downwelling_photosynthetic_photon_radiance_in_sea_water',
                 'bbp': 'particle_backscatter_at_443_nm'
                 }

    bandNames = {'mass_concentration_of_chlorophyll_a_in_sea_water': 'algal_1',
                 'sea_surface_temperature': 'SST',
                 'instantaneous_downwelling_photosynthetic_photon_radiance_in_sea_water': 'ipar',
                 'particle_backscatter_at_443_nm': 'bbp_443'
                 }

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' Ocean Productivity website VRT '''

        try:
            assert 'IDL' in gdalMetadata['Projection Category']
            assert '-9999' in gdalMetadata['Hole Value']
        except:
            raise WrongMapperError

        print('Ocean Productivity website data')
        # get list of similar (same date) files in the directory
        iDir, iFile = os.path.split(filename)
        iFileName, iFileExt = os.path.splitext(iFile)
        simFilesMask = os.path.join(iDir, '*' + iFileName[4:11] + iFileExt)
        #print 'simFilesMask', simFilesMask
        simFiles = glob.glob(simFilesMask)
        #print 'simFiles', simFiles

        metaDict = []
        for simFile in simFiles:
            #print 'simFile',simFile
            # open subdataset with GDAL
            tmpSourceFilename = simFile
            tmpGdalDataset = gdal.Open(tmpSourceFilename)

            # get metadata, get 'Parameter'
            tmpGdalMetadata = tmpGdalDataset.GetMetadata()
            iDir, ifilename = os.path.split(tmpSourceFilename)
            #print 'ifilename',ifilename
            simParameter = ifilename[0:3]

            # set params of the similar file
            simSourceFilename = tmpSourceFilename
            simGdalDataset = tmpGdalDataset
            simGdalMetadata = tmpGdalMetadata

            # get WKV from the similar file
            for param in self.param2wkv:
                #print 'param', param
                if param in simParameter:
                    simWKV = self.param2wkv[param]
                    break
            #print 'simWKV', simWKV
            # generate entry to metaDict
            metaEntry = {
                'src': {'SourceFilename': simSourceFilename,
                        'SourceBand':  1,
                        'ScaleRatio': float(simGdalMetadata['Slope']),
                        'ScaleOffset': float(simGdalMetadata['Intercept'])},
                'dst': {'wkv': simWKV,
                        'name': self.bandNames[simWKV],
                        'Parameter': simParameter}}
            #print 'metaEntry', metaEntry
            # append entry to metaDict
            metaDict.append(metaEntry)

        #get array with data and make 'mask'
        a = simGdalDataset.ReadAsArray()
        mask = np.zeros(a.shape, 'uint8') + 128
        mask[a < -9990] = 1
        self.band_vrts = {'maskVRT': VRT(array=mask)}

        metaDict.append({'src': {'SourceFilename': (self.band_vrts['maskVRT'].
                                                    filename),
                                 'SourceBand': 1},
                         'dst': {'name': 'mask'}})

        # create empty VRT dataset with geolocation only
        # print 'simGdalMetadata', simGdalMetadata
        latitudeStep = 0.08333334
        longitudeStep = 0.08333334
        numberOfColumns = 4320
        numberOfLines = 2160
        #longitudeStep = float(simGdalMetadata['Longitude Step'])
        self._init_from_dataset_params(numberOfColumns, numberOfLines,
                                       (-180.0, longitudeStep, 0.0, 90.0, 0.0, -longitudeStep),
                                       'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # Add valid time
        startYear = int(iFile[4:8])
        startDay = int(iFile[8:11])
        self.dataset.SetMetadataItem('time_coverage_start',
            (datetime.datetime(startYear, 1, 1) +
                       datetime.timedelta(startDay)).isoformat())
