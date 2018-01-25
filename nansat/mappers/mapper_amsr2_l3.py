# Name:        mapper_masr2_l3
# Purpose:     Mapping for L3 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse
import datetime
import os.path
import glob
import json

import pythesint as pti

import numpy as np

from nansat.vrt import VRT
from nansat.nsr import NSR
from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError


class Mapper(VRT):
    ''' Mapper for Level-3 AMSR2 data from https://gcom-w1.jaxa.jp'''

    freqs = [6, 7, 10, 18, 23, 36, 89]

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' OBPG L3 VRT '''

        # test the product
        try:
            assert gdalMetadata['PlatformShortName'] == 'GCOM-W1'
            assert gdalMetadata['SensorShortName'] == 'AMSR2'
            assert gdalMetadata['ProductName'] == 'AMSR2-L3'
        except:
            raise WrongMapperError

        # get list of similar (same date, A/D orbit) files in the directory
        iDir, iFile = os.path.split(filename)
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
                    metaEntry = {
                        'src': {'SourceFilename': simSubDataset[0],
                                'SourceBand': 1,
                                'ScaleRatio': 0.0099999998,
                                'ScaleOffset': 0},
                        'dst': {'wkv': 'brightness_temperature',
                                'frequency': '%02d' % freq,
                                'polarisation': simSubDataset[0][-2:-1],
                                'suffix': ('%02d%s' %
                                           (freq, simSubDataset[0][-2:-1]))}}
                    metaDict.append(metaEntry)

        # initiate VRT for the NSIDC 10 km grid
        self._init_from_dataset_params(760, 1120, (-3850000, 10000, 0.0, 5850000, 0.0, -10000), NSR(3411).wkt)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start', parse(gdalMetadata['ObservationStartDateTime']).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', parse(gdalMetadata['ObservationStartDateTime']).isoformat())

        mm = pti.get_gcmd_instrument('AMSR2')
        ee = pti.get_gcmd_platform('GCOM-W1')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

