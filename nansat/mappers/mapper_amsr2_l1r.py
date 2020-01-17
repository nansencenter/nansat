# Name:        mapper_obpg_l2
# Purpose:     Mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
from datetime import datetime, timedelta
from math import ceil
import json

import pythesint as pti

from nansat.utils import gdal, ogr, parse_time
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT
from nansat.nsr import NSR


class Mapper(VRT):
    ''' Mapper for AMSR-2 L1 data

    '''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 GCP_STEP=20, MAX_LAT=90, MIN_LAT=50, resolution='low',
                 **kwargs):
        ''' Create VRT
        Parameters
        ----------
        GCP_COUNT : int
            number of GCPs along each dimention
        '''
        ifile = os.path.split(filename)[1]
        if not ifile.startswith('GW1AM2_') or not ifile.endswith('.h5'):
            raise WrongMapperError
        try:
            ProductName = gdalMetadata['ProductName']
            PlatformShortName = gdalMetadata['PlatformShortName']
            SensorShortName = gdalMetadata['SensorShortName']
        except:
            raise WrongMapperError

        if (not ProductName == 'AMSR2-L1R' or
            not PlatformShortName == 'GCOM-W1' or
            not SensorShortName == 'AMSR2'):
            raise WrongMapperError

        if resolution == 'low':
            subDatasetWidth = 243
        else:
            subDatasetWidth = 486

        # get GCPs from lon/lat grids
        latGrid = gdal.Open('HDF5:"%s"://Latitude_of_Observation_Point_for_89A' % filename).ReadAsArray()
        lonGrid = gdal.Open('HDF5:"%s"://Longitude_of_Observation_Point_for_89A' % filename).ReadAsArray()
        if subDatasetWidth == 243:
            latGrid = latGrid[:, ::2]
            lonGrid = lonGrid[:, ::2]

        dx = .5
        dy = .5
        gcps = []
        k = 0
        maxY = 0
        minY = latGrid.shape[0]
        for i0 in range(0, latGrid.shape[0], GCP_STEP):
            for i1 in range(0, latGrid.shape[1], GCP_STEP):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(lonGrid[i0, i1])
                lat = float(latGrid[i0, i1])
                if (lon >= -180 and
                    lon <= 180 and
                    lat >= MIN_LAT and
                    lat <= MAX_LAT):
                    gcp = gdal.GCP(lon, lat, 0, i1 + dx, i0 + dy)
                    gcps.append(gcp)
                    k += 1
                    maxY = max(maxY, i0)
                    minY = min(minY, i0)
        yOff = minY
        ySize = maxY - minY

        # remove Y-offset from gcps
        for gcp in gcps:
            gcp.GCPLine -= yOff

        metaDict = []

        subDatasets = gdalDataset.GetSubDatasets()
        metadata = gdalDataset.GetMetadata()
        for subDataset in subDatasets:
            # select subdatasets fro that resolution (width)
            if (subDatasetWidth == int(subDataset[1].split(']')[0].split('x')[-1]) and
                'Latitude' not in subDataset[0] and 'Longitude' not in subDataset[0]):
                name = subDataset[0].split('/')[-1]
                # find scale
                scale = 1
                for meta in metadata:
                    if name + '_SCALE' in meta:
                        scale = float(metadata[meta])
                # create meta entry
                metaEntry = {'src': {'SourceFilename': subDataset[0],
                                     'sourceBand':  1,
                                     'ScaleRatio': scale,
                                     'ScaleOffset': 0,
                                     'yOff': yOff,
                                     'ySize': ySize,},
                             'dst': {'name': name}
                             }
                metaDict.append(metaEntry)

        # create VRT from one of the subdatasets
        gdalSubDataset = gdal.Open(metaEntry['src']['SourceFilename'])
        self._init_from_dataset_params(subDatasetWidth, ySize, (1,0,0,ySize,0,-1), NSR().wkt)
        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        self.dataset.SetMetadataItem('time_coverage_start',
                    parse_time(gdalMetadata['ObservationStartDateTime']).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end',
                    parse_time(gdalMetadata['ObservationEndDateTime']).isoformat())
        # append GCPs and lat/lon projection to the vsiDataset
        self.dataset.SetGCPs(gcps, NSR().wkt)
        self.reproject_gcps('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
        self.tps = True

        mm = pti.get_gcmd_instrument('AMSR2')
        ee = pti.get_gcmd_platform('GCOM-W1')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
