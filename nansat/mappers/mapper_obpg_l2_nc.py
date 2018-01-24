# Name:        mapper_obpg_l2
# Purpose:     Mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
from dateutil.parser import parse
import json

import numpy as np

import pythesint as pti

from nansat.tools import gdal, ogr
from nansat.vrt import VRT
from nansat.nsr import NSR
from nansat.mappers.obpg import OBPGL2BaseClass

from nansat.exceptions import WrongMapperError

class Mapper(OBPGL2BaseClass):
    ''' Mapper for SeaWIFS/MODIS/MERIS/VIIRS L2 data from OBPG in NC4 format
    '''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 GCP_COUNT=10, **kwargs):
        ''' Create VRT
        Parameters
        ----------
        GCP_COUNT : int
            number of GCPs along each dimention
        '''

        # extension must be .nc
        if os.path.splitext(filename)[1] != '.nc':
            raise WrongMapperError

        # file must contain navigation_data/longitude
        try:
            ds = gdal.Open('HDF5:"%s"://navigation_data/longitude' % filename)
        except RuntimeError:
            raise WrongMapperError
        else:
            dsMetadata = ds.GetMetadata()

        # title value must be known
        if dsMetadata.get('title', '') not in self.titles:
            raise WrongMapperError

        # get geophysical data variables
        subDatasets = gdal.Open(filename).GetSubDatasets()
        metaDict = []
        for subDataset in subDatasets:
            groupName = subDataset[0].split('/')[-2]
            if groupName not in ['geophysical_data', 'navigation_data']:
                continue
            varName = subDataset[0].split('/')[-1]
            subds = gdal.Open(subDataset[0])
            b = subds.GetRasterBand(1)
            bMetadata = b.GetMetadata()

            # set SRC/DST parameters
            metaEntry = {'src': {'SourceFilename': subDataset[0],
                                 'sourceBand': 1,
                                 'DataType': b.DataType},
                         'dst': {'name': varName}}

            # replace datatype for l2_flags
            if varName == 'l2_flags':
                metaEntry['src']['DataType'] = 4
                metaEntry['src']['SourceType'] = 'SimpleSource'

            # set scale if exist
            metaKey = '%s_%s_scale_factor' % (groupName, varName)
            if metaKey in bMetadata:
                metaEntry['src']['ScaleRatio'] = bMetadata[metaKey]

            # set offset if exist
            metaKey = '%s_%s_add_offset' % (groupName, varName)
            if metaKey in bMetadata:
                metaEntry['src']['ScaleOffset'] = bMetadata[metaKey]

            # set standard_name if exists
            metaKey = '%s_%s_standard_name' % (groupName, varName)
            if metaKey in bMetadata:
                metaEntry['dst']['wkv'] = bMetadata[metaKey]

            # set other metadata
            for metaKey in bMetadata:
                newMetaKey = metaKey.replace('%s_%s_' %  (groupName, varName), '')
                if newMetaKey not in ['scale_factor', 'add_offset', 'DIMENSION_LIST', '_FillValue']:
                    metaEntry['dst'][newMetaKey] = bMetadata[metaKey]
            metaDict.append(metaEntry)

        # make GCPs
        # get lat/lon grids
        longitude = gdal.Open('HDF5:"%s"://navigation_data/longitude' % filename).ReadAsArray()
        latitude = gdal.Open('HDF5:"%s"://navigation_data/latitude' % filename).ReadAsArray()
        rasterYSize, rasterXSize = longitude.shape

        step0 = max(1, int(float(latitude.shape[0]) / GCP_COUNT))
        step1 = max(1, int(float(latitude.shape[1]) / GCP_COUNT))

        gcps = []
        k = 0
        center_lon = 0
        center_lat = 0
        for i0 in range(0, latitude.shape[0], step0):
            for i1 in range(0, latitude.shape[1], step1):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(longitude[i0, i1])
                lat = float(latitude[i0, i1])

                if (lon >= -180 and lon <= 180 and lat >= -90 and lat <= 90):
                    gcp = gdal.GCP(lon, lat, 0, i1+0.5, i0+0.5)
                    gcps.append(gcp)
                    center_lon += lon
                    center_lat += lat
                    k += 1

        time_coverage_start = dsMetadata['time_coverage_start']
        time_coverage_end = dsMetadata['time_coverage_end']

        # create VRT
        # x_size, y_size, geo_transform, projection, gcps=None, gcp_projection='', **kwargs
        self._init_from_dataset_params(rasterXSize, rasterYSize, (0, 1, 0, rasterYSize, 0, -1),
                                        NSR().wkt, gcps, NSR().wkt)
        # add bands
        self.create_bands(metaDict)

        # reproject GCPs
        center_lon /= k
        center_lat /= k
        srs = '+proj=stere +datum=WGS84 +ellps=WGS84 +lon_0=%f +lat_0=%f +no_defs' % (center_lon, center_lat)
        self.reproject_GCPs(srs)

        ### BAD, BAd, bad ...
        self.dataset.SetProjection(self.dataset.GetGCPProjection())

        # use TPS for reprojection
        self.tps = True

        # add NansenCloud metadata
        self.dataset.SetMetadataItem('time_coverage_start',
                                     str(time_coverage_start))
        self.dataset.SetMetadataItem('time_coverage_end',
                                     str(time_coverage_end))
        self.dataset.SetMetadataItem('source_type', 'Satellite')
        self.dataset.SetMetadataItem('mapper', 'obpg_l2_nc')


        platform  = {'Orbview-2': 'SEASTAR'}.get(dsMetadata.get('platform'),
                                                 dsMetadata.get('platform'))


        mm = pti.get_gcmd_instrument(dsMetadata.get('instrument'))
        ee = pti.get_gcmd_platform(platform)
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
