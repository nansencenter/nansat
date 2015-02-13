#-------------------------------------------------------------------------------
# Name:		mapper_oceansat2_l2c.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	13.02.2015
# Last modified:13.02.2015 11:38
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
from osgeo import gdal

from nansat.vrt import VRT

class Mapper(VRT):

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        if 'Title' not in gdalMetadata.keys():
            raise WrongMapperError
        if not 'Oceansat OCM2 Level-2C' in gdalMetadata['Title']:
            raise WrongMapperError

        subDatasets = gdalDataset.GetSubDatasets()
        gds = gdal.Open(subDatasets[0][0])
        meta = gds.GetMetadata()
        data = gds.GetRasterBand(1).ReadAsArray()
        cols = data.shape[1]
        rows = data.shape[0]

        ullon = float(meta['Upper Left Longitude'])
        urlon = float(meta['Upper Right Longitude'])
        pwidth = (urlon - ullon)/cols

        ullat = float(meta['Upper Left Latitude'])
        lllat = float(meta['Lower Left Latitude'])
        pheight = (ullat - lllat)/rows

        srcGeoTransform = ( 
                ullon, pwidth, 0,
                ullat, 0, pheight,
            )
        import ipdb
        ipdb.set_trace()
        pass
