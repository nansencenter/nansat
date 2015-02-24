#-------------------------------------------------------------------------------
# Name:		mapper_amsre_UHAM_lead_fraction.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	18.02.2015
# Last modified:24.02.2015 09:26
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import datetime
from osgeo import gdal, osr
from nansat.nsr import NSR
from nansat.vrt import VRT

from nansat.tools import WrongMapperError

class Mapper(VRT):

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        
        title_correct = False
        if not gdalMetadata:
            raise WrongMapperError
        for key, val in gdalMetadata.iteritems():
            if 'title' in key:
                if not val == \
                        'Daily AMSR-E Arctic lead area fraction [in percent]':
                    raise WrongMapperError
                else:
                    title_correct = True

        if not title_correct:
            raise WrongMapperError

        # initiate VRT for the NSIDC 10 km grid
        VRT.__init__(self,
                     srcGeoTransform=(-3850000, 6250, 0.0,
                                      5850000, 0.0, -6250),
                     srcProjection=NSR(3411).wkt,
                     srcRasterXSize=1216,
                     srcRasterYSize=1792)

        src = {
            'SourceFilename': 'NETCDF:"%s":lf'%fileName,
            'SourceBand': 1,
        }
        dst = {
            'name': 'leadFraction',
            'long_name': 'AMSRE sea ice lead fraction',
        }

        self._create_band(src, dst)
        self.dataset.FlushCache()


