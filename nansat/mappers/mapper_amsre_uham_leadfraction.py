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

from nansat.exceptions import WrongMapperError

class Mapper(VRT):

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):

        title_correct = False
        if not gdalMetadata:
            raise WrongMapperError
        for key, val in list(gdalMetadata.items()):
            if 'title' in key:
                if not val == 'Daily AMSR-E Arctic lead area fraction [in percent]':
                    raise WrongMapperError
                else:
                    title_correct = True

        if not title_correct:
            raise WrongMapperError

        # initiate VRT for the NSIDC 10 km grid
        self._init_from_dataset_params(1216, 1792, (-3850000, 6250, 0.0, 5850000, 0.0, -6250),
                                       NSR(3411).wkt)

        src = {
            'SourceFilename': 'NETCDF:"%s":lf' % filename,
            'SourceBand': 1,
        }
        dst = {
            'name': 'leadFraction',
            'long_name': 'AMSRE sea ice lead fraction',
        }

        self.create_band(src, dst)
        self.dataset.FlushCache()


