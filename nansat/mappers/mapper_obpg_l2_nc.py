# Name:        mapper_obpg_l2
# Purpose:     Mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    raise ImportError('''
         Cannot import Dataset from netCDF4.
         You cannot access Oceancolor NC4 data but
         Nansat will work.
         ''')

from nansat.tools import gdal, ogr, WrongMapperError
from nansat.vrt import GeolocationArray, VRT
from nansat.nsr import NSR
from nansat.mappers.obpg import OBPGL2BaseClass

class Mapper(OBPGL2BaseClass):
    ''' Mapper for SeaWIFS/MODIS/MERIS/VIIRS L2 data from OBPG in NC4 format
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 GCP_COUNT=10, **kwargs):
        ''' Create VRT
        Parameters
        ----------
        GCP_COUNT : int
            number of GCPs along each dimention
        '''
        # file must be opened by netCDF4 lib
        try:
            ds = Dataset(fileName)
        except RuntimeError:
            raise WrongMapperError

        # 'title' should be in attributes
        try:
            title = ds.title
        except:
            raise WrongMapperError

        # title value should be prelisted
        if title not in self.titles:
            raise WrongMapperError

        # get geophysical data variables
        gdGroup = ds.groups['geophysical_data']
        metaDict = []
        for varName in gdGroup.variables:
            var = gdGroup.variables[varName]
            dataType = {'uint8': gdal.GDT_Byte,
                        'int8':  gdal.GDT_Byte,
                        'uint16': gdal.GDT_UInt16,
                        'int16':  gdal.GDT_Int16,
                        'uint32': gdal.GDT_UInt32,
                        'int32':  gdal.GDT_Int32,
                        'float32': gdal.GDT_Float32,
                        'float64': gdal.GDT_Float64,
                        'complex64': gdal.GDT_CFloat32,
                        'complex128': gdal.GDT_CFloat64}.get(str(var.dtype))

            # set SRC/DST parameters
            metaEntry = {'src': {'SourceFilename': 'HDF5:"%s"://geophysical_data/%s' % (fileName, varName),
                                 'sourceBand': 1,
                                 'DataType': dataType},
                         'dst': {'name': varName}}
            # set standard_name if exists
            if 'standard_name' in var.ncattrs():
                metaEntry['dst']['wkv'] = var.standard_name
            # set scale if exist
            if 'scale_factor' in  var.ncattrs():
                metaEntry['src']['ScaleRatio'] = var.scale_factor
            # set offset if exist
            if 'add_offset' in  var.ncattrs():
                metaEntry['src']['ScaleOffset'] = var.add_offset
            metaDict.append(metaEntry)

        # make GCPs
        # get control points columns
        cntl_pt_cols = ds.groups['navigation_data'].variables['cntl_pt_cols'][:]
        rasterYSize, rasterXSize = var.shape
        gcpCols, gcpRows = np.meshgrid(cntl_pt_cols, range(rasterYSize))
        longitude = ds.groups['navigation_data'].variables['longitude'][:]
        latitude = ds.groups['navigation_data'].variables['latitude'][:]

        step0 = max(1, int(float(latitude.shape[0]) / GCP_COUNT))
        step1 = max(1, int(float(latitude.shape[1]) / GCP_COUNT))

        gcps = []
        k = 0
        for i0 in range(0, longitude.shape[0], step0):
            for i1 in range(0, latitude.shape[1], step1):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(longitude[i0, i1])
                lat = float(latitude[i0, i1])
                pixel = float(gcpCols[i0, i1]) + .5
                line = float(gcpRows[i0, i1]) + .5

                if (lon >= -180 and lon <= 180 and lat >= -90 and lat <= 90):
                    gcp = gdal.GCP(lon, lat, 0, pixel, line)
                    gcps.append(gcp)
                    k += 1

        # clean up netCDF4.Dataset
        gdGroup = None
        ds = None
        var = None

        VRT.__init__(self, srcProjection=NSR().wkt,
                     srcGCPs=gcps,
                     srcGCPProjection=NSR().wkt,
                     srcRasterXSize=rasterXSize,
                     srcRasterYSize=rasterYSize)

        self._create_bands(metaDict)
