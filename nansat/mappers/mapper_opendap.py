# Name:         mapper_opendap.py
# Purpose:      Mapping for generic OpenDAP online data
# Author:       Anton Korosov
# Created:      01.01.2014
# Copyright:    (c) NERSC
# License:      This file is part of NANSAT. NANSAT is free software: you can
#               redistribute it and/or modify it under the terms of the GNU
#               General Public License as published by the Free Software
#               Foundation, version 3 of the License.
#               http://www.gnu.org/licenses/gpl-3.0.html This program is
#               distributed in the hope that it will be useful, but WITHOUT ANY
#               WARRANTY without even the implied warranty of MERCHANTABILITY
#               or FITNESS FOR A PARTICULAR PURPOSE.
import warnings

import numpy as np

from nansat.tools import gdal, ogr, WrongMapperError
from nansat.vrt import VRT
from nansat.nsr import NSR

try:
    from netCDF4 import Dataset
except ImportError:
    raise ImportError('''
         Cannot import Dataset from netCDF4.
         You cannot access data thought opendap but
         Nansat will work.
         ''')

#fileName = 'http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh/osisaf-nh_aggregated_ice_concentration_nh_polstere-100_197810010000.nc'
#fileName = 'http://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-nat-myoceanv2-20111026'
#fileName = 'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/siw-metno-glo-osisaf/conc/2014/06/ice_conc_sh_polstere-100_multi_201406051200.nc'
#fileName = 'http://thredds.nersc.no/thredds/dodsC/normap/arctic12km_seaice/arctic12km_seaice_20110301_20110331.nc'


class Mapper(VRT):
    def get_proj4_from_ncvar(self, var):
        projDict = {
            'albers_conical_equal_area': {
                0: '+proj=aea',
                'standard_parallel': '+lat_1',
                'longitude_of_central_meridian': '+lon_0',
                'latitude_of_projection_origin': '+lat_0',
                'false_easting': '+x_0',
                'false_northing': '+y_0',
            },
            'polar_stereographic': {
                0: '+proj=stere ',
                'straight_vertical_longitude_from_pole': '+lon_0',
                'standard_parallel': '+lat_1',
                'latitude_of_projection_origin': '+lat_0',
                'scale_factor_at_projection_origin': '+k_0',
                'false_easting': '+x_0',
                'false_northing': '+y_0',
            },
            'stereographic': {
                0: '+proj=stere ',
                'standard_parallel': '+lat_1',
                'longitude_of_projection_origin': '+lon_0',
                'latitude_of_projection_origin': '+lat_0',
                'scale_factor_at_projection_origin': '+k_0',
                'false_easting': '+x_0',
                'false_northing': '+y_0',
            },
            'latitude_longitude': {
                0: '+proj=longlat ',
            }
        }

        attrs = var.ncattrs()
        gmName = str(var.getncattr('grid_mapping_name'))

        if gmName in projDict:
            projSubDict = projDict[gmName]
            proj4 = projSubDict[0]
            for projKey in projSubDict:
                if projKey in attrs:
                    proj4 += (projSubDict[projKey] + '='
                              + str(var.getncattr(projKey)) + ' ')
        return proj4

    def __init__(self, fileName, gdalDataset, gdalMetadata, varName=None, **kwargs):
        ''' Create VRT from OpenDAP dataset'''
        raise WrongMapperError
        # quit if file is not online
        if fileName[:7] not in  ['http://', 'https:/']:
            raise WrongMapperError

        if bandName is None:
            WrongMapperError('Please specify band name')

        # open file through OpenDAP using netCDF4 library
        f = Dataset(fileName)

        # assume CF-compatibility:
        # compulsory grid_mapping_name

        # find grid_mapping_name
        # and get all parameters
        # generate proj4 and WKT strings
        srcProjection = ''
        for varName in f.variables:
            var = f.variables[varName]
            attrs = var.ncattrs()
            if 'grid_mapping_name' in attrs:
                proj4str = self.get_proj4_from_ncvar(var)
                srcProjection = NSR(proj4str).wkt
                break

        print 'WKT:', srcProjection

        xDim = None
        yDim = None
        validVars = []
        validDims = []

        # if grid_mapping_name is found: find all bands that have grid_mapping
        # else: find lon/lat variables that have 1D
        import ipdb; ipdb.set_trace()

        if srcProjection != '':
            # find bands with grid_mapping
            var = f.variables[varName]
            attrs = var.ncattrs()
            if 'grid_mapping' in attrs:
                validVars.append(str(varName))
                validDims += [str(dim) for dim in var.dimensions]

            # assume NORMAP compatibility:
            # dimensions should be
            # x, or xc, or lon, or longitude
            # and
            # y, or yc, or lat, or latitude
            # and/or
            # time
            # and/or
            # depth
            # and/or
            # etc

            # assign x, y dimension names
            for dim in validDims:
                if 'x' in dim or 'lon' in dim or 'east' in dim:
                    xDim = dim
                if 'y' in dim or 'lat' in dim or 'north' in dim:
                    yDim = dim

        else:
            # if input file is not CF-compliant but has 1D lon/lat dimensions
            # find 1D lon, lat (longitude, latitude) dims
            var = f.variables[varName]
            if 'lon' in varName and var.ndim == 1:
                xDim = varName
            if 'lat' in varName and var.ndim == 1:
                yDim = varName

            # find all datasets with lon/lat dimensions
            for varName in f.variables:
                var = f.variables[varName]
                if xDim in var.dimensions and yDim in var.dimensions:
                    validVars.append(str(varName))
                    validDims += [str(dim) for dim in var.dimensions]

            if len(validVars) > 1:
                srcProjection = NSR().wkt
            else:
                # cancel usage of this mapper a input file seem unappropriate
                raise


        # get X/Y size
        var0 = f.variables[validVars[0]]
        xDimI = var0.dimensions.index(xDim)
        srcRasterXSize = var0.shape[xDimI]
        yDimI = var0.dimensions.index(yDim)
        srcRasterYSize = var0.shape[yDimI]

        print 'x/ySize', srcRasterXSize, srcRasterYSize

        # get GDAL GeoTransform
        xdata = f.variables[xDim][:]
        ydata = f.variables[yDim][:]
        x0 = xdata[0]
        dx = xdata[1] - xdata[0]
        y0 = ydata[0]
        dy = ydata[1] - ydata[0]
        srcGeoTransform = (x0, dx, 0, y0, 0, dy)
        print srcGeoTransform

        # make list of metadata dictionary entries
        # for each var make [k][i][j][x][y] depending on order of dims
        metaDict = []
        for varName in validVars:
            var = f.variables[varName]
            dims = var.dimensions

            # find nor X neither Y dimensions (e.g. time, depth, etc)
            nonxyDims = []
            nonxyDimValues = []
            for dim in dims:
                if dim != xDim and dim != yDim:
                    nonxyDims.append(dim)
                    dimValues = f.variables[dim][:]
                    dimValuesStr = [str(val) for val in dimValues]
                    nonxyDimValues.append(dimValues)

            # calculate total dimensionality in addition to X/Y
            # make vector of nonXY dimensions
            nonXYShape = []
            for nonxyDim in nonxyDims:
                dimVar = f.variables[nonxyDim]
                nonXYShape.append(dimVar.shape[0])

            if len(nonXYShape) == 0:
                nonXYShape = [1]
            totNonXYDims = np.cumprod(nonXYShape)[-1]

            urls = []
            # generate bands for each additional dimension
            for nonxyi in range(totNonXYDims):
                dstVarName = str(varName)
                url = fileName + '?%s.%s' % (varName, varName)
                #assert varName != 'v'
                # vector of nonX/Y indeces
                iVec = np.unravel_index(nonxyi, nonXYShape)
                # dim metadata keeps name of dimension and index
                # in this dimension (value)
                dimMetadata = {}
                # add either [x], or [y], or respective index in each dimension
                for dim in dims:
                    if dim == xDim:
                        url = url + '[x]'
                    elif dim == yDim:
                        url += '[y]'
                    else:
                        # get index of that dimension
                        dimN = nonxyDims.index(dim)
                        dimV = iVec[dimN]
                        url += '[%d]' % dimV
                        dimMetadata[dim] = nonxyDimValues[dimN][dimV]
                        dstVarName += '%03d' % dimV

                # get band metadata
                attrs = var.ncattrs()

                metaEntry = {'src': {'SourceFilename': url, 'sourceBand':  1},
                             'dst': {'name': dstVarName}
                             }

                # put band metadata
                for attr in attrs:
                    attrKey = attr.encode('ascii', 'ignore')
                    attrVal = var.getncattr(attr)
                    if type(attrVal) in [str, unicode]:
                        attrVal = attrVal.encode('ascii', 'ignore')
                    else:
                        attrVal = str(attrVal)
                    metaEntry['dst'][attrKey] = attrVal

                # add wkv
                if 'standard_name' in attrs:
                    metaEntry['dst']['wkv'] = metaEntry['dst']['standard_name']

                # add dim metadata (location within each dimension)
                for dimKey in dimMetadata:
                    metaEntry['dst'][str(dimKey)] = dimMetadata[dimKey]

                metaDict.append(metaEntry)

        # read global metadata
        srcMetadata = {}
        for attr in f.ncattrs():
            srcMetadata[str(attr)] = str(attr)

        # create VRT with bands
        VRT.__init__(self, srcGeoTransform=srcGeoTransform,
                     srcProjection=srcProjection,
                     srcRasterXSize=srcRasterXSize,
                     srcRasterYSize=srcRasterYSize,
                     srcMetadata=srcMetadata)
        self._create_bands(metaDict)
