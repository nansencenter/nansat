from vrt import VRT, gdal, osr, latlongSRS, np
from netCDF4 import Dataset

#f = Dataset()
#furl = 'http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh/osisaf-nh_aggregated_ice_concentration_nh_polstere-100_197810010000.nc'
#furl = 'http://thredds.nersc.no/thredds/dodsC/normap/arctic12km_seaice/arctic12km_seaice_19870801_19870831.nc'

class Mapper(VRT):

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 GCP_COUNT=10, **kwargs):

        assert fileName[:7] == 'http://'

        f = Dataset(fileName)

        # assume CF-compatibility:
        # compulsory grid_mapping_name
        
        # find grid_mapping_name
        # and get all parameters
        for varName in f.variables:
            var = f.variables[varName]
            attrs = var.ncattrs()
            if 'grid_mapping_name' in attrs:
                gmName = str(varName)
                gmVal = str(var.getncattr('grid_mapping_name'))
                break

        # find bands with grid_mapping
        validVars = []
        validDims = []
        for varName in f.variables:
            var = f.variables[varName]
            attrs = var.ncattrs()
            if 'grid_mapping' in attrs:
                validVars.append(str(varName))
                for dim in var.dimensions:
                    validDims.append(str(dim))
        
        validDims = list(set(validDims))
        print validVars, validDims
        
        # assume NORMAP compatibility:
        # diemnsions should be
        # x, or xc, or lon, or longitude
        # and
        # y, or yc, or lat, or latitude
        # and/or
        # time
        # and/or
        # depth
        # and/or
        # etc

        # assign time, x, y dimension names
        xDim = 'None'
        yDim = 'None'
        for dim in validDims:
            if 'x' in dim or 'lon' in dim or 'east' in dim:
                xDim = dim
            if 'y' in dim or 'lat' in dim or 'north' in dim:
                yDim = dim

        # make list of metadata dictionary entries
        # for each var make [k][i][j][x][y] depending on order of dims
        metaDict = []
        for varName in validVars:
            var = f.variables[varName]
            dims = var.dimensions
            
            # find nor X neither Y dimensions (e.g. time, depth, etc)
            nonxyDims = []
            for dim in dims:
                if dim != xDim and dim != yDim:
                    nonxyDims.append(dim)
            
            # calculate total dimensionality in addition to X/Y
            # make vector of nonXY dimensions
            nonXYShape = []
            for nonxyDim in nonxyDims:
                dimVar = f.variables[nonxyDim]
                nonXYShape.append(dimVar.shape[0])
        
            if len(nonXYShape) == 0:
                nonXYShape = [1]
            #nonXYShape = [3, 5]
            totNonXYDims = np.cumprod(nonXYShape)[-1]
            
            urls = []
            # generate bands for each additional dimension
            for nonxyi in range(totNonXYDims):
                url = fileName + '?%s.%s' % (varName, varName)
                # vector of nonX/Y indeces
                iVec = np.unravel_index(nonxyi, nonXYShape)
                # add either [x], or [y], or respective index in each dimension
                for dim in dims:
                    if dim == xDim:
                        url = url + '[x]'
                    elif dim == yDim:
                        url += '[y]'
                    else:
                        dimN = nonxyDims.index(dim)
                        dimV = iVec[dimN]
                        url += '[%d]' % dimV

                # get band metadata
                attrs = var.ncattrs()
                
                metaEntry = {'src': {'SourceFilename': url, 'sourceBand':  1},
                             'dst': {'name': varName}
                             }
            
                # put band metadata
                for attr in attrs:
                    print attr
                    metaEntry['dst'][str(attr)] = str(var.getncattr(attr))

                # add wkv
                if 'standard_name' in attrs:
                    metaEntry['dst']['wkv'] = metaEntry['dst']['standard_name']
            
                print metaEntry
                metaDict.append(metaEntry)

        print metaDict[0]['src']['SourceFilename']
        gdalDataset = gdal.Open(metaDict[0]['src']['SourceFilename'])
        VRT.__init__(self, gdalDataset)
        self._create_bands(metaDict)
