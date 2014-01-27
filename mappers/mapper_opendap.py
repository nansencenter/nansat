from vrt import GeolocationArray, VRT, gdal, osr, latlongSRS
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from osgeo import gdal

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
        
        print gmName, gmVal
        

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
        # for each var make [0][x][y] depending on order of dims
        metaDict = []
        for varName in validVars:
            var = f.variables[varName]
            dims = var.dimensions
            url = fileName + '?%s.%s' % (varName, varName)
            
            for dim in dims:
                if xDim in dim:
                    url += '[x]'
                elif yDim in dim:
                    url += '[y]'
                else:
                    url += '[0]'
            
            attrs = var.ncattrs()
            
            metaEntry = {'src': {'SourceFilename': url, 'sourceBand':  1},
                         'dst': {'name': varName}
                         }
        
            for attr in attrs:
                print attr
                metaEntry['dst'][str(attr)] = str(var.getncattr(attr))
            if 'standard_name' in attrs:
                metaEntry['dst']['wkv'] = metaEntry['dst']['standard_name']
        
            print metaEntry
            metaDict.append(metaEntry)

        print metaDict[0]['src']['SourceFilename']
        gdalDataset = gdal.Open(metaDict[0]['src']['SourceFilename'])
        VRT.__init__(self, gdalDataset)
        self._create_bands(metaDict)
