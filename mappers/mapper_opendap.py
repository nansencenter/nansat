from pydap.client import open_url
from vrt import VRT, Geolocation
import matplotlib.pyplot as plt

import osr

class Mapper(VRT):
    ''' VRT with mapping for OpenDAP '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        '''Create VRT from OpenDAP url
        '''
        #fileName = 'http://thredds.nersc.no/thredds/dodsC/normap/test1/westernnorway300m_oceancolor_20111001-20111031.nc'
        #fileName = 'http://thredds.met.no/thredds/dodsC/cryoclim/met.no/osisaf-nh/osisaf-nh_aggregated_ice_concentration_nh_polstere-100_200712010000.nc'
        # OpenDAP object
        od = open_url(fileName)
        
        self.odVRTDict = {}
        metaDict = []
        gridMappingName = ''
        dim0 = 0
        dim1 = 0
        # choose 2D (3D) arrays
        for bandName in od.keys():
            odBand = od[bandName]
            bandShape = odBand.shape
            bandMetadata = {}
            # convert metadata to str
            for metaName in odBand.attributes:
                bandMetadata[metaName] = str(odBand.attributes[metaName])
                if metaName == 'grid_mapping':
                    gridMappingName = odBand.attributes[metaName]
            # set the dims of the dataset to the dims of the first 2D band
            if len(bandShape) == 2 and bandShape[0] > 1 and bandShape[1] > 1 and dim0 == 0 and dim1 == 0:
                    dim0 = bandShape[0]
                    dim1 = bandShape[1]
                    dim0Name = odBand.dimensions[0]
                    dim1Name = odBand.dimensions[1]
            if len(bandShape) == 3 and bandShape[1] > 1 and bandShape[2] > 1 and dim0 == 0 and dim1 == 0:
                    dim0 = bandShape[1]
                    dim1 = bandShape[2]
                    dim0Name = odBand.dimensions[1]
                    dim1Name = odBand.dimensions[2]

            # get the band if dims are equal to the dims of the first band
            if len(bandShape) == 2 and bandShape[0] == dim0 and bandShape[1] == dim1:
                subBandName = bandName
                bandMetadata['band_name'] = subBandName
                # read array from thredds and create VRT
                self.odVRTDict[subBandName] = VRT(array = odBand.array[:].astype('float32'))
                metaDict.append(
                    {'source': self.odVRTDict[subBandName].fileName,
                     'sourceBand': 1,
                     'wkv': '',
                     'parameters': bandMetadata})
            if len(bandShape) == 3 and bandShape[1] == dim0 and bandShape[2] == dim1:
                for i in range(0, bandShape[0]):
                    subBandName = bandName + '_%03d' % i
                    bandMetadata['band_name'] = subBandName
                    # read array from thredds and create VRT                    
                    self.odVRTDict[subBandName] = VRT(array = odBand.array[i, :, :].reshape(bandShape[1], bandShape[2]).astype('float32'))
                    metaDict.append(
                        {'source': self.odVRTDict[subBandName].fileName,
                         'sourceBand': 1,
                         'wkv': '',
                         'parameters': bandMetadata})
                    
        # get variable with grid_mapping 
        try:
            gridMappingVar = od[gridMappingName]
        except:
            self.logger.warning('No grid_mapping!')
        else:
            gmMetadata = gridMappingVar.attributes
        
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, self.odVRTDict[subBandName].dataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set projection
        sr = osr.SpatialReference()
        if 'proj4_string' in gmMetadata:
            sr.ImportFromProj4(gmMetadata['proj4_string'])
        elif 'spatial_ref' in gmMetadata:
            sr.ImportFromWkt(gmMetadata['spatial_ref'])
        else:
            sr.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
        srs = sr.ExportToWkt()
        self.dataset.SetProjection(srs)

        # set geoTransform (from grid_mapping metadata or from variable vectors)
        if 'GeoTransform' in gmMetadata:
            geoTransform = tuple(map(int, gmMetadata['spatial_ref'].split(' ')))
        else:
            d0Vec = od[dim0Name][:].astype('float32') # Y-axis
            d1Vec = od[dim1Name][:].astype('float32') # X-axis
            # convert from km to m
            if od[dim1Name].attributes.get('units', '') == 'km':
                d0Vec *= 1000.
                d1Vec *= 1000.
            geoTransform = (d1Vec[0], d1Vec[1] - d1Vec[0], 0, d0Vec[0], 0, d0Vec[1]-d0Vec[0])
        self.dataset.SetGeoTransform(geoTransform)
        
