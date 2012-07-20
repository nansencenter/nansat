#-------------------------------------------------------------------------------
# Name:        nansat_mapper_asar
# Purpose:     Mapping for Envisat ASAR-L1 data
#
# Author:      Knut-Frode
#
# Created:     13.12.2011
# Copyright:
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import VRT, Geolocation
from envisat import Envisat

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        product = gdalMetadata.get("MPH_PRODUCT", "Not_ASAR")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")
        # prepare parameters to create a dataset with a small band from geolocation array
        gadsDSName = 'DS_NAME="GEOLOCATION GRID ADS        "\n'
        parameters = [{"DS_OFFSET": 3,
                       "substream" : {"incidentAngle": 25+(4*2)*11}},
                      {"NUM_DSR" : 5},
                      {"DSR_SIZE" : 6}]
        dataType = "float32"

        # create a data set with a small band for incident Angle
        incAngleDataset = self.create_geoDataset(fileName, gadsDSName, parameters, dataType)

        geoDataset.dataset.FlushCache()
        # Enlarge the small band to the size of the underlying data
        self.incAngleDataset = incAngleDataset.resized(gdalDataset.RasterXSize, gdalDataset.RasterYSize)

        # Create a dictionary for ASAR band and incident angle
        metaDict = [{'source': fileName, 'sourceBand': 1, 'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', 'parameters':{'band_name': 'sigma0', 'polarisation': gdalMetadata['SPH_MDS1_TX_RX_POLAR'], 'pass': gdalMetadata['SPH_PASS'], 'warning': 'fake sigma0, not yet calibrated'}},
                    {'source': self.incAngleDataset.fileName, 'sourceBand': 1, 'wkv': '', 'parameters':{'band_name': 'Incidence Angle'}}]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)
        '''
        # prepare parameters to create a dataset from longitude and latitude
        gadsDSName = 'DS_NAME="GEOLOCATION GRID ADS        "\n'
        parameters = [{"DS_OFFSET": 3,
                       "substream" : {"longitude": 25+(4*4)*11, "latitude": 25+(4*3)*11}},
                      {"NUM_DSR" : 5},
                      {"DSR_SIZE" : 6}]
        dataType = "int32"

        # create dataset for longitude and
        self.geoDataset = self.create_geoDataset(fileName, gadsDSName, parameters, dataType)
        band = self.geoDataset.dataset.GetRasterBand(1)
        if band.GetMetadata()["band_name"] == "longitude":
            xBand = 1
            yBand = 2
        else:
            xBand = 2
            yBand = 1

        self.add_geolocation(Geolocation(xVRT=self.geoDataset.fileName,
                                  yVRT=self.geoDataset.fileName,
                                  xBand=xBand, yBand=yBand,
                                  srs = gdalDataset.GetGCPProjection()))
        '''
