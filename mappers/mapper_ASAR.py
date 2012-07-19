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
from vrt import VRT
from envisat import Envisat

class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1 '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):

        product = gdalMetadata.get("MPH_PRODUCT", "Not_ASAR")

        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        # Get offset for incidence angle, longitude and latitude
        offsetDict = self.read_GeolocationGrid_Offset(fileName)
        # Create a small empty VRT for incident angle
        geoDataset = VRT(srcRasterXSize=11, srcRasterYSize=offsetDict["numOfDSR"])
        # Create a dictionary for VRT RawRasterband
        # dataType "GDT_Float32=6","GDT_Int16=3","GDT_Int32=5"
        parameters = {"ImageOffset" : offsetDict["incidentAngleOffset"],
                       "PixelOffset" : 4,
                       "LineOffset" : offsetDict["DSRsize"],
                       "ByteOrder" : "MSB", "dataType": 6}
        # Add VRT RawRasterband
        geoDataset._create_band(fileName, 0, "", parameters)
        # Enlarge the small band to the size of the underlying data
        self.geoDataset = geoDataset.resized(gdalDataset.RasterXSize, gdalDataset.RasterYSize)

        # Create a dictionary for ASAR band and incident angle
        metaDict = [{'source': fileName, 'sourceBand': 1, 'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', 'parameters':{'band_name': 'sigma0', 'polarisation': gdalMetadata['SPH_MDS1_TX_RX_POLAR'], 'pass': gdalMetadata['SPH_PASS'], 'warning': 'fake sigma0, not yet calibrated'}},
                    {'source': self.geoDataset.fileName, 'sourceBand': 1, 'wkv': '', 'parameters':{'band_name': 'Incidence Angle'}}]

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)
