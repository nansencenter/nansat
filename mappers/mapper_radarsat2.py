#-------------------------------------------------------------------------------
# Name:        nansat_mapper_radarsat2
# Purpose:     Mapping for Radarsat2 data
#
# Author:      knutfd
# Modified:    mortenwh
#
# Created:     29.11.2011
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from numpy import mod

from vrt import *
from domain import Domain
import tarfile

try:
    from osgeo import gdal
except ImportError:
    import gdal

class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Radarsat2 '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create Radarsat2 VRT '''
        fPathName, fExt = os.path.splitext(fileName)
        if fExt == '.ZIP' or fExt == '.zip':
            fPath, fName = os.path.split(fPathName)
            fileName = '/vsizip/%s/%s' % (fileName, fName)
            gdalDataset = gdal.Open(fileName)
            gdalMetadata = gdalDataset.GetMetadata()

        product = gdalMetadata.get("SATELLITE_IDENTIFIER", "Not_RADARSAT-2")

        #if it is not RADARSAT-2, return
        if product != 'RADARSAT-2':
            raise AttributeError("RADARSAT-2 BAD MAPPER");

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset);

        #define dictionary of metadata and band specific parameters
        pol = []
        metaDict = []

        for i in range(1, gdalDataset.RasterCount+1):
            polString = gdalDataset.GetRasterBand(i).GetMetadata()['POLARIMETRIC_INTERP']
            pol.append(polString)
            metaDict.append(
                {'src': {'SourceFilename':
                         'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml',
                         'SourceBand': i,
                         'DataType': 6},
                 'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', 
                         'BandName': 'sigma0_' + polString,
                         'polarization': polString}})

        ##############################################################
        # Adding derived band (incidence angle) calculated
        # using pixel function "BetaSigmaToIncidence":
        #      incidence = arcsin(sigma0/beta0)*180/pi
        ##############################################################
        metaDict.append(
            {'src': [
                {'SourceFilename': 'RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml',
                 'SourceBand':  1,
                 'DataType': 6},
                {'SourceFilename': 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml',
                 'SourceBand':  1,
                 'DataType': 6}],
            'dst': {'wkv': 'sensor_zenith_angle',
                    'PixelFunctionType': 'BetaSigmaToIncidence',
                    'BandName': 'incidence_angle'}})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

                    
        ###################################################################
        # Add sigma0_VV - pixel function of incidence angle and sigma0_HH
        ###################################################################
        if 'VV' not in pol and 'HH' in pol:        
            sourceBandHH = pol.index('HH')+1
            sourceBandInci = len(metaDict)
            src = [{'SourceFilename': 'RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml',
                    'SourceBand':  sourceBandHH,
                    'DataType': 6},
                   {'SourceFilename': self.fileName,
                    'SourceBand':  sourceBandInci,
                    'DataType': 6}]
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sigma0HHIncidenceToSigma0VV',
                   'polarisation': 'VV',
                   'BandName': 'sigma0_VV'}
            self._create_band(src, dst)
            self.dataset.FlushCache()

        ############################################
        # Add SAR look direction to metadata domain
        ############################################
        self.dataset.SetMetadataItem('SAR_look_direction', str(mod(
            Domain(ds=gdalDataset).upwards_azimuth_direction()
            + 90, 360)))


        # Set time
        validTime = gdalDataset.GetMetadata()['ACQUISITION_START_TIME']
        self.logger.info('Valid time: %s', str(validTime))
        self._set_time(parse(validTime))
