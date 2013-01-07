#-------------------------------------------------------------------------------
# Name:        nansat_mapper_radarsat2
# Purpose:     Mapping for Radarsat2 data
#
# Author:      knutfd
# Modified:	Morten Wergeland Hansen
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

import pdb

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

        # Get the subdataset with calibrated sigma0 only
        for dataset in gdalDataset.GetSubDatasets():
            if dataset[1]=='Sigma Nought calibrated':
                s0dataset = gdal.Open(dataset[0])
                s0datasetName = dataset[0][:]
                s0datasetPol = s0dataset.GetRasterBand(1).GetMetadata()['POLARIMETRIC_INTERP']
                print 's0datasetName', s0datasetName
                print 's0datasetPol', s0datasetPol
                for i in range(1, s0dataset.RasterCount+1):
                    polString = s0dataset.GetRasterBand(i).GetMetadata()['POLARIMETRIC_INTERP']
                    ''' 
                    Make complex bands if the SAR data is of type 10, and 
                    adjust methods to also work for complex numbers... 
                    '''
                    #if s0dataset.GetRasterBand(i).DataType==10:
                    #    suffix = polString+'_complex'
                    #else:
                    suffix = polString
                    pol.append(polString)
                    metaDict.append(
                        {'src': {'SourceFilename':
                            'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml',
                            'SourceBand': i,
                            'DataType': s0dataset.GetRasterBand(i).DataType},
                        'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', 
                            'suffix': suffix,
                            'polarization': polString}})
            if dataset[1]=='Beta Nought calibrated':
                b0dataset = gdal.Open(dataset[0])
                b0datasetName = dataset[0][:]
                for j in range(1, b0dataset.RasterCount+1):
                    polString = b0dataset.GetRasterBand(j).GetMetadata()['POLARIMETRIC_INTERP']
                    #pdb.set_trace()
                    if polString==s0datasetPol:
                        b0datasetBand = j
        
        # add Sigma0 and Beta0 bands with metadata
        # and corresponding values to the empty VRT
        self._create_bands(metaDict)

        #This doesn't work - the resulting array is still complex, but with
        #real part equal to its absolute value and the imaginary part equal to
        #zero....
        ##use ModulePixelFunc to calculate absolute sigma0 if data is complex!
        #metaDict = []
        #if s0dataset.GetRasterBand(1).DataType==10:
        #    # make new bands
        #    for i in range(1, s0dataset.RasterCount+1):
        #        polString = s0dataset.GetRasterBand(i).GetMetadata()['POLARIMETRIC_INTERP']
        #        suffix = polString
        #        metaDict.append(
        #                {'src': {'SourceFilename': self.fileName,
        #                    'SourceBand': i,
        #                    'DataType': s0dataset.GetRasterBand(i).DataType},
        #                'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', 
        #                    'suffix': suffix,
        #                    'PixelFunctionType': 'ModulePixelFunc',
        #                    #'dataType': 6, # src and dst datatypes must be the same and complex
        #                    'polarization': polString}})

        #self._create_bands(metaDict)


        ##############################################################
        # Adding derived band (incidence angle) calculated
        # using pixel function "BetaSigmaToIncidence":
        #      incidence = arcsin(sigma0/beta0)*180/pi
        ##############################################################

        # bruk enten sigma0 fra pixelfunksjon eller en omregning (ny
        # pikselfunksjon) med komplekse tall fra originalfil  - isafall blir
        # inci og kompleks, men det trenger ikke bety noe sa lenge de andre
        # bandene ogsa er komplekse...

        src = [{'SourceFilename': b0datasetName,
                'SourceBand':  b0datasetBand,
                'DataType': 6},
                {'SourceFilename': s0datasetName,
                'SourceBand': 1, 
                'DataType': 1}]
        dst = {'wkv': 'angle_of_incidence',
                    'PixelFunctionType': 'BetaSigmaToIncidence',
                    'dataType': 6,
                    'name': 'incidence_angle'}
            
        self._create_band(src,dst)
        self.dataset.FlushCache()

        ###################################################################
        # Add sigma0_VV - pixel function of sigma0_HH and beta0_HH
        # incidence angle is calculated within pixel function
        # It is assummed that HH is the first band in sigma0 and beta0 sub datasets
        ###################################################################
        if 'VV' not in pol and 'HH' in pol:        
            s0datasetNameHH = pol.index('HH')+1
            src = [{'SourceFilename': s0datasetName,
                    'SourceBand':  s0datasetNameHH,
                    'DataType': 6},
                   {'SourceFilename': b0datasetName,
                    'SourceBand':  b0datasetBand,
                    'DataType': 6}]
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sigma0HHBetaToSigma0VV',
                   'polarisation': 'VV',
                   'suffix': 'VV'}
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
