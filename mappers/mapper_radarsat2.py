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
                inciCalcIndex_SIGMA0 = 1
                inciCalcPol_SIGMA0 = s0dataset.GetRasterBand(inciCalcIndex_SIGMA0).GetMetadata()['POLARIMETRIC_INTERP']
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
                inciCalcIndex_BETA0 = 0
                b0dataset = gdal.Open(dataset[0])
                for j in range(1, b0dataset.RasterCount+1):
                    polString = b0dataset.GetRasterBand(j).GetMetadata()['POLARIMETRIC_INTERP']
                    #pdb.set_trace()
                    if polString==inciCalcPol_SIGMA0:
                        inciCalcIndex_BETA0 = i+1
                        break
                if inciCalcIndex_BETA0 == 0:
                    # Throw an error...
                    print 'This should be an error...'
                #if b0dataset.GetRasterBand(i).DataType==10:
                #    suffix = polString+'_complex'
                #else:
                suffix = polString
                #pol.append(polString)
                metaDict.append(
                        {'src': {
                                'SourceFilename': 'RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml',
                                'SourceBand': j,
                                'DataType': b0dataset.GetRasterBand(j).DataType
                            },
                            'dst': {
                                'wkv': 'radar_brightness_coefficient', 
                                'suffix': suffix,
                                'polarization': polString
                            }})
        
        #for i in range(1, gdalDataset.RasterCount+1):
        #    polString = gdalDataset.GetRasterBand(i).GetMetadata()['POLARIMETRIC_INTERP']
        #    pol.append(polString)
        #    metaDict.append(
        #        {'src': {'SourceFilename':
        #                 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml',
        #                 'SourceBand': i,
        #                 'DataType': 6}, # don't want hardcoding, so wrote the
        #                 # above and now it works with complex data as well
        #         'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', 
        #                 'suffix': polString,
        #                 'polarization': polString}})

        # add bands with metadata and corresponding values to the empty VRT
        #pdb.set_trace()
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
        # we check the datatype within the pixel function
        ##############################################################
        #pdb.set_trace()

        # bruk enten sigma0 fra pixelfunksjon eller en omregning (ny
        # pikselfunksjon) med komplekse tall fra originalfil  - isafall blir
        # inci og kompleks, men det trenger ikke bety noe sa lenge de andre
        # bandene ogsa er komplekse...

        #print 'Source band beta0: '+str(inciCalcIndex_BETA0)
        #print 'Source band sigma0: '+str(inciCalcIndex_SIGMA0)
        src = [{'SourceFilename': self.fileName,
                'SourceBand':  inciCalcIndex_BETA0,
                'DataType': b0dataset.GetRasterBand(j).DataType},
                {'SourceFilename': self.fileName,
                'SourceBand': inciCalcIndex_SIGMA0, 
                'DataType': s0dataset.GetRasterBand(1).DataType}
                ]
        dst = {'wkv': 'angle_of_incidence',
                    'PixelFunctionType': 'BetaSigmaToIncidence',
                    'dataType': s0dataset.GetRasterBand(inciCalcIndex_SIGMA0).DataType,
                    'name': 'incidence_angle'}
            
        self._create_band(src,dst)
        self.dataset.FlushCache()


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
