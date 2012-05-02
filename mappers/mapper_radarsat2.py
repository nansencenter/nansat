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

from datetime import datetime
from numpy import mod

from vrt import *
from domain import Domain


class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Radarsat2 '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        ''' Create Radarsat2 VRT '''
        product = gdalMetadata.get("SATELLITE_IDENTIFIER", "Not_RADARSAT-2")
        #if it is not RADARSAT-2, return
        if product!= 'RADARSAT-2':
            raise AttributeError("RADARSAT-2 BAD MAPPER");


        #define dictionary of metadata and band specific parameters
        pol = []        
        metaDict = []

        for i in range(1, gdalDataset.RasterCount+1):
            polString = gdalDataset.GetRasterBand(i).GetMetadata()['POLARIMETRIC_INTERP'] 
            pol.append(polString)
            metaDict.append( {'source': 'RADARSAT_2_CALIB:SIGMA0:' + fileName +
                '/product.xml', 'sourceBand': i, 'wkv':
                'normalized_radar_cross_section', 'parameters':
                {'band_name':'sigma0_'+polString, 'polarization': polString}} );

        ##############################################################
        # Adding derived band (incidence angle) calculated
        # using pixel function "BetaSigmaToIncidence":
        #      incidence = arcsin(sigma0/beta0)*180/pi 
        ##############################################################
        incidence_angle_source = [
            'RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml',
            'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml']
        metaDict.append(
            {'source': incidence_angle_source, 'sourceBand': [1, 1], 'wkv': 'sensor_zenith_angle',
                'parameters': {'pixel_function': 'BetaSigmaToIncidence', 'band_name': 'incidence_angle'}})
        
        
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)


        ##################################
        # Add time to metadata domain
        ##################################
        validTime = gdalDataset.GetMetadata()['ACQUISITION_START_TIME']
        self.dataset.SetMetadataItem('time', str(datetime.strptime(validTime, '%Y-%m-%dT%H:%M:%S.%fZ')))

        ############################################
        # Add SAR look direction to metadata domain
        ############################################
        self.dataset.SetMetadataItem('SAR_look_direction', str(mod(
            Domain(gdalDataset, logLevel=self.logger.level).upwards_azimuth_direction()
            + 90, 360)))


        # Experimental feature for the Radarsat2-mapper:
        # Rationale:
        # - to convert sigma0 from HH-pol to VV-pol we can use the
        #   pixelfunction Sigma0HHIncidenceToSigma0VV which takes as input
        #   sigma0HH and incidence_angle. However, incidence_angle is itself a
        #   pixelfunction, so here we need a pixelfunction of a pixelfunction!
        # Issue:
        # - this second pixelfunction cannot access the first pixelfunction 
        #   through the regular VRT/VSI-file because if we e.g. downscale the 
        #   nansat object, then we end up with a double downscaling of the 
        #   second pixelfunction band
        # Solution:
        # - to copy the vsiDataset into another vsimem-file which remains
        #   unchanged under reprojection and downscaling
        if 'VV' not in pol and 'HH' in pol:
            # Write the vrt to a VSI-file
            vrtDatasetCopy_temp = self.vrtDriver.CreateCopy(
                    '/vsimem/vsi_original.vrt', self.dataset)
            options = ['subClass=VRTDerivedRasterBand', 
                    'PixelFunctionType=Sigma0HHIncidenceToSigma0VV']
            self.dataset.AddBand(datatype=gdal.GDT_Float32, options=options)
            md = {}
            for i in range(len(pol)):
                if pol[i]=='HH':
                    sourceBandHH = i+1
            sourceBandInci = len(pol)+1

            BlockXSize, BlockYSize = gdalDataset.GetRasterBand(1).GetBlockSize()
            md['source_0'] = self.SimpleSource.substitute(
                    XSize=self.dataset.RasterXSize,
                    YSize=self.dataset.RasterYSize, BlockXSize=BlockXSize, 
                    BlockYSize=BlockYSize, DataType=gdal.GDT_Float32,
                    SourceBand=sourceBandHH,
                    Dataset='RADARSAT_2_CALIB:BETA0:'+fileName+'/product.xml')
            md['source_1'] = self.SimpleSource.substitute(
                    XSize=self.dataset.RasterXSize,
                    YSize=self.dataset.RasterYSize, BlockXSize=BlockXSize,
                    BlockYSize=BlockYSize, DataType=gdal.GDT_Float32, 
                    SourceBand=sourceBandInci, 
                    Dataset='/vsimem/vsi_original.vrt')
            self.dataset.GetRasterBand(self.dataset.RasterCount).SetMetadata(md, 'vrt_sources')
            self.dataset.GetRasterBand(self.dataset.RasterCount).SetNoDataValue(0)
            
            self._put_metadata(self.dataset.GetRasterBand(self.dataset.RasterCount),
                           self._get_wkv('normalized_radar_cross_section'))
            self._put_metadata(self.dataset.GetRasterBand(self.dataset.RasterCount),
                {'band_name': 'sigma0_VV',
                 'polarisation': 'VV',
                'pixelfunction': 'Sigma0HHIncidenceToSigma0VV'})
            self.dataset.FlushCache()
                  
        return
