#-------------------------------------------------------------------------------
# Name:        nansat_mapper_radarsat2
# Purpose:     Mapping for Radarsat2 data
#
# Author:      asumak
#
# Created:     29.11.2011
# Copyright:   (c) asumak 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from datetime import datetime

try:
    from osgeo import gdal
except ImportError:
    import gdal

from vrt import *

class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Radarsat2 '''

    def __init__(self, rawVRTFileName, fileName, dataset, metadata, vrtBandList):        
        ''' Create Radarsat2 VRT '''
        VRT.__init__(self, dataset, metadata, rawVRTFileName);
        
        if vrtBandList == None:
            vrtBandList = [1,2];

        product = metadata.get("SATELLITE_IDENTIFIER", "Not_RADARSAT-2")

        #if it is not RADARSAT-2, return
        if product!= 'RADARSAT-2':
            raise AttributeError("RADARSAT-2 BAD MAPPER");

        #define dictionary of metadata and band specific parameters
        # NB: it should be checked dynamically what is the polarizations of band 1 and band2, it should not be hardcoded as it is now
        metaDict = [ {'source': 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml', 'sourceBand': 1, 'wkv': 'normalized_radar_cross_section', 'parameters': {'band_name': 'sigma0_HH', 'polarization': 'HH'}},\
                     {'source': 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml', 'sourceBand': 2, 'wkv': 'normalized_radar_cross_section', 'parameters': {'band_name': 'sigma0_HV', 'polarization': 'HV'}} ];

        self._createVRT(metaDict, vrtBandList);

        ##################################
        # Add time to metadata domain
        ##################################
        validTime = dataset.GetMetadata()['ACQUISITION_START_TIME']
        self.metadata['Valid Time'] = datetime.strptime(validTime, '%Y-%m-%dT%H:%M:%S.%fZ')

        ##############################################################
        # Adding derived band (incidence angle) calculated
        # using pixel function "BetaSigmaToIncidence":
        #      incidence = arcsin(sigma0/beta0)*180/pi 
        ##############################################################
        from gdal import GDT_Float32
        from string import Template
        
        options = ['subClass=VRTDerivedRasterBand', 'PixelFunctionType=BetaSigmaToIncidence']
        self.vsiDataset.AddBand(datatype=GDT_Float32, options=options)  
        
        md = {}
        BlockXSize, BlockYSize = dataset.GetRasterBand(1).GetBlockSize()
        md['source_0'] = self.SimpleSource.substitute(XSize=self.vsiDataset.RasterXSize, YSize=self.vsiDataset.RasterYSize, 
                                    BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=GDT_Float32, SourceBand=1,
                                    Dataset='RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml')
        md['source_1'] = self.SimpleSource.substitute(XSize=self.vsiDataset.RasterXSize, YSize=self.vsiDataset.RasterYSize, 
                                    BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=GDT_Float32, SourceBand=1,
                                    Dataset='RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml')
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadata(md, 'vrt_sources');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetNoDataValue(0);
        # Should ideally use WKV-class for setting the metadata below
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('long_name', 'incidence_angle');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('standard_name', 'incidence_angle');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('band_name', 'incidence_angle');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('unit', 'degrees');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('pixelfunction', 'BetaSigmaToIncidence');
        self.vsiDataset.FlushCache()
        
        # Experimental feature for the Radarsat2-mapper:
        # Rationale:
        # - to convert sigma0 from HH-pol to VV-pol we can use the pixelfunction 
        #     Sigma0HHIncidenceToSigma0VV which takes as input sigma0HH and incidence_angle.
        #     However, incidence_angle is itself a pixelfunction, so here we need a pixelfunction
        #     of a pixelfunction!
        # Issue:
        # - this second pixelfunction cannot access the first pixelfunction though the regular VRT/VSI-file
        #     because if we e.g. downscale the nansat object, then we end up with a double downscaling 
        #     of the second pixelfunction band
        # Solution:
        # - to duplicate the vsiDataset in another vsimem-file to remain unchanged (under reprojection and downscaling)
        vrtDriver = gdal.GetDriverByName("VRT")
        # Write the vrt to a VSI-file
        vrtDatasetCopy_temp = vrtDriver.CreateCopy('/vsimem/vsi_original.vrt',
                                                   self.vsiDataset)
        
        options = ['subClass=VRTDerivedRasterBand', 'PixelFunctionType=Sigma0HHIncidenceToSigma0VV']
        self.vsiDataset.AddBand(datatype=GDT_Float32, options=options)  
        md = {}
        BlockXSize, BlockYSize = dataset.GetRasterBand(1).GetBlockSize()
        md['source_0'] = self.SimpleSource.substitute(XSize=self.vsiDataset.RasterXSize, YSize=self.vsiDataset.RasterYSize, 
                                    BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=GDT_Float32, SourceBand=1,
                                    Dataset='RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml')
        md['source_1'] = self.SimpleSource.substitute(XSize=self.vsiDataset.RasterXSize, YSize=self.vsiDataset.RasterYSize, 
                                    BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=GDT_Float32, SourceBand=3,
                                    Dataset='/vsimem/vsi_original.vrt')
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadata(md, 'vrt_sources');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetNoDataValue(0);
        # Should ideally use WKV-class for setting the metadata below
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('long_name', 'normalized_radar_cross_section');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('standard_name', 'normalized_radar_cross_section');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('band_name', 'sigma0_VV');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('polarisation', 'VV');
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).SetMetadataItem('pixelfunction', 'Sigma0HHIncidenceToSigma0VV');
        self.vsiDataset.FlushCache()
                  
        return
