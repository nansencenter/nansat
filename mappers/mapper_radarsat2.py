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

try:
    from osgeo import gdal
except ImportError:
    import gdal

from vrt import *

class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Radarsat2 '''

    def __init__(self, ds, fileName, metadata, vrtBandList, rawVRTName):
        ''' Create Radarsat2 VRT '''
        VRT.__init__(self, metadata, rawVRTName);

        if vrtBandList == None:
            vrtBandList = [1,2];

        product = metadata.get("SATELLITE_IDENTIFIER", "Not_RADARSAT-2")

        #if it is not RADARSAT-2, return
        if product!= 'RADARSAT-2':
            raise AttributeError("RADARSAT-2 BAD MAPPER");

        #define dictionary of metadata and band specific parameters
        # NB: it should be checked dynamically what is the polarizations of band 1 and band2, it should not be hardcoded as it is now
        metaDict = [ {'source': 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml', 'sourceBand': 1, 'wkv': 'sigma0', 'parameters': {'polarization': 'HH'}},\
                     {'source': 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml', 'sourceBand': 2, 'wkv': 'sigma0', 'parameters': {'polarization': 'HV'}} ];

        self.createVRT_and_add_bands(ds, metaDict, vrtBandList);

        ##############################################################
        # Adding derived band (incidence angle) calculated
        # using pixel function "BetaSigmaToIncidence":
        #      incidence = arcsin(sigma0/beta0)*180/pi 
        ##############################################################
        from gdal import GDT_Float32
        from string import Template
        
        options = ['subClass=VRTDerivedRasterBand', 'PixelFunctionType=BetaSigmaToIncidence']
        self.vsiDs.AddBand(datatype=GDT_Float32, options=options)  
        
        md = {}
        BlockXSize, BlockYSize = ds.GetRasterBand(1).GetBlockSize()
        md['source_0'] = self.SimpleSource.substitute(XSize=self.vsiDs.RasterXSize, YSize=self.vsiDs.RasterYSize, 
                                    BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=GDT_Float32, SourceBand=1,
                                    Dataset='RADARSAT_2_CALIB:BETA0:' + fileName + '/product.xml')
        md['source_1'] = self.SimpleSource.substitute(XSize=self.vsiDs.RasterXSize, YSize=self.vsiDs.RasterYSize, 
                                    BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=GDT_Float32, SourceBand=1,
                                    Dataset='RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml')
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadata(md, 'vrt_sources');
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetNoDataValue(0);
        # Should ideally use WKV-class for setting the metadata below
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadataItem('longname', 'incidence_angle');
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadataItem('name', 'incidence_angle');
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadataItem('unit', 'degrees');
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadataItem('pixelfunction', 'BetaSigmaToIncidence');
        self.vsiDs.FlushCache()
        
        return
