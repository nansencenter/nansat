#-------------------------------------------------------------------------------
# Name:    nansat
# Purpose: main of nansat module.
#          Reference nansat_open, nansat_transform and nansat_write
#
# Author:      asumak
#
# Created:     29.06.2011
# Copyright:   (c) asumak 2011
# Licence:
#-------------------------------------------------------------------------------
import os
import sys
from string import Template
from xml.etree.ElementTree import *

try:
    from osgeo import gdal
    from gdal import GDT_Float32
except ImportError:
    import gdal

class VRT():
    '''    Top class of Nansat readers    '''
    
    def __init__(self, metadata, rawVRTName):
        ''' Set attributes common for all mappers'''
        self.metadata = metadata;
        self.rawVRTName = rawVRTName;

    def createVRT_and_add_bands(self, ds, metaDict, vrtBandList):
        ''' Create VSI VRT dataset and set geo-metadata from source file '''
        self.vsiDs, srcRasterXSize, srcRasterYSize = self.createVRT(ds, len(vrtBandList));
        self.ds = ds # To be able to access bands of dataset, to know Blocksize and DataType
        #add bands with metadata and corresponding values
        self.vsiDs = self.addAllBands(self.vsiDs, vrtBandList, metaDict, srcRasterXSize, srcRasterYSize);

    def createVRT(self, srcDs, length):
        '''    Create VRT with common parameters '''
        #get geo-metadata from source dataset (or subdataset)
        srcGeoTransform = srcDs.GetGeoTransform();
        srcProjection = srcDs.GetProjection();
        srcProjectionRef = srcDs.GetProjectionRef();
        srcGCPCount = srcDs.GetGCPCount();
        srcGCPs = srcDs.GetGCPs();
        srcGCPProjection = srcDs.GetGCPProjection();

        srcRasterXSize = srcDs.GetRasterBand(1).XSize;
        srcRasterYSize = srcDs.GetRasterBand(1).YSize;

        #create VSI VRT dataset
        vrtDrv = gdal.GetDriverByName("VRT");
        vsiDs = vrtDrv.Create(self.rawVRTName,\
                              srcRasterXSize, srcRasterYSize,\
                              length, gdal.GDT_Float32 );

        # set geo-metadata in the VSI VRT dataset
        vsiDs.SetGCPs(srcGCPs, srcGCPProjection);
        vsiDs.SetProjection(srcProjection);
        vsiDs.SetGeoTransform(srcGeoTransform);
        
        #set metadata
        vsiDs.SetMetadata(self.metadata)

        return vsiDs, srcRasterXSize, srcRasterYSize

    def addAllBands(self, vsiDs, vrtBandList, metaDict, srcRasterXSize, srcRasterYSize):
            '''Loop through all bands and add metadata and band XML source'''
            # for bandNo in vrtBandList:
            for iBand in range(len(vrtBandList)):
                bandNo = vrtBandList[iBand];
                # check if
                if int(bandNo) > int(metaDict.__len__()):
                    print "vrt.addAllBands(): an element in the bandList is unproper";
                    return;
                # add metadata
                xBlockSize, yBlockSize = self.ds.GetRasterBand(bandNo).GetBlockSize()
                srcDataType = self.ds.GetRasterBand(bandNo).DataType
                wkv_name = metaDict[bandNo-1]["wkv"];
                wkvDict = self.get_wkv(wkv_name);
                for key in wkvDict:
                    vsiDs.GetRasterBand(iBand + 1).SetMetadataItem(key, wkvDict[key]);
                if "parameters" in metaDict[bandNo - 1]:
                    for metaName in metaDict[bandNo - 1]["parameters"]: # Warning: zero-indexing for medaDict, but 1-indexing for band numbers!
                        vsiDs.GetRasterBand(iBand + 1).SetMetadataItem(metaName, metaDict[bandNo-1]["parameters"][metaName]);
                # add band
                bandSource = self.SimpleSource.substitute(XSize=srcRasterXSize, \
                                YSize=srcRasterYSize, Dataset=metaDict[bandNo-1]['source'], \
                                SourceBand=metaDict[bandNo-1]['sourceBand'], \
                                BlockXSize=xBlockSize, BlockYSize=yBlockSize, DataType=srcDataType)
                vsiDs.GetRasterBand(iBand +1).SetMetadataItem("source_0", bandSource, "new_vrt_sources");
            vsiDs.FlushCache();
            return vsiDs;

    def get_wkv(self, wkv_name):
        ''' Get wkv from wkv.xml '''
        # fetch band information corresponding to the fileType
        fileName_wkv = os.path.join(os.path.dirname(os.path.realpath( __file__ )), "wkv.xml");
        fd = file(fileName_wkv, "rb");
        element = ElementTree(file=fd).getroot();

        for e1 in list(element):
            if e1.find("name").text == wkv_name:
                wkvDict = {"name" : wkv_name}
                for e2 in list(e1):
                    wkvDict[e2.tag] = e2.text;

        return wkvDict;
    
    def addPixelFunction(self, pixelFunction, bands, dataset, parameters):
        ''' Generic function for mappers to add PixelFunctions from bands in the same dataset 
        
        Warning: so far input parameter dataset refers to the original file on disk 
        and not the Nansat dataset. Correspondingly, bands refers to bands of the original gdal dataset. 
        This will however soon be changed, as it is more convenient to refer to band numbers
        of the current Nansat object. Then dataset-name may also be omitted '''
        newBand = self.vsiDs.GetRasterBand(self.vsiDs.RasterCount)
        options = ['subClass=VRTDerivedRasterBand', 'PixelFunctionType=' + pixelFunction]
        self.vsiDs.AddBand(datatype=GDT_Float32, options=options)
        md = {}
        for i, bandNo in enumerate(bands):
            BlockXSize, BlockYSize = self.ds.GetRasterBand(bandNo).GetBlockSize()
            DataType = self.ds.GetRasterBand(bandNo).DataType
            md['source_' + str(i)] = self.SimpleSource.substitute(XSize=self.vsiDs.RasterXSize, \
                                BlockXSize=BlockXSize, BlockYSize=BlockYSize, DataType=DataType,
                                YSize=self.vsiDs.RasterYSize, Dataset=dataset, SourceBand=bandNo)
        
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadata(md, 'vrt_sources')
    
        for parameter, value in parameters.items():
            self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadataItem(parameter, value)
        self.vsiDs.GetRasterBand(self.vsiDs.RasterCount).SetMetadataItem('pixelfunction', pixelFunction)
        self.vsiDs.FlushCache() # Took 5 hours of debugging to find this one!!!
               
    SimpleSource = Template('''
            <SimpleSource>
                <SourceFilename relativeToVRT="0">$Dataset</SourceFilename>
                <SourceBand>$SourceBand</SourceBand> \
                <SourceProperties RasterXSize="$XSize" RasterYSize="$YSize" \
                        DataType="$DataType" BlockXSize="$BlockXSize" BlockYSize="$BlockYSize"/>\
                <SrcRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
                <DstRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
            </SimpleSource> ''')
