# Name:    vrt.py
# Purpose: Top class of Nansat mappers
#
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad
#
# Created:     29.06.2011
# Copyright:   (c) NERSC 2012
# Licence:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details:
# http://www.gnu.org/licenses/

import os
from string import Template, ascii_uppercase, digits
from random import choice

import logging

from xml.etree.ElementTree import ElementTree

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal
    import osr

from nansat_tools import add_logger

class VRT():
    '''VRT dataset management

    Used in Domain and Nansat
    Perfroms all peration on VRT datasets: creation, copying, modification,
    writing, etc.
    All mapper inherit from VRT
    '''

    def __init__(self, gdalDataset=None, vrtDataset=None,
                                         srcGeoTransform=None,
                                         srcProjection=None,
                                         srcRasterXSize=None,
                                         srcRasterYSize=None,
                                         srcMetadata=None,
                                         logLevel=30):
        ''' Create VRT dataset from GDAL dataset, or from given parameters
        
        If vrtDataset is given, creates full copy of VRT content
        Otherwise takes reprojection parameters (geotransform, projection, etc)
        either from given GDAL dataset or from seperate parameters.
        Create VRT dataset (self.dataset) based on these parameters
        Adds logger
        
        Parameters
        ----------
            gdalDataset: GDAL Dataset
                source dataset of geo-reference
            vrtDataset: GDAL VRT Dataset
                source dataset of all content (geo-reference and bands)
            srcGeoTransform: GDALGeoTransform
                parameter of geo-reference
            srcProjection, GDALProjection
                parameter of geo-reference
            srcRasterXSize, INT
                parameter of geo-reference
            srcRasterYSize, INT
                parameter of geo-reference
            srcMetadata: GDAL Metadata
                additional parameter
            logLevel: int
        
        Modifies:
        ---------
            self.dataset: GDAL VRT dataset
            self.logger: logging logger
            self.vrtDriver: GDAL Driver
        
        '''
        # essential attributes
        self.logger = add_logger('Nansat', logLevel=logLevel)
        self.fileName = self._make_filename()
        self.vrtDriver = gdal.GetDriverByName("VRT")
        self.logger.debug('input vrtDataset: %s' % str(vrtDataset))
        # copy content of the provided VRT dataset
        if vrtDataset is not None:
            self.logger.debug('Making copy of %s ' % str(vrtDataset))
            self.dataset = self.vrtDriver.CreateCopy(self.fileName, vrtDataset)
        else:
            # get geo-metadata from given GDAL dataset
            if gdalDataset is not None:
                srcGeoTransform = gdalDataset.GetGeoTransform()
                srcProjection = gdalDataset.GetProjection()
                srcProjectionRef = gdalDataset.GetProjectionRef()
                srcGCPCount = gdalDataset.GetGCPCount()
                srcGCPs = gdalDataset.GetGCPs()
                srcGCPProjection = gdalDataset.GetGCPProjection()
        
                srcRasterXSize = gdalDataset.RasterXSize
                srcRasterYSize = gdalDataset.RasterYSize
                
                srcMetadata = gdalDataset.GetMetadata()
                self.logger.debug('RasterXSize %d' % gdalDataset.RasterXSize)
                self.logger.debug('RasterYSize %d' % gdalDataset.RasterYSize)
            else:
                srcGCPs=[]
                srcGCPProjection=None
    
            # create VRT dataset
            self.dataset = self.vrtDriver.Create(self.fileName,
                                            srcRasterXSize, srcRasterYSize, bands=0)
    
            # set geo-metadata in the VRT dataset
            self.dataset.SetGCPs(srcGCPs, srcGCPProjection)
            self.dataset.SetProjection(srcProjection)
            self.dataset.SetGeoTransform(srcGeoTransform)
    
            # set metadata
            self.dataset.SetMetadata(srcMetadata)
            # write file contents
            self.dataset.FlushCache()

        self.logger.debug('VRT self.dataset: %s' % self.dataset)
        self.logger.debug('VRT description: %s ' % self.dataset.GetDescription())
        self.logger.debug('VRT metadata: %s ' % self.dataset.GetMetadata())
        self.logger.debug('VRT RasterXSize %d' % self.dataset.RasterXSize)
        self.logger.debug('VRT RasterYSize %d' % self.dataset.RasterYSize)        
    
    def _make_filename(self):
        '''Create random VSI file name'''
        allChars=ascii_uppercase + digits
        randomChars = ''.join(choice(allChars) for x in range(10))
        
        return '/vsimem/%s.vrt' % randomChars

    def _create_bands(self, metaDict):
        ''' Generic function called from the mappers to create bands
        in the VRT dataset from an input dictionary of metadata
        
        Keys and values in the metaDict dictionary:
        -------------------------------------------
        source: string
            name of the source dataset (e.g. filename)
        sourceBand: integer
            band number of the source band in the given source dataset
        wkv: string
            refers to the "standard_name" of some of the 
            "well known variables" listed in wkv.xml
            The corresponding parameters are added as metadata to the band
        parameters: dictionary
            metadata to be added to the band: {key: value}
        
        If one of the latter parameter keys is "pixel_function", this 
        band will be a pixel function defined by the corresponding name/value.
        In this case sourceBands and source may be lists of integers/strings 
        (i.e. possibly several bands and several sources as input).
        If source is a single string, this source will 
        be used for all source bands.
         
        '''
        for bandDict in metaDict:
            self._create_band(bandDict["source"], bandDict["sourceBand"],
                              bandDict["wkv"], bandDict["parameters"])
        return

    def _create_band(self, source, sourceBands, wkv, parameters):
        ''' Function to add a band to the VRT from a source.
        See function _create_bands() for explanation of the input parameters
        
        ''' 
        # Make sure sourceBands and source are lists, ready for loop
        # There will be a single sourceBand for regular bands, 
        # but several for bands which are pixel functions 
        if isinstance(sourceBands, int):
            sourceBands = [sourceBands]
        if isinstance(source, str):
            source = [source]*len(sourceBands)

        # Find datatype and blocksizes
        srcDataset = gdal.Open(source[0])
        srcRasterBand = srcDataset.GetRasterBand(sourceBands[0])
        blockXSize, blockYSize = srcRasterBand.GetBlockSize()
        dataType = srcRasterBand.DataType
        
        # Create band
        if parameters.has_key("pixel_function"):
            options = ['subClass=VRTDerivedRasterBand',
                       'PixelFunctionType=' + parameters["pixel_function"]]
        else:
            options = []
        self.dataset.AddBand(dataType, options=options)
        dstRasterBand = self.dataset.GetRasterBand(
                                        self.dataset.RasterCount)
        
        # Prepare sources 
        # (only one item for regular bands, several for pixelfunctions) 
        md = {}
        for i in range(len(sourceBands)):
            bandSource = self.SimpleSource.substitute(
                                XSize=self.dataset.RasterXSize,
                                YSize=self.dataset.RasterYSize,
                                BlockXSize=blockXSize,
                                BlockYSize=blockYSize,
                                DataType=dataType,
                                Dataset=source[i], SourceBand=sourceBands[i])
            
            if parameters.has_key("pixel_function"):
                md['source_' + str(i)] = bandSource
        
        # Append sources to destination dataset
        if parameters.has_key("pixel_function"):
            dstRasterBand.SetMetadata(md, 'vrt_sources')
        else:
            dstRasterBand.SetMetadataItem("source_0", bandSource,
                                "new_vrt_sources")
 
        # set metadata from WKV and from provided parameters
        dstRasterBand = self._put_metadata(dstRasterBand, self._get_wkv(wkv))
        dstRasterBand = self._put_metadata(dstRasterBand, parameters)

        self.dataset.FlushCache()

        return
    
    def _add_pixel_function(self, pixelFunction, bands, fileName, metaDict):
        ''' Generic function for mappers to add PixelFunctions
        from bands in the same dataset

        Warning: so far input parameter dataset refers to the original
        file on disk and not the Nansat dataset. Correspondingly,
        bands refers to bands of the original gdal dataset.
        This will however soon be changed, as it is more convenient to refer
        to band numbers of the current Nansat object.
        Then dataset-name may also be omitted

        Parameters
        ----------
        pixelFunction: string
            value of 'pixelfunction' attribute for each band. Name of
            the actual pixel function.
        bands: list
            input band numbers
        fileName: string
            name of the file with input bands
        metaDict: dictionary
            metadata to be included into a band

        Modifies
        --------
        self.dataset: VRT dataset
            add PixelFunctions from bands in the same dataset
        '''

        newBand = self.dataset.GetRasterBand(self.dataset.RasterCount)
        options = ['subClass=VRTDerivedRasterBand',
                   'PixelFunctionType=' + pixelFunction]
        self.dataset.AddBand(datatype=gdal.GDT_Float32, options=options)
        md = {}
        srcDataset = gdal.Open(fileName)
        for i, bandNo in enumerate(bands):
            srcRasterBand = srcDataset.GetRasterBand(bandNo)
            blockXSize, blockYSize = srcRasterBand.GetBlockSize()
            dataType = srcRasterBand.DataType
            md['source_' + str(i)] = self.SimpleSource.substitute(
                                        XSize=self.dataset.RasterXSize,
                                        BlockXSize=blockXSize,
                                        BlockYSize=blockYSize,
                                        DataType=dataType,
                                        YSize=self.dataset.RasterYSize,
                                        Dataset=fileName, SourceBand=bandNo)

        # set metadata for each destination raster band
        dstRasterBand = self.dataset.GetRasterBand(self.dataset.\
                                                        RasterCount)

        dstRasterBand.SetMetadata(md, 'vrt_sources')

        # set metadata from WKV
        wkvName = metaDict["wkv"]
        dstRasterBand = self._put_metadata(dstRasterBand,
                                            self._get_wkv(wkvName))
        # set metadata from parameters (if exist)
        if "parameters" in metaDict:
            dstRasterBand = self._put_metadata(dstRasterBand,
                                               metaDict["parameters"])

        dstRasterBand.SetMetadataItem('pixelfunction', pixelFunction)
        # Took 5 hours of debugging to find this one!!!
        self.dataset.FlushCache()

    SimpleSource = Template('''
            <SimpleSource>
                <SourceFilename relativeToVRT="0">$Dataset</SourceFilename>
                <SourceBand>$SourceBand</SourceBand>
                <SourceProperties RasterXSize="$XSize" RasterYSize="$YSize"
                        DataType="$DataType" BlockXSize="$BlockXSize"
                        BlockYSize="$BlockYSize"/>
                <SrcRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
                <DstRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
            </SimpleSource> ''')

    def _add_all_bands(self, metaDict):
        '''Loop through all bands and add metadata and band XML source

        Parameters
        -----------
        metaDict: list
            incldes some dictionaries.
            The number of dictionaries is same as number of bands.
            Each dictionary represents metadata for each band.

        Modifies
        --------
        self.dataset: VRT dataset
            Band data and metadata is added to the VRT dataset
        '''
        self.logger.debug('self.dataset: %s' % str(self.dataset))
        for bandNo in range(metaDict.__len__()):
            srcRasterBand = gdal.Open(metaDict[bandNo]['source']).\
                       GetRasterBand(metaDict[bandNo]['sourceBand'])

            xBlockSize, yBlockSize = srcRasterBand.GetBlockSize()
            srcDataType = srcRasterBand.DataType
            bandMetaParameters = metaDict[bandNo].get('parameters', {})
            
            # add a band to the VRT dataset
            self.dataset.AddBand(srcDataType)
            # set metadata for each destination raster band
            dstRasterBand = self.dataset.GetRasterBand(bandNo+1)
            # set metadata from WKV
            wkvName = metaDict[bandNo]["wkv"]
            dstRasterBand = self._put_metadata(dstRasterBand,
                                                self._get_wkv(wkvName))
            # set metadata from 'parameters'
            dstRasterBand = self._put_metadata(dstRasterBand, bandMetaParameters)

            # create band source metadata
            bandSource = self.SimpleSource.substitute(
                    XSize=self.dataset.RasterXSize,
                    YSize=self.dataset.RasterYSize,
                    Dataset=metaDict[bandNo]['source'],
                    SourceBand=metaDict[bandNo]['sourceBand'],
                    BlockXSize=xBlockSize, BlockYSize=yBlockSize,
                    DataType=srcDataType)

            # set band source metadata
            dstRasterBand.SetMetadataItem("source_0", bandSource,
                                          "new_vrt_sources")
        self.dataset.FlushCache()

    def _get_wkv(self, wkvName):
        ''' Get wkv from wkv.xml

        Parameters
        ----------
        wkvName: string
            value of 'wkv' key in metaDict

        Returns
        -------
        wkvDict: dictionay
            WKV corresponds to the given wkv_name

        '''
        # fetch band information corresponding to the fileType
        fileName_wkv = os.path.join(os.path.dirname(
                                    os.path.realpath(__file__)), "wkv.xml")
        fd = file(fileName_wkv, "rb")
        element = ElementTree(file=fd).getroot()

        for e1 in list(element):
            if e1.find("standard_name").text == wkvName:
                wkvDict = {"standard_name": wkvName}
                for e2 in list(e1):
                    wkvDict[e2.tag] = e2.text

        return wkvDict

    def _put_metadata(self, rasterBand, metadataDict):
        ''' Put all metadata into a raster band

        Take metadata from metadataDict and put to the GDAL Raster Band

        Parameters:
        ----------
        rasterBand: GDALRasterBand
            destination band without metadata

        metadataDict: dictionary
            keys are names of metadata, values are values

        Returns:
        --------
        rasterBand: GDALRasterBand
            destination band with metadata
        '''
        self.logger.debug('Put: %s ' % str(metadataDict))
        for key in metadataDict:
            rasterBand.SetMetadataItem(key, metadataDict[key])

        return rasterBand

    def read_xml(self):
        '''Read XML content of the VRT dataset

        Returns:
            vsiFileContent: string
                XMl Content which is read from the VSI file
        '''

        # write dataset content into VRT-file
        self.dataset.FlushCache()
        #read from the vsi-file
        # open
        vsiFile = gdal.VSIFOpenL(self.fileName, "r")
        # get file size
        gdal.VSIFSeekL(vsiFile, 0, 2)
        vsiFileSize = gdal.VSIFTellL(vsiFile)
         # fseek to start again
        gdal.VSIFSeekL(vsiFile, 0, 0)
        # read
        vsiFileContent = gdal.VSIFReadL(vsiFileSize, 1, vsiFile)
        gdal.VSIFCloseL(vsiFile)
        return vsiFileContent

    def write_xml(self, vsiFileContent=None):
        '''Write XML content into a VRT dataset

        Parameters:
            vsiFileContent: string, optional
                XML Content of the VSI file to write
        Modifies:
            self.dataset
                If XML content was written, self.dataset is re-opened
        '''
        #write to the vsi-file
        vsiFile = gdal.VSIFOpenL(self.fileName, 'w')
        gdal.VSIFWriteL(vsiFileContent,
                        len(vsiFileContent), 1, vsiFile)
        gdal.VSIFCloseL(vsiFile)
        # re-open self.dataset with new content
        self.dataset = gdal.Open(self.fileName)
                
    def export(self, fileName):
        '''Export VRT file as XML into <fileName>'''
        self.vrtDriver.CreateCopy(fileName, self.dataset)

    def copy(self):
        '''Creates full copy of VRT dataset'''
        try:
            # deep copy (everything including bands)
            vrt = VRT(vrtDataset=self.dataset, logLevel=self.logger.level)
        except:
            # shallow copy (only geometadata)
            vrt = VRT(gdalDataset=self.dataset, logLevel=self.logger.level)
        
        return vrt
