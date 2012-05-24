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
import datetime
import dateutil.parser

import logging

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal
    import osr

from nansat_tools import add_logger, Node, openXML


class VRT():
    '''VRT dataset management

    Used in Domain and Nansat
    Perfroms all peration on VRT datasets: creation, copying, modification,
    writing, etc.
    All mapper inherit from VRT
    '''
    DatasetSource = Template('''
            <VRTDataset rasterXSize="$XSize" rasterYSize="$YSize">
              <SRS>$SRS</SRS>
              <GeoTransform>$GeoTransform</GeoTransform>
            </VRTDataset> ''')

    RasterBandSource = Template('''
            <VRTRasterBand dataType="$DataType" band="$BandNum">
            </VRTRasterBand> ''')

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

    RawRasterBandSource = Template('''
            <VRTRasterBand dataType="$DataType" band="$BandNum" subClass="VRTRawRasterBand">
              <SourceFilename relativeToVRT="0">$SrcFileName</SourceFilename>
              <ImageOffset>0</ImageOffset>
              <PixelOffset>$PixelOffset</PixelOffset>
              <LineOffset>$LineOffset</LineOffset>
            </VRTRasterBand> ''')

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
            self.dataset = self.vrtDriver.CreateCopy(self.fileName,
                                                     vrtDataset)
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
                srcGCPs = []
                srcGCPProjection = None

            # create VRT dataset
            self.dataset = self.vrtDriver.Create(self.fileName,
                                                 srcRasterXSize,
                                                 srcRasterYSize,
                                                 bands=0)

            # set geo-metadata in the VRT dataset
            self.dataset.SetGCPs(srcGCPs, srcGCPProjection)
            self.dataset.SetProjection(srcProjection)
            self.dataset.SetGeoTransform(srcGeoTransform)

            # set metadata
            self.dataset.SetMetadata(srcMetadata)
            # write file contents
            self.dataset.FlushCache()

        self.logger.debug('VRT self.dataset: %s' % self.dataset)
        self.logger.debug('VRT description: %s ' %
                                             self.dataset.GetDescription())
        self.logger.debug('VRT metadata: %s ' % self.dataset.GetMetadata())
        self.logger.debug('VRT RasterXSize %d' % self.dataset.RasterXSize)
        self.logger.debug('VRT RasterYSize %d' % self.dataset.RasterYSize)

    def _make_filename(self):
        '''Create random VSI file name'''
        allChars = ascii_uppercase + digits
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
            source = [source] * len(sourceBands)

        # Find datatype and blocksizes
        srcDataset = gdal.Open(source[0])
        srcRasterBand = srcDataset.GetRasterBand(sourceBands[0])
        blockXSize, blockYSize = srcRasterBand.GetBlockSize()
        dataType = srcRasterBand.DataType

        # Create band
        if "pixel_function" in parameters:
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

            if "pixel_function" in parameters:
                md['source_' + str(i)] = bandSource

        # Append sources to destination dataset
        if "pixel_function" in parameters:
            dstRasterBand.SetMetadata(md, 'vrt_sources')
        else:
            dstRasterBand.SetMetadataItem("source_0", bandSource,
                                "new_vrt_sources")

        # set metadata from WKV and from provided parameters
        dstRasterBand = self._put_metadata(dstRasterBand, self._get_wkv(wkv))
        dstRasterBand = self._put_metadata(dstRasterBand, parameters)

        self.dataset.FlushCache()

        return

    def _set_time(self, time):
        ''' Set time of dataset and/or its bands

        Parameters
        ----------
        time: datetime

        If a single datetime is given, this is stored in
        all bands of the dataset as a metadata item "time".
        If a list of datetime objects is given, different
        time can be given to each band.

        '''
        # Make sure time is a list with one datetime element per band
        numBands = self.dataset.RasterCount
        if isinstance(time, datetime.datetime):
            time = [time]
        if len(time) == 1:
            time = time * numBands
        if len(time) != numBands:
            self.logger.error("Dataset has " + str(numBands) +
                    " elements, but given time has "
                    + str(len(time)) + " elements.")

        # Store time as metadata key "time" in each band
        for i in range(numBands):
            self.dataset.GetRasterBand(i + 1).SetMetadataItem(
                    'time', str(time[i].isoformat(" ")))

        return

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
        node0 = Node.create(openXML(fileName_wkv))
        for iNode in node0.nodeList("wkv"):
            tagsList = iNode.tagList()
            if iNode.node("standard_name").value == wkvName:
                wkvDict = {"standard_name": wkvName}
                for itag in tagsList:
                    wkvDict[itag] = str(iWkv.node(itag).value)

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

    def create_band_from_array(self, array):
        '''Create a band from an array

        Parameters:
        -----------
            array: numpy arrayvrt

        Modifies:
        ---------
            a new band which is created from the array is added to vrt.
        '''

        """ HOW to set file name???"""
        binaryTestFile = "c:/Users/asumak/Data/output/binaryTest.raw"

        ofile = open(binaryTestFile, 'wb')
        ofile.write(array.tostring())
        ofile.close()

        dataType = {"uint8": "Byte", "int8": "Byte",
                    "uint16": "UInt16", "int16": "Int16",
                    "uint32": "UInt32", "int32": "Int32",
                    "float32": "Float32","float64": "Float64",
                    "complex64": "CFloat64"}.get(str(array.dtype))

        pixelOffset = {"Byte": "1",
                    "UInt16": "2", "Int16": "2",
                    "UInt32": "4", "Int32": "4",
                    "Float32": "4","Float64": "8",
                    "CFloat64": "8"}.get(dataType)

        lineOffset = str(int(pixelOffset)*array.shape[1])

        datasetContents = self.read_xml()
        rasterBandSource = self.RawRasterBandSource.substitute(
                                        DataType=dataType,
                                        BandNum=1,
                                        SrcFileName=binaryTestFile,
                                        PixelOffset=pixelOffset,
                                        LineOffset=lineOffset)

        contents = Node.create(datasetContents).insert(rasterBandSource)

        self.write_xml(contents)

    def merge_vrt(self, vrt2):
        '''Merge two vrts

        If band size and projections are same, the bands are merged.
        If vrt has 'subClass' or 'simpleSource', they are copied.
        Otherwise 'simpleSource' is created and added for each band.

        Parameters:
        -----------
            vrt2: vrt

        Modifies:
        ---------
            self.vrt has "VRTDataset", "SRS", "GeoTransform", and each band info.
            it does not include global metadata.

        '''
        vrt1 = self.copy()

        projection1 = vrt1.dataset.GetProjection()
        if projection1 =="":
            projection1 = vrt1.dataset.GetGCPProjection()

        projection2 = vrt2.dataset.GetProjection()
        if projection2 =="":
            projection2 = vrt2.dataset.GetGCPProjection()
        '''
        print vrt1.dataset.RasterXSize ," : ",vrt2.dataset.RasterXSize
        print vrt1.dataset.RasterYSize ," : ",vrt2.dataset.RasterYSize
        print geotrans1," : ",geotrans2
        print projection1
        print projection2
        '''
        if vrt1.dataset.RasterXSize==vrt2.dataset.RasterXSize and \
           vrt1.dataset.RasterYSize==vrt2.dataset.RasterYSize and \
           projection1 == projection2 : ## and vrt1.dataset.GetGeoTransform() == vrt2.dataset.GetGeoTransform():

            datasetContents = self.DatasetSource.substitute(
                                    XSize=vrt1.dataset.RasterXSize,
                                    YSize=vrt1.dataset.RasterYSize,
                                    SRS= projection1,
                                    GeoTransform= str(vrt1.dataset.GetGeoTransform()).strip("()"))

            for iVrt in range(2):
                if iVrt == 0:
                    tmpVrt = vrt1
                else:
                    tmpVrt = vrt2

                for iBandNum, iBand in enumerate (Node.create(tmpVrt.read_xml()).nodeList('VRTRasterBand')):
                    nameList, valList = iBand.getAttributeList()
                    if 'subClass' in nameList and valList[nameList.index('subClass')] == "VRTRawRasterBand":
                        iBand.replaceAttribute("band", str(iVrt*vrt1.dataset.RasterCount+iBandNum+1))
                        bandContents = iBand.rawxml()
                    else:
                        if "SimpleSource" in iBand.tagList():
                            iBand.replaceAttribute("band", str(iVrt*vrt1.dataset.RasterCount+iBandNum+1))
                            bandContents = iBand.rawxml()
                        else:
                            band = tmpVrt.dataset.GetRasterBand(iBandNum + 1)
                            dataType = {1: "Byte",
                                        2: "UInt16", 3: "Int16",
                                        4: "UInt32", 5: "Int32",
                                        6: "Float32",7: "Float64",
                                        8: "CFloat64"}.get(band.DataType)

                            rasterBandSource = self.RasterBandSource.substitute(
                                               DataType=dataType,
                                               BandNum= iVrt*vrt1.dataset.RasterCount+iBandNum+1)

                            simpleSource = self.SimpleSource.substitute(
                                                XSize=band.XSize,
                                                YSize=band.YSize,
                                                BlockXSize=band.GetBlockSize()[0],
                                                BlockYSize=band.GetBlockSize()[1],
                                                DataType=dataType,
                                                Dataset=tmpVrt.fileName,
                                                SourceBand=band.GetBand())

                            bandContents = Node.create(rasterBandSource).insert(simpleSource)
                    datasetContents = Node.create(datasetContents).insert(bandContents)
            self.write_xml(datasetContents)
        else:
            self.logger.warning("Unable to merge vrt files")

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
