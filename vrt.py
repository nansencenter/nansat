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
from dateutil.parser import parse
import logging

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal, osr

from nansat_tools import add_logger, Node

class Geolocation():
    '''Container for GEOLOCATION data'''
    def __init__(self, xVRT=None, yVRT=None, xBand=1, yBand=1,
                        srs="",
                        lineOffset=0, lineStep=1, pixelOffset=0, pixelStep=1,
                        dataset=None):
        '''Create Geolocation object from input parameters
        Parameters:
            xVRT: VRT-object or str
                VRT with array of x-coordinates OR string with dataset source
            yVRT: VRT-object or str
                VRT with array of y-coordinates OR string with dataset source
            xBand: number of band in the xDataset
            xBand: number of band in the yDataset
            srs: str, WKT
            lineOffset: int, offset of first line
            lineStep: int, step of lines
            pixelOffset: int: offset of first pixel
            pixelStep: step of pixels
            dataset: GDAL dataset to take geolocation from
        Modifies:
            All input parameters are copied to self
        '''
        # dictionary with all metadata
        self.d = {}
        # VRT objects
        self.xVRT = None
        self.yVRT = None

        # make object from GDAL dataset
        if dataset is not None:
            self.d = dataset.GetMetadata('GEOLOCATION')
            return

        # make empty object
        if xVRT is None or yVRT is None:
            return

        if isinstance(xVRT, str):
            # make object from strings
            self.d['X_DATASET'] = xVRT
            self.d['Y_DATASET'] = yVRT
        else:
            # make object from VRTs
            self.xVRT = xVRT
            self.d['X_DATASET'] = xVRT.fileName
            self.yVRT = yVRT
            self.d['Y_DATASET'] = yVRT.fileName

        # proj4 to WKT
        if srs == "":
            sr = osr.SpatialReference()
            sr.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
            srs = sr.ExportToWkt()
        self.d['SRS'] = srs
        self.d['X_BAND'] = str(xBand)
        self.d['Y_BAND'] = str(yBand)
        self.d['LINE_OFFSET'] = str(lineOffset)
        self.d['LINE_STEP'] = str(lineStep)
        self.d['PIXEL_OFFSET'] = str(pixelOffset)
        self.d['PIXEL_STEP'] = str(pixelStep)



class VRT():
    '''VRT dataset management

    Used in Domain and Nansat
    Perfroms all peration on VRT datasets: creation, copying, modification,
    writing, etc.
    All mapper inherit from VRT
    '''
    ComplexSource = Template('''
            <$SourceType>
                <SourceFilename relativeToVRT="0">$Dataset</SourceFilename>
                <SourceBand>$SourceBand</SourceBand>
                <SourceProperties RasterXSize="$XSize" RasterYSize="$YSize"
                        DataType="$DataType" BlockXSize="$BlockXSize"
                        BlockYSize="$BlockYSize"/>
                <NODATA>$NODATA</NODATA>
                <LUT>$LUT</LUT>
                <SrcRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
                <DstRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
            </$SourceType> ''')

    RawRasterBandSource = Template('''
            <VRTDataset rasterXSize="$XSize" rasterYSize="$YSize">
              <VRTRasterBand dataType="$DataType" band="$BandNum" subClass="VRTRawRasterBand">
                <SourceFilename relativeToVRT="0">$SrcFileName</SourceFilename>
                <ImageOffset>0</ImageOffset>
                <PixelOffset>$PixelOffset</PixelOffset>
                <LineOffset>$LineOffset</LineOffset>
              </VRTRasterBand>
            </VRTDataset> ''')

    def __init__(self, gdalDataset=None, vrtDataset=None,
                                         array=None,
                                         srcGeoTransform=(0,1,0,0,0,1),
                                         srcProjection="",
                                         srcRasterXSize=None,
                                         srcRasterYSize=None,
                                         srcGCPs=[],
                                         srcGCPProjection="",
                                         srcMetadata="",
                                         geolocation=None):
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
            array: Numpy array
                source matrix with data
            srcGeoTransform: GDALGeoTransform
                parameter of geo-reference
            srcProjection, GDALProjection
                parameter of geo-reference
            srcRasterXSize, INT
                parameter of geo-reference
            srcRasterYSize, INT
                parameter of geo-reference
            srcMetadata: GDAL Metadata
                all global metadata
            geolocation: Geolocation
                object with info on geolocation and VRTs with x/y datasets

        Modifies:
        ---------
            self.dataset: GDAL VRT dataset
            self.logger: logging logger
            self.vrtDriver: GDAL Driver

        '''
        # essential attributes
        self.logger = add_logger('Nansat')
        self.fileName = self._make_filename()
        self.vrtDriver = gdal.GetDriverByName("VRT")
        if vrtDataset is not None:
            # copy content of the provided VRT dataset using CreateCopy
            self.logger.debug('copy content of the provided VRT dataset using CreateCopy')
            self.dataset = self.vrtDriver.CreateCopy(self.fileName,
                                                     vrtDataset)
            # add geolocation from vrt dataset
            self.add_geolocation(Geolocation(dataset=vrtDataset))
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

            # create VRT dataset (empty or with a band from array)
            if array is None:
                self.dataset = self.vrtDriver.Create(self.fileName,
                                                     srcRasterXSize,
                                                     srcRasterYSize,
                                                     bands=0)
            else:
                self.create_dataset_from_array(array)

            # set geo-metadata in the VRT dataset
            self.dataset.SetGCPs(srcGCPs, srcGCPProjection)
            self.dataset.SetProjection(srcProjection)
            self.dataset.SetGeoTransform(srcGeoTransform)

            # set metadata
            self.dataset.SetMetadata(srcMetadata)

            # add geolocation from gdal dataset
            self.add_geolocation(Geolocation(dataset=gdalDataset))

        # add geolocation from input geolocation
        # if not None: overwrite geoloc from vrt or gdal datasets
        if geolocation is not None:
            self.add_geolocation(geolocation)

        # add self.fileName to metadata
        self.dataset.SetMetadataItem('fileName', self.fileName)

        # write file contents
        self.dataset.FlushCache()

        self.logger.debug('VRT self.dataset: %s' % self.dataset)
        self.logger.debug('VRT description: %s ' %
                                             self.dataset.GetDescription())
        #self.logger.debug('VRT metadata: %s ' % self.dataset.GetMetadata())
        self.logger.debug('VRT RasterXSize %d' % self.dataset.RasterXSize)
        self.logger.debug('VRT RasterYSize %d' % self.dataset.RasterYSize)

    def __del__(self):
        ''' Destructor deletes VRT and RAW files'''
        try:
            gdal.Unlink(self.fileName)
            gdal.Unlink(self.fileName.replace('vrt', 'raw'))
        except:
            pass

    def _make_filename(self, extention="vrt"):
        '''Create random VSI file name'''
        allChars = ascii_uppercase + digits
        randomChars = ''.join(choice(allChars) for x in range(10))
        return '/vsimem/%s.%s' % (randomChars, extention)

    def _create_bands(self, metaDict):
        ''' Generic function called from the mappers to create bands
        in the VRT dataset from an input dictionary of metadata

        Keys and values in the metaDict dictionary:
        -------------------------------------------
        source: string
            name of the source dataset (e.g. filename)
        sourceBand: integer
            band number of the source band in the given source dataset
            if it is 0, it means VRTRawRasterband
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
            # get values or set default
            NODATA = bandDict.get("NODATA", "")
            LUT = bandDict.get("LUT", "")
            SourceType = bandDict.get('SourceType', 'ComplexSource')
            self.logger.debug('SourceType %s' % SourceType)
            self.logger.debug('Creating band %s', str(bandDict))
            self._create_band(bandDict["source"], bandDict["sourceBand"],
                    bandDict["wkv"], bandDict.get("parameters", {}), NODATA, LUT, SourceType)
            self.logger.debug('OK!')
        self.dataset.FlushCache()

    def _create_band(self, source, sourceBands, wkv, parameters, NODATA="", LUT="", SourceType='ComplexSource'):
        ''' Function to add a band to the VRT from a source.
        See function _create_bands() for explanation of the input parameters

        Returns:
        --------
            band_name: string, name of the added band

        '''
        self.logger.debug('INPUTS: %s, %s %s %s" ' % (str(source), str(sourceBands), str(wkv), str(parameters)))
        # Make sure sourceBands and source are lists, ready for loop
        # There will be a single sourceBand for regular bands,
        # but several for bands which are pixel functions
        if isinstance(sourceBands, int):
            sourceBands = [sourceBands]
        if isinstance(source, str):
            source = [source] * len(sourceBands)

        # add source and sourceBand parameters for Band metadata
        parameters["sourceBands"] = str(sourceBands[0])
        parameters["source"] = source[0]

        # create list of available bands (to prevent duplicate names)
        bandNames = []
        for iBand in range(self.dataset.RasterCount):
            bandNames.append(self.dataset.GetRasterBand(iBand+1).GetMetadataItem("band_name"))

        # if band_name is not given add 'band_00N'
        if "band_name" not in parameters:
            for n in range(999):
                bandName = 'band_%03d' % n
                if bandName not in bandNames:
                    parameters['band_name'] = bandName
                    break
        # if band_name already exist add '_00N'
        elif parameters["band_name"] in bandNames:
            for n in range(999):
                bandName = parameters["band_name"] + '_%03d' % n
                if bandName not in bandNames:
                    parameters['band_name'] = bandName
                    break

        # if add VRT RawRasterBand
        if sourceBands[0] == 0:
            options = ["subclass=VRTRawRasterBand",
                       "SourceFilename=%s" % parameters.pop("source"),
                       "ImageOffset=%f" % parameters.pop("ImageOffset"),
                       "PixelOffset=%f" % parameters.pop("PixelOffset"),
                       "LineOffset=%f" % parameters.pop("LineOffset"),
                       "ByteOrder=%s" % parameters.pop("ByteOrder")]
            dataType = parameters.pop("dataType")

        # else
        else:
            # Find datatype and blocksizes
            srcDataset = gdal.Open(source[0])
            srcRasterBand = srcDataset.GetRasterBand(sourceBands[0])
            blockXSize, blockYSize = srcRasterBand.GetBlockSize()
            if "datatype" in parameters:
                dataType = parameters.pop("dataType")
                dataType = {"uint8": gdal.GDT_Byte, "int8": gdal.GDT_Byte,
                    "uint16": gdal.GDT_UInt16, "int16":  gdal.GDT_Int16,
                    "uint32": gdal.GDT_UInt32, "int32":  gdal.GDT_Int32,
                    "float32": gdal.GDT_Float32,"float64": gdal.GDT_Float64,
                    "complex64":  gdal.GDT_CFloat64}.get(str(arrayDType), gdal.GDT_Float32)
            else:
                dataType = srcRasterBand.DataType

            dataType = 6
            # If we apply LUT, we must allow Byte values to be mapped into floats
            if LUT <> "" and dataType == gdal.GDT_Byte:
                dataType = gdal.GDT_Float32

            # Create band
            if "pixel_function" in parameters:
                options = ['subClass=VRTDerivedRasterBand',
                           'PixelFunctionType=' + parameters["pixel_function"]]
            else:
                options = []
        self.dataset.AddBand(dataType, options=options)
        dstRasterBand = self.dataset.GetRasterBand(self.dataset.RasterCount)

        # Prepare sources
        # (only one item for regular bands, several for pixelfunctions)
        if SourceType == "ComplexSource" or SourceType == "SimpleSource":
            md = {}
            rasterXSize=self.dataset.RasterXSize
            rasterYSize=self.dataset.RasterYSize
            for i in range(len(sourceBands)):
                bandSource = self.ComplexSource.substitute(
                                    SourceType=SourceType,
                                    XSize=rasterXSize,
                                    YSize=rasterYSize,
                                    BlockXSize=blockXSize,
                                    BlockYSize=blockYSize,
                                    NODATA=NODATA,
                                    LUT=LUT,
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

        # return name of the created band
        return parameters['band_name']

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
        node0 = Node.create(fileName_wkv)
        wkvDict = {}
        for iNode in node0.nodeList("wkv"):
            tagsList = iNode.tagList()
            if iNode.node("standard_name").value == wkvName:
                wkvDict = {"standard_name": wkvName}
                for iTag in tagsList:
                    wkvDict[iTag] = str(iNode.node(iTag).value)
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

    def create_dataset_from_array(self, array):
        '''Create a dataset with a band from an array

        Writes contents of the array into flat binary file (VSI)
        Writes VRT file with RawRastesrBand, which points to the binary file
        Opens the VRT file as self.dataset with GDAL

        Parameters:
        -----------
            array: numpy array

        Modifies:
        ---------
            binary file is written (VSI)
            VRT file is written (VSI)
            self.dataset is opened
        '''
        arrayDType = array.dtype.name
        arrayShape = array.shape
        # create flat binary file from array (in VSI)
        binaryFile = self.fileName.replace(".vrt", ".raw")
        ofile = gdal.VSIFOpenL(binaryFile, "wb")
        gdal.VSIFWriteL(array.tostring(), len(array.tostring()), 1, ofile)
        gdal.VSIFCloseL(ofile)
        array = None

        self.logger.debug('arrayDType: %s', arrayDType)

        #create conents of VRT-file pointing to the binary file
        dataType = {"uint8": "Byte", "int8": "Byte",
                    "uint16": "UInt16", "int16": "Int16",
                    "uint32": "UInt32", "int32": "Int32",
                    "float32": "Float32","float64": "Float64",
                    "complex64": "CFloat64"}.get(str(arrayDType))

        pixelOffset = {"Byte": "1",
                    "UInt16": "2", "Int16": "2",
                    "UInt32": "4", "Int32": "4",
                    "Float32": "4","Float64": "8",
                    "CFloat64": "8"}.get(dataType)

        self.logger.debug('DataType: %s', dataType)

        lineOffset = str(int(pixelOffset)*arrayShape[1])
        contents = self.RawRasterBandSource.substitute(
                                        XSize=arrayShape[1],
                                        YSize=arrayShape[0],
                                        DataType=dataType,
                                        BandNum=1,
                                        SrcFileName=binaryFile,
                                        PixelOffset=pixelOffset,
                                        LineOffset=lineOffset)
        #write XML contents to
        self.write_xml(contents)

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
            vrt = VRT(vrtDataset=self.dataset, geolocation=self.geoloc)
        except:
            # shallow copy (only geometadata)
            vrt = VRT(gdalDataset=self.dataset, geolocation=self.geoloc)

        return vrt

    def add_geolocation(self, geoloc=Geolocation()):
        ''' Add GEOLOCATION to the VRT

        Parameters:
            geoloc: Geolocation object

        Modifes:
            add geoloc to self
            Sets GEOLOCATION metadata
        '''
        self.geoloc = geoloc

        # add GEOLOCATION metadata (empty if geoloc is empty)
        self.dataset.SetMetadata(geoloc.d, 'GEOLOCATION')

    def resized(self, xSize, ySize, resamplingAlg=1):
        '''Resize VRT

        Create Warped VRT with modidied RasterXSize, RasterYSize, GeoTransform
        Parameters:
        -----------
            xSize, ySize: int
                new size of the VRT object
            resamplingAlg:
                0, 1, 2 stands for nearest, bilinear, cubic

        Returns:
        --------
            Resized VRT object
        '''
        # modify GeoTransform: set resolution from new X/Y size
        geoTransform = (0,
                        float(self.dataset.RasterXSize) / float(xSize),
                        0,
                        self.dataset.RasterYSize,
                        0,
                        - float(self.dataset.RasterYSize) / float(ySize))

        # update size and GeoTranform in XML of the warped VRT object
        warpedVRT = self.create_warped_vrt(xSize=xSize, ySize=ySize,
                                           geoTransform=geoTransform,
                                           use_geoloc=False,
                                           use_gcps=False,
                                           use_geotransform=False,
                                           resamplingAlg=resamplingAlg)
        # add source VRT (self) to the warpedVRT
        # in order not to loose RAW file from self
        warpedVRT.srcVRT = self

        return warpedVRT


    def _modify_warped_XML(self, rasterXSize=0, rasterYSize=0, geoTransform=None):
        ''' Modify rasterXsize, rasterYsize and geotranforms in the warped VRT

        Parameters
        ----------
            rasterXSize: integer
                desired X size of warped image
            rasterYSize: integer
                desired Y size of warped image
            geoTransform: tuple of 6 integers
                desired GeoTransform size of the warped image

        Modifies
        --------
            XML of the self VRT file: size and geotranform is updated

        '''
        warpedXML = self.read_xml()

        node0 = Node.create(warpedXML)

        if rasterXSize > 0:
            node0.replaceAttribute("rasterXSize", str(rasterXSize))
        if rasterYSize > 0:
            node0.replaceAttribute("rasterYSize", str(rasterYSize))


        if geoTransform is not None:
            invGeotransform = gdal.InvGeoTransform(geoTransform)
            # convert proper string style and set to the GeoTransform element
            node0.node("GeoTransform").value = str(geoTransform).strip("()")
            node0.node("DstGeoTransform").value = str(geoTransform).strip("()")
            node0.node("DstInvGeoTransform").value = str(
                                                invGeotransform[1]).strip("()")

            if node0.node("SrcGeoLocTransformer"):
                node0.node("BlockXSize").value = str(rasterXSize)
                node0.node("BlockYSize").value = str(rasterYSize)

        self.write_xml(str(node0.rawxml()))

    def _remove_geotransform(self):
        '''Remove GeoTransfomr from VRT Object'''
        # read XML content from VRT
        tmpVRTXML = self.read_xml()
        # find and remove GeoTransform
        node0 = Node.create(tmpVRTXML)
        node1 = node0.delNode("GeoTransform")
        # Write the modified elemements back into temporary VRT
        self.write_xml(str(node0.rawxml()))

    def _add_gcp_metadata(self):
        '''Add GCPs to metadata (required e.g. by Nansat.export())
        Creates string representation of GCPs line/pixel/X/Y
        Add these string to metadata

        Modifies:
            Adds self.vrd.dataset.Metadata
        '''
        gcpNames = ['GCPPixel', 'GCPLine', 'GCPX', 'GCPY']
        gcps = self.dataset.GetGCPs()
        srs = self.dataset.GetGCPProjection()
        chunkLength = 5000

        # if GCPs exist
        if len(gcps) > 0:
            # add GCP Projection
            self.dataset.SetMetadataItem('NANSAT_GCPProjection', srs.replace(",", "|").replace('"', "&"))

            # make empty strings
            gspStrings  = ['', '', '', '']

            # fill string with values
            for gcp in gcps:
                gspStrings[0] = '%s%05d| '    % (gspStrings[0], int(gcp.GCPPixel))
                gspStrings[1] = '%s%05d| '    % (gspStrings[1], int(gcp.GCPLine))
                gspStrings[2] = '%s%012.8f| ' % (gspStrings[2], gcp.GCPX)
                gspStrings[3] = '%s%012.8f| ' % (gspStrings[3], gcp.GCPY)

            for i, gspString in enumerate(gspStrings):
                #split string into chunks
                numberOfChunks = int(float(len(gspString)) / chunkLength)
                chunki = 0
                for chunki in range(0, numberOfChunks+1):
                    chunk = gspString[(chunki * chunkLength):min(((chunki+1) * chunkLength), len(gspString))]
                    # add chunk to metadata
                    self.dataset.SetMetadataItem('NANSAT_%s_%03d' % (gcpNames[i], chunki), chunk)

    def create_warped_vrt(self, dstSRS=None, resamplingAlg=0,
                                xSize=0, ySize=0, geoTransform=None,
                                use_geoloc=True, use_gcps=True, use_geotransform=True,
                                dstGCPs=[], dstGeolocation=None):
        ''' Create VRT object with WarpedVRT

        Modifies the input VRT according to the input options
        Creates simple WarpedVRT with AutoCreateWarpedVRT
        Modifies the WarpedVRT according to the input options

        The function tries to use geolocation by default; if not present (or
        canceled) tries to use GCPs; if not present (or canceled) tries to use
        GeoTransform (either from input dataset or calculates a new one with
        dx=1,dy=-1). Three switches (use_geoloc, use_gcps, use_geotransform)
        allow to select which method to apply for warping. E.g.:
        # #1: srcVRT has Geolocation, geolocation is used
        warpedVRT = srcVRT.create_warped_vrt(dstSRS, xSize, ySize, geoTransform)
        # #2: srcVRT has Geolocation, geolocation is not used, either GCPs (if
        # present) or GeoTransform is used
        warpedVRT = srcVRT.create_warped_vrt(dstSRS, xSize, ySize, geoTransform,
                                             use_geoloc=False)
        # #3: srcVRT has Geolocation or GCPs, geolocation is not used, and GCPs
        # are not used either. Only input GeoTranform is used
        warpedVRT = srcVRT.create_warped_vrt(dstSRS, xSize, ySize, geoTransform,
                                             use_geoloc=False, use_gcps=False)

        # #4: srcVRT has whatever georeference, geolocation is not used, GCPs
        # are not used, GeoTransform is not used either.
        # Artificial GeoTranform is calculated: (0, 1, 0, srcVRT.xSize, -1)
        # Warping becomes pure affine resize
        warpedVRT = srcVRT.create_warped_vrt(dstSRS, xSize, ySize, geoTransform,
                                             use_geoloc=False, use_gcps=False.,
                                             use_geotransform=false)

        If destination image has GCPs (provided in <dstGCPs>): fake GCPs for
        referencing line/piex of SRC image and X/Y of DST image are created and
        added to the SRC image. After warping dstGCPs are added to the WarpedVRT

        If destination image has geolocation (provided in <dstGeolocation>):
        this geolocation is added to the WarpedVRT


        Parameters:
        -----------
        dstSRS: string
            WKT of the destination projection
        resamplingAlg: int
            0, 1, 2 stands for: NearestNeigbour, Bilinear, Cubic
        xSize, ySize: int
            width and height of the destination rasetr
        geoTransform: tuple with 6 floats
            destination GDALGeoTransfrom
        dstGCPs: list with GDAL GCPs
            GCPs of the destination image
        dstGeolocation: Geolocation object
            Geolocation of the destination object
        use_geoloc: Boolean (True)
            Use geolocation in input dataset (if present) for warping
        use_gcps: Boolean (True)
            Use GCPs in input dataset (if present) for warping
        use_geotransform: Boolean (True)
            Use GeoTransform in input dataset for warping or make artificial
            GeoTransform: (0, 1, 0, srcVRT.xSize, -1)

        Returns:
        --------
        warpedVRT: VRT object with WarpedVRT
        '''
        # VRT to be warped
        srcVRT = self.copy()

        # if destination GCPs are given: create and add fake GCPs
        if len(dstGCPs) > 0:
            fakeGCPs = srcVRT._create_fake_gcps(dstGCPs)
            srcVRT.dataset.SetGCPs(fakeGCPs['gcps'], fakeGCPs['srs'])
            # don't use geolocation
            use_geoloc = False
            dstSRS=None

        # prepare VRT.dataset for warping.
        # Select if GEOLOCATION, or GCPs, or GeoTransform from the original
        # dataset are used
        if len(self.geoloc.d)>0 and use_geoloc:
            # use GEOLOCATION by default
            # (remove GCP and GeoTransform)
            srcVRT.dataset.SetGCPs([], '')
            srcVRT._remove_geotransform()
        elif len(srcVRT.dataset.GetGCPs()) > 0 and use_gcps:
            # fallback to GCPs
            # (remove Geolocation and GeoTransform)
            srcVRT.dataset.SetMetadata('', 'GEOLOCATION')
            srcVRT._remove_geotransform()
        elif use_geotransform:
            # fallback to GeoTransform in input VRT
            # (remove Geolocation and GCP)
            srcVRT.dataset.SetMetadata('', 'GEOLOCATION')
            srcVRT.dataset.SetGCPs([], '')
        else:
            # fallback to simplest GeoTransform
            # (remove Geolocation and GCP and replace GeoTransform)
            srcVRT.dataset.SetMetadata('', 'GEOLOCATION')
            srcVRT.dataset.SetGCPs([], '')
            srcVRT.dataset.SetGeoTransform((0, 1, 0, srcVRT.dataset.RasterYSize, 0, -1))

        # create Warped VRT GDAL Dataset
        self.logger.debug('Run AutoCreateWarpedVRT...')
        warpedVRT = gdal.AutoCreateWarpedVRT(srcVRT.dataset, None, dstSRS, resamplingAlg)

        # check if Warped VRT was created
        if warpedVRT is None:
            raise AttributeError("Cannot create warpedVRT!")

        # create VRT object from Warped VRT GDAL Dataset
        self.logger.debug('create VRT object from Warped VRT GDAL Dataset')
        warpedVRT = VRT(vrtDataset=warpedVRT)

        # set x/y size, geoTransform
        self.logger.debug('set x/y size, geoTransform')
        warpedVRT._modify_warped_XML(xSize, ySize, geoTransform)

        # if given, add dst GCPs
        self.logger.debug('if given, add dst GCPs')
        if len(dstGCPs) > 0:
            warpedVRT.dataset.SetGCPs(dstGCPs, dstSRS)
            warpedVRT._remove_geotransform()
            warpedVRT.dataset.SetProjection('')

        # if given, add dst Geolocation
        self.logger.debug('# if given, add dst Geolocation')
        if dstGeolocation is not None:
            warpedVRT._remove_geotransform()
            warpedVRT.add_geolocation(dstGeolocation)
            warpedVRT.dataset.SetProjection('')

        # replace the reference from srcVRT to self
        self.logger.debug('replace the reference from srcVRT to self')
        rawFileName = str(os.path.basename(self.fileName))
        srcFileName = str(os.path.basename(srcVRT.fileName))
        xmlString = str(warpedVRT.read_xml())
        xmlString = str(xmlString.replace(srcFileName, rawFileName))
        warpedVRT.write_xml(xmlString)

        return warpedVRT

    def _create_fake_gcps(self, gcps):
        '''Create GCPs with reference self.pixel/line ==> dst.pixel/line

        GCPs from a destination image (dstGCP) are converted to a gcp of source
        image (srcGCP) this way:

        srcGCPPixel = srcPixel
        srcGCPLine = srcLine
        srcGCPX = dstGCPPixel = f(srcSRS, dstGCPX, dstGCPY)
        srcGCPY = dstGCPLine = f(srcSRS, dstGCPX, dstGCPY)

        Parameters:
        -----------
        gcps: list
            GDAL GCPs
        Returns:
        --------
        gcps: dict
            {'gcps': list with GDAL GCPs, 'srs': fake stereo WKT}
        '''

        # prepare pure lat/lon WKT
        sr = osr.SpatialReference()
        sr.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
        latlongWKT = sr.ExportToWkt()

        # get source SRS (either Projection or GCPProjection)
        srcWKT = self.dataset.GetProjection()
        if srcWKT == '':
            srcWKT = self.dataset.GetGCPProjection()

        # the transformer converts lat/lon to pixel/line of SRC image
        srcTransformer = gdal.Transformer(
                             self.dataset, None,
                             ['SRC_SRS=' + srcWKT,
                             'DST_SRS=' + latlongWKT])

        # create 'fake' GCPs
        for g in gcps:
            # transform DST lat/lon to SRC pixel/line
            succ, point = srcTransformer.TransformPoint(1, g.GCPX, g.GCPY)
            srcPixel = point[0]
            srcLine = point[1]

            # swap coordinates in GCPs:
            # pix1/line1 -> lat/lon  =>=>  pix2/line2 -> pix1/line1
            g.GCPX = g.GCPPixel
            g.GCPY = g.GCPLine
            g.GCPPixel = srcPixel
            g.GCPLine = srcLine

        # create 'fake' STEREO projection for 'fake' GCPs of SRC image
        srsString = ("+proj=stere +lon_0=0 +lat_0=0 +k=1 "
                     "+ellps=WGS84 +datum=WGS84 +no_defs ")
        stereoSRS = osr.SpatialReference()
        stereoSRS.ImportFromProj4(srsString)
        stereoSRSWKT = stereoSRS.ExportToWkt()

        return {'gcps': gcps, 'srs': stereoSRSWKT}
