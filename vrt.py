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

from osgeo import gdal
from osgeo import osr

from nansat_tools import add_logger, Node

class Geolocation():
    '''Container for GEOLOCATION data'''
    def __init__(self, xVRT=None, yVRT=None, xBand=1, yBand=1,
                        srs="+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs",
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
            srs: str, projection WKT
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
            self.yVRT = xVRT
            self.d['Y_DATASET'] = yVRT.fileName
        
        # proj4 to WKT
        sr = osr.SpatialReference()
        sr.ImportFromProj4(srs)
        self.d['SRS'] = sr.ExportToWkt()
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
        self.logger.debug('input vrtDataset: %s' % str(vrtDataset))
        # copy content of the provided VRT dataset
        if vrtDataset is not None:
            self.logger.debug('Making copy of %s ' % str(vrtDataset))
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
            # write file contents
            self.dataset.FlushCache()

        # add geolocation from input geolocation 
        # if not None: overwrite geoloc from vrt or gdal datasets
        if geolocation is not None:
            self.add_geolocation(geolocation)

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
            if it is 0, it means VRT RawRasterband
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
            self.logger.debug('Creating band %s', str(bandDict))
            self._create_band(bandDict["source"], bandDict["sourceBand"],
                    bandDict["wkv"], bandDict.get("parameters", {}), NODATA, LUT, SourceType)
            self.logger.debug('OK!')
        self.dataset.FlushCache()

    def _create_band(self, source, sourceBands, wkv, parameters, NODATA="", LUT="", SourceType='ComplexSource'):
        ''' Function to add a band to the VRT from a source.
        See function _create_bands() for explanation of the input parameters
        '''
        self.logger.info('INPUTS: %s, %s %s %s" ' % (str(source), str(sourceBands), str(wkv), str(parameters)))
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
                       "SourceFilename=%s" % source[0],
                       "ImageOffset=%f" % parameters["ImageOffset"],
                       "PixelOffset=%f" % parameters["PixelOffset"],
                       "LineOffset=%f" % parameters["LineOffset"],
                       "ByteOrder=%s" % parameters["ByteOrder"]]
            self.dataset.AddBand(parameters["dataType"], options)
            return
        # else
        else:
            # Find datatype and blocksizes
            srcDataset = gdal.Open(source[0])
            srcRasterBand = srcDataset.GetRasterBand(sourceBands[0])
            blockXSize, blockYSize = srcRasterBand.GetBlockSize()
            dataType = srcRasterBand.DataType

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

    def resized(self, xSize, ySize):
        '''Resize VRT

        Create Warped VRT with modidied RasterXSize, RasterYSize, GeoTransform
        Parameters:
        -----------
            xSize, ySize: int
                new size of the VRT object

        Returns:
        --------
            Resized VRT object
        '''
        # create a temporary copy of original VRT and remove GCPs in order to
        # apply only affine tranformation in warping
        tmp = self.copy()
        tmp.dataset.SetGCPs([], None)
        tmp.dataset.SetGeoTransform((0, 1, 0, self.dataset.RasterYSize, 0, -1))
        tmp.dataset.SetProjection('')

        # create simplest Warped VRT GDAL Dataset
        warpedVRT = gdal.AutoCreateWarpedVRT(tmp.dataset, None, None, gdal.GRA_Bilinear)

        # create VRT object from warped VRT
        warpedVRT = VRT(vrtDataset=warpedVRT)

        # modify GeoTransform: set resolution from new X/Y size
        geoTransform = (0,
                        float(self.dataset.RasterXSize) / float(xSize),
                        0,
                        self.dataset.RasterYSize,
                        0,
                        - float(self.dataset.RasterYSize) / float(ySize))

        # update size and GeoTranform in XML of the warped VRT object
        warpedVRT._modify_warped_XML(xSize, ySize, geoTransform)

        # append temporary VRT to the warped VRT
        warpedVRT.rawVRT = tmp

        return warpedVRT


    def _modify_warped_XML(self, rasterXSize, rasterYSize, geoTransform):
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

        invGeotransform = gdal.InvGeoTransform(geoTransform)
        node0 = Node.create(warpedXML)

        node0.replaceAttribute("rasterXSize", str(rasterXSize))
        node0.replaceAttribute("rasterYSize", str(rasterYSize))

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
        ''' remove GeoTransfomr from VRT Object'''
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
