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
from string import Template

from xml.etree.ElementTree import ElementTree

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal
    import osr


class VRT():
    '''Top class of Nansat readers

    Set attributes to VRT.
    self.vsiDataset include all band information specified in each mapper.

    '''

    def __init__(self, dataset, metadata, rawVRTName):
        ''' Set attributes common for all mappers

        Parameters
        ----------
        datset: GDAL Dataset
            dataset from Nansat
        metadata: GDAL metadata
            metadata from Nansat
        rawVRTName: string
            file name from Nansat ('/vsimem/vsiFile.vrt')
        '''

        self.metadata = metadata
        self.rawVRTName = rawVRTName
        self.dataset = dataset

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
        self.vsiDataset: VRT dataset
            add PixelFunctions from bands in the same dataset
        '''

        newBand = self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount)
        options = ['subClass=VRTDerivedRasterBand',
                   'PixelFunctionType=' + pixelFunction]
        self.vsiDataset.AddBand(datatype=gdal.GDT_Float32, options=options)
        md = {}
        srcDataset = gdal.Open(fileName)
        for i, bandNo in enumerate(bands):
            srcRasterBand = srcDataset.GetRasterBand(bandNo)
            blockXSize, blockYSize = srcRasterBand.GetBlockSize()
            dataType = srcRasterBand.DataType
            md['source_' + str(i)] = self.SimpleSource.substitute(
                                        XSize=self.vsiDataset.RasterXSize,
                                        BlockXSize=blockXSize,
                                        BlockYSize=blockYSize,
                                        DataType=dataType,
                                        YSize=self.vsiDataset.RasterYSize,
                                        Dataset=fileName, SourceBand=bandNo)

        # set metadata for each destination raster band
        dstRasterBand = self.vsiDataset.GetRasterBand(self.vsiDataset.\
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
        self.vsiDataset.FlushCache()

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

    ComplexSource = Template('''
            <ComplexSource>
                <SourceFilename relativeToVRT="0">$Dataset</SourceFilename>
                <SourceBand>$SourceBand</SourceBand>
                <ScaleOffset>$ScaleOffset</ScaleOffset>
                <ScaleRatio>$ScaleRatio</ScaleRatio>
                <SourceProperties RasterXSize="$XSize" RasterYSize="$YSize"
                        DataType="$DataType" BlockXSize="$BlockXSize"
                        BlockYSize="$BlockYSize"/>
                <SrcRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
                <DstRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
            </ComplexSource> ''')

    def _add_all_bands(self, vrtBandList, metaDict,
                     srcRasterXSize, srcRasterYSize):
        '''Loop through all bands and add metadata and band XML source

        Parameters
        -----------
        vrtBandList: list
            band numbers to fetch
        metaDict: list
            incldes some dictionaries.
            The number of dictionaries is same as number of bands.
            Each dictionary represents metadata for each band.
        srcRasterXSize, srcRasterYSize: int
            raster XSize and rasterYSize

        Returns 0
        --------

        Modifies
        --------
        self.vsiDataset: VRT dataset
            VSI VRT dataset added band metadata
        '''
        for iBand, bandNo in enumerate(vrtBandList):
            # check if the band in the list exist
            if int(bandNo) > int(metaDict.__len__()):
                print ("vrt.addAllBands(): "
                       "an element in the bandList is improper")
                break

            srcRasterBand = gdal.Open(metaDict[bandNo - 1]['source']).\
                       GetRasterBand(metaDict[bandNo - 1]['sourceBand'])

            xBlockSize, yBlockSize = srcRasterBand.GetBlockSize()
            srcDataType = srcRasterBand.DataType

            # set metadata for each destination raster band
            dstRasterBand = self.vsiDataset.GetRasterBand(iBand + 1)
            # set metadata from WKV
            wkvName = metaDict[bandNo - 1]["wkv"]
            dstRasterBand = self._put_metadata(dstRasterBand,
                                                self._get_wkv(wkvName))
            # set metadata from parameters (if exist)
            if "parameters" in metaDict[bandNo - 1]:
                dstRasterBand = self._put_metadata(dstRasterBand,
                                     metaDict[bandNo - 1]["parameters"])

            # set statistics
            vmin, vmax, vmean, vstd = self.dataset.GetRasterBand(\
                                        metaDict[bandNo - 1]['sourceBand']).\
                                        GetStatistics(True, True)
            self.vsiDataset.GetRasterBand(iBand+1).\
                            SetStatistics(vmin, vmax, vmean, vstd)

            # get scale/offset from metaDict (or set default 1/0)
            if 'scale' in metaDict[bandNo - 1]:
                scaleRatio = metaDict[bandNo - 1]['scale']
            else:
                scaleRatio = 1
            if 'offset' in metaDict[bandNo - 1]:
                scaleOffset = metaDict[bandNo - 1]['offset']
            else:
                scaleOffset = 0

            # create band source metadata
            bandSource = self.ComplexSource.\
                              substitute(XSize=srcRasterXSize,
                              YSize=srcRasterYSize,
                              Dataset=metaDict[bandNo-1]['source'],
                              SourceBand=metaDict[bandNo-1]['sourceBand'],
                              BlockXSize=xBlockSize, BlockYSize=yBlockSize,
                              DataType=srcDataType,
                              ScaleOffset=scaleOffset, ScaleRatio=scaleRatio)

            # set band source metadata
            dstRasterBand.SetMetadataItem("source_0", bandSource,
                                          "new_vrt_sources")

        self.vsiDataset.FlushCache()

        return 0

    def _createVRT(self, metaDict, vrtBandList):
        ''' Create VSI VRT dataset, set geo-metadata, add all bands

        Parameters
        -----------
        metaDict: list
            includes some dictionaries for mapping bands in input Datset
            and well-known-variables. The number of dictionaries is same
            as number of bands. Each dictionary represents metadata for
            each band.
        vrtBandList: list
            band numbers to fetch.
            If it is None, all bands in the file are fetched.

        Modifies
        --------
            Creates self.vsiDataset
            Sets in vsiDataset: srcRasterXSize, srcRasterYSize,
            GeoTransform, GCPS, ...
            Add all bands to the vsiDataset
        '''
        # get geo-metadata from source dataset (or subdataset)
        srcGeoTransform = self.dataset.GetGeoTransform()
        srcProjection = self.dataset.GetProjection()
        srcProjectionRef = self.dataset.GetProjectionRef()
        srcGCPCount = self.dataset.GetGCPCount()
        srcGCPs = self.dataset.GetGCPs()
        srcGCPProjection = self.dataset.GetGCPProjection()

        srcRasterXSize = self.dataset.GetRasterBand(1).XSize
        srcRasterYSize = self.dataset.GetRasterBand(1).YSize

        # create VSI VRT dataset
        vrtDrv = gdal.GetDriverByName("VRT")
        self.vsiDataset = vrtDrv.Create(self.rawVRTName,
                                   srcRasterXSize, srcRasterYSize,
                                   len(vrtBandList), gdal.GDT_Float32)

        # set geo-metadata in the VSI VRT dataset
        self.vsiDataset.SetGCPs(srcGCPs, srcGCPProjection)
        self.vsiDataset.SetProjection(srcProjection)
        self.vsiDataset.SetGeoTransform(srcGeoTransform)

        # set metadata
        self.vsiDataset.SetMetadata(self.metadata)

        # add bands with metadata and corresponding values
        self._add_all_bands(vrtBandList, metaDict,
                            srcRasterXSize, srcRasterYSize)

        return 0

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
        for key in metadataDict:
            rasterBand.SetMetadataItem(key, metadataDict[key])

        return rasterBand
