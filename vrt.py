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

from xml.etree.ElementTree import *

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

    def __init__(self, metadata, rawVRTName):
        ''' Set attributes common for all mappers

        Parameters
        ----------
        metadata: metadata
        rawVRTName: file name
            '/vsimem/vsiFile.vrt'

        '''
        self.metadata = metadata
        self.rawVRTName = rawVRTName

    def addPixelFunction(self, pixelFunction, bands, dataset, parameters):
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
            value of 'pixelfunction' attribute for each band
        bands: list
            elements in the list represent band number
        dataset: string
            file name
        parameters: dictionary
            "longname", "units" (+ some optional) keys and their values

        Modifies
        --------
        self.vsiDataset: VRT dataset
            add PixelFunctions from bands in the same dataset

        '''
        newBand = self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount)
        options = ['subClass=VRTDerivedRasterBand',
                   'PixelFunctionType='+pixelFunction]
        self.vsiDataset.AddBand(datatype=gdal.GDT_Float32, options=options)
        md = {}
        for i, bandNo in enumerate(bands):
            BlockXSize, BlockYSize = self.dataset.GetRasterBand(bandNo).\
                                                  GetBlockSize()
            DataType = self.dataset.GetRasterBand(bandNo).DataType
            md['source_'+str(i)] = self.SimpleSource.substitute(
                                        XSize=self.vsiDataset.RasterXSize,
                                        BlockXSize=BlockXSize,
                                        BlockYSize=BlockYSize,
                                        DataType=DataType,
                                        YSize=self.vsiDataset.RasterYSize,
                                        Dataset=dataset, SourceBand=bandNo)

        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).\
                        SetMetadata(md, 'vrt_sources')

        for parameter, value in parameters.items():
            self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).\
                            SetMetadataItem(parameter, value)
        self.vsiDataset.GetRasterBand(self.vsiDataset.RasterCount).\
                        SetMetadataItem('pixelfunction', pixelFunction)
        # Took 5 hours of debugging to find this one!!!
        self.vsiDataset.FlushCache()

    SimpleSource = Template('''
            <SimpleSource>
                <SourceFilename relativeToVRT="0">$Dataset</SourceFilename>
                <SourceBand>$SourceBand</SourceBand> \
                <SourceProperties RasterXSize="$XSize" RasterYSize="$YSize" \
                        DataType="$DataType" BlockXSize="$BlockXSize" \
                        BlockYSize="$BlockYSize"/>\
                <SrcRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
                <DstRect xOff="0" yOff="0" xSize="$XSize" ySize="$YSize"/>
            </SimpleSource> ''')

    def createVRT_and_add_bands(self, dataset, metaDict, vrtBandList):
        ''' Create VSI VRT dataset and set geo-metadata from source file

        Parameters
        -----------
        dataset : dataset
            dataset or subdataset if the data has subdataset.
        metaDict: list
            incldes some dictionaries.
            The number of dictionaries is same as number of bands.
            Each dictionary represents metadata for each band.
        vrtBandList: list
            band numbers to fetch.
            If it is None, all bands in the file are fetched.

        Modifies
        --------
        set attributes (vsiDataset, srcRasterXSize, srcRasterYSize,
        dataset, vsiDataset)

        '''
        self.vsiDataset, srcRasterXSize, srcRasterYSize = \
                self._createVRT(dataset, len(vrtBandList))
        self.dataset = dataset
        # To be able to access bands of dataset, to know Blocksize and DataType
        # add bands with metadata and corresponding values
        self.vsiDataset = self._addAllBands(self.vsiDataset, vrtBandList,
                                            metaDict, srcRasterXSize,
                                            srcRasterYSize)

    def _addAllBands(self, vsiDataset, vrtBandList, metaDict,
                     srcRasterXSize, srcRasterYSize):
        '''Loop through all bands and add metadata and band XML source

        Parameters
        -----------
        vsiDataset: VRT dataset
            VSI VRT dataset with common parameters
        vrtBandList: list
            band numbers to fetch
        metaDict: list
            incldes some dictionaries.
            The number of dictionaries is same as number of bands.
            Each dictionary represents metadata for each band.
        srcRasterXSize, srcRasterYSize: int
            raster XSize and rasterYSize

        Returns
        --------
        vsiDataset: VRT dataset
            VSI VRT dataset added band metadata

        '''
        # for bandNo in vrtBandList:
        for iBand in range(len(vrtBandList)):
            bandNo = vrtBandList[iBand]
            # check if
            if int(bandNo) > int(metaDict.__len__()):
                print ("vrt.addAllBands(): "
                       "an element in the bandList is improper")
                return
            # add metadata
            # !!! This (GetRasterBand(1)) is a just temporary solution
            # because self.dataset means the first subdataset
            # if the data has subdatasets. should be fixed !!!
            ##xBlockSize, yBlockSize = self.dataset.\
            ## 	GetRasterBand(bandNo).GetBlockSize()
            xBlockSize, yBlockSize = self.dataset.\
                                          GetRasterBand(1).GetBlockSize()
            ##srcDataType = self.dataset.GetRasterBand(bandNo).DataType
            srcDataType = self.dataset.GetRasterBand(1).DataType
            wkv_name = metaDict[bandNo-1]["wkv"]
            wkvDict = self._get_wkv(wkv_name)
            for key in wkvDict:
                vsiDataset.GetRasterBand(iBand+1).\
                           SetMetadataItem(key, wkvDict[key])
            if "parameters" in metaDict[bandNo-1]:
                # Warning: zero-indexing for medaDict,
                # but 1-indexing for band numbers!
                for metaName in metaDict[bandNo-1]["parameters"]:
                    vsiDataset.GetRasterBand(iBand+1).\
                               SetMetadataItem(metaName,
                               metaDict[bandNo-1]["parameters"][metaName])
            # add band
            bandSource = self.SimpleSource.\
                              substitute(XSize=srcRasterXSize,
                              YSize=srcRasterYSize,
                              Dataset=metaDict[bandNo-1]['source'],
                              SourceBand=metaDict[bandNo-1]['sourceBand'],
                              BlockXSize=xBlockSize, BlockYSize=yBlockSize,
                              DataType=srcDataType)
            vsiDataset.GetRasterBand(iBand+1).\
                       SetMetadataItem("source_0", bandSource,
                       "new_vrt_sources")
        vsiDataset.FlushCache()

        return vsiDataset

    def _createVRT(self, srcDataset, length):
        ''' Create VRT with common parameters

        Parameters
        ----------
        srcDataset: dataset
            GDAL dataset / subdataset
        length: int
            number of elements in vrtBandList

        Returns
        -------
        vsiDataset: dataset
            VRT dataset which includes common information
        srcRasterXSize, srcRasterYSize: int
            raster Xsize and raster Ysize

        '''
        # get geo-metadata from source dataset (or subdataset)
        srcGeoTransform = srcDataset.GetGeoTransform()
        srcProjection = srcDataset.GetProjection()
        srcProjectionRef = srcDataset.GetProjectionRef()
        srcGCPCount = srcDataset.GetGCPCount()
        srcGCPs = srcDataset.GetGCPs()
        srcGCPProjection = srcDataset.GetGCPProjection()

        srcRasterXSize = srcDataset.GetRasterBand(1).XSize
        srcRasterYSize = srcDataset.GetRasterBand(1).YSize

        # create VSI VRT dataset
        vrtDrv = gdal.GetDriverByName("VRT")
        vsiDataset = vrtDrv.Create(self.rawVRTName,
                                   srcRasterXSize, srcRasterYSize,
                                   length, gdal.GDT_Float32)

        # set geo-metadata in the VSI VRT dataset
        vsiDataset.SetGCPs(srcGCPs, srcGCPProjection)
        vsiDataset.SetProjection(srcProjection)
        vsiDataset.SetGeoTransform(srcGeoTransform)

        # set metadata
        vsiDataset.SetMetadata(self.metadata)

        return vsiDataset, srcRasterXSize, srcRasterYSize

    def _get_wkv(self, wkv_name):
        ''' Get wkv from wkv.xml

        Parameters
        ----------
        wkv_name: string
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
            if e1.find("name").text == wkv_name:
                wkvDict = {"name": wkv_name}
                for e2 in list(e1):
                    wkvDict[e2.tag] = e2.text

        return wkvDict
