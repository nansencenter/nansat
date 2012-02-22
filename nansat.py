# Name:    nansat.py
# Purpose: main file of the NANSAT module.
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

import cStringIO
from os import path, listdir
from string import maketrans
import sys
import time

import fnmatch
import Image
import ImageDraw
import ImageFont
import ImageOps
import matplotlib.cm as cm
import numpy as np
from scipy.misc import toimage, pilutil
from scipy.misc import imsave
from scipy.stats import cumfreq
from xml.etree.ElementTree import XML, ElementTree, tostring
from pylab import *

import matplotlib.pyplot as plt

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal
    import osr

from domain import Domain
from vrt import *


class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass


class GDALError(Error):
    '''Error from GDAL '''
    pass


class ProjectionError(Error):
    '''Cannot get the projection'''
    pass


class DataError(Error):
    '''Error for data.
        e.g. : empty pixel value array in get_pixelValueRange()'''
    pass


class OptionError(Error):
    '''Error for unproper options (arguments) '''
    pass


class Nansat():
    '''Main of Nansat

    Construct Nansat object that consist of
        basic dataset information (fileName, dataset, metadata etc..),
        VRT file which points to orignal file with satellite data and
        is saved in an XML format in memory (GDAL VSI).
    '''
    def __init__(self, fileName, mapperName='', bandList=None):
        '''Construct Nansat object

        Open GDAL dataset,
        Read metadata,
        Generate GDAL VRT file with mapping of variables in memory

        Parameters
        ----------
        fileName : string
            location of the file
        mapperName : string, optional
            "ASAR", "hurlam", "merisL1", "merisL2", "ncep", "radarsat2",
            "seawifsL2" are currently available.  (27.01.2012)
        bandList : list, optional
            band numbers to fetch.
            If it is None, all bands in the file are fetched.

        Modifies
        --------
        self.fileName : file name
            set file name given by the argument
        self.dataset : GDAL dataset
            set GDAL dataset
        self.metadata : metadata
            set metadata of the dataset
        self.rawVRTFileName : file name
            set '/vsimem/vsiFile.vrt'
        self.warpedVRTFileName : file name
            set '/vsimem/vsi_warped.vrt'
        self.vrtDriver : VRT driver
            set GDAL VRT driver
        self.rawVRT : VRT dataset
            set VRT dataset with mapping of variables
        self.warpedVRT : VRT dataset
            None
        self.vrt : VRT dataset
            Copy of self.rawVRT

        Raises
        ------
            GDALError: occurs when the dataset is None or "".
            GDALError: occurs when the metadata is None or "".

        '''
        # location of the data
        self.fileName = fileName

        # dataset
        self.dataset = gdal.Open(self.fileName)
        if (self.dataset is None) or (self.dataset == ""):
            raise GDALError("Nansat._init_(): Cannot get the dataset from "
                            + self.fileName)

        # metadata
        self.metadata = self.dataset.GetMetadata()
        if (self.metadata is None) or (self.metadata == ""):
            raise GDALError("Nansat._init_(): Cannot get the metdadata")

        # names of raw and warped VRT files in memory
        self.rawVRTFileName = '/vsimem/vsiFile.vrt'
        self.warpedVRTFileName = '/vsimem/vsi_warped.vrt'
        self.vrtDriver = gdal.GetDriverByName("VRT")

        # VRT with mapping of variables
        self.rawVRT = self._get_mapper(mapperName, bandList)
        # Warped VRT
        self.warpedVRT = None
        # Current VRT
        self.vrt = self.rawVRT

    def __getitem__(self, bandNo):
        ''' Returns the band as a NumPy array, by overloading []

        Returns
        -------
            self.get_GDALRasterBand(bandNo).ReadAsArray(): NumPy array

        '''
        return self.get_GDALRasterBand(bandNo).ReadAsArray()

    def __repr__(self):
        '''Prints basic info about the Nansat object to the terminal

        '''
        print '-' * 40
        print self.fileName
        print '-' * 40
        self.list_bands()
        print self.get_domain()
        return ''

    def dereproject(self):
        '''Cancel reprojection

        Modifies
        --------
            self.vrt : VRT dataset
                replaced the raw/underlaying dataset
        '''
        self.vrt = self.rawVRT

    def downscale(self, factor=1, method="average"):
        '''Downscale the size of the data.

        The size of data is downscaled as (xSize/factor, ySize/factor).
        self.vrt is rewritten to the the downscaled sizes.
        If GCPs are given, they are also rewritten.

        Parameters
        ----------
            factor: int, optional
            method: "average" (default) or "subsample" (= nearest neighbor),
                    optional

        Modifies
        --------
            self.vrt : VRT dataset
                raster size are modified to downscaled size.
                If GCPs are given in the dataset, they are also overwritten.

        Raises
        ------
            OptionError: occurs when method is not "average" or "subsample"

        '''
        if not (method == "average" or method == "subsample"):
            raise OptionError("method should be 'average' or 'subsample'")

        # Write the vrt to a VSI-file
        vrtDatasetCopy = self.vrtDriver.CreateCopy(self.rawVRTFileName,
                                                   self.vrt)

        # Get XML content from VSI-file
        vsiFileContent = self._read_write_vsi_file(self.rawVRTFileName)

        # Get element from the XML content and modify some elements
        # using Domain object parameters
        element = XML(vsiFileContent)
        rasterXSize = int(float(element.get("rasterXSize")) / factor)
        rasterYSize = int(float(element.get("rasterYSize")) / factor)
        element.set("rasterXSize", str(rasterXSize))
        element.set("rasterYSize", str(rasterYSize))

        for elem in element.iter("DstRect"):
            elem.set("xSize", str(rasterXSize))
            elem.set("ySize", str(rasterYSize))

        # if method = "average", overwrite "SimpleSource" to "AveragedSource"
        if method == "average":
            for elem in element.iter("SimpleSource"):
                elem.tag = "AveragedSource"

        # Edit GCPs to correspond to the downscaled size
        for elem in element.iter("GCP"):
            pxl = float(elem.get("Pixel")) / factor
            if pxl > float(rasterXSize):
                pxl = rasterXSize
            lin = float(elem.get("Line")) / factor
            if lin > float(rasterYSize):
                lin = rasterYSize
            elem.set("Pixel", str(pxl))
            elem.set("Line", str(lin))

        # Overwrite element
        # Write the modified elemements into VSI-file
        self._read_write_vsi_file(self.rawVRTFileName, tostring(element))

        self.vrt = gdal.Open(self.rawVRTFileName)

    def export_VRT(self, fileName=None):
        '''Export in-memory VRT dataset to a physical file

        If fileName is None, this method is skipped.
        Otherwise, open VSI-file and copy it to a physical file
        whose location is given by the argument.

        Parameters
        ----------
            fileName: string, optional
                location for an output VRT file

        '''
        if fileName is None:
            # GDAL special name to flush output
            fileName = "/vsistdout/"
            # to console. Unfortunately an error message about
            # non-existing file is reported this must be a bug in GDAL.
        vrtDatasetCopy = self.vrtDriver.CreateCopy(fileName, self.vrt)
        vrtDatasetCopy = None

    def get_GDALRasterBand(self, bandNo=1, bandID=None):
        '''Get a GDALRasterBand of a given Nansat object.

        Get a GDALRasterBand specified by the argument.

        If a bandID is given, secify a bandNo based on it.
        Otherwise check if the given bandNo is proper.
        Get a GDALRasterBand from vrt.

        Parameters
        ----------
            bandNo: serial number or string, optional (default is 1)
                if number - a band number of the band to fetch
                if string bandID = {'band_name': bandNo}
            bandID: a dictionary with metadata unique for one band

        Returns
        -------
            self.vrt.GetRasterBand: a GDAL RasterBand

        Raises
        ------
            OptionError: occurs when the bandNo is not a proper number.
        Example
        -------
            b = get_GDALRasterBand(1)
            b = get_GDALRasterBand('sigma0')
            b = get_GDALRasterBand({"short_name":"radiance",
                                    "wavelength":"1240"})


        See Also
        --------
            _specify_bandNo: specify a band number based on the bandID list

        '''

        print "raster num" ,self.rawVRT.RasterCount
        # If bandID is given, bandNo is specified here.
        if bandID is not None:
            bandNo = self._specify_bandNo(bandID)
        # if bandNo is given and it is string fetch band which has
        # band_name == bandNo
        elif isinstance(bandNo, str):
            bandNo = self._specify_bandNo({'band_name': bandNo})
        # if given bandNo is over the existing bands, give error message
        elif (bandNo < 1 or bandNo > self.rawVRT.RasterCount):
            raise OptionError("Nansat.get_GDALRasterBand(): "
                             "bandNo takes from 1 to",
                             self.rawVRT.RasterCount)

        # Based on bandNo,
        # the GDAL RasterBand of the corresponding band is returned
        return self.vrt.GetRasterBand(bandNo)

    def list_bands(self):
        '''Show band information of the given Nansat object

        Show serial number, longName, name and all parameters
        for each band in the metadata of the given Nansat object.

        '''
        for iBand in range(self.rawVRT.RasterCount):
            metadata = self.rawVRT.GetRasterBand(iBand + 1).GetMetadata()
            print "Band :", iBand + 1
            for i in metadata:
                if i != "units":
                    print "  %s: %s" % (i,
                          self.rawVRT.GetRasterBand(iBand + 1).\
                          GetMetadataItem(i))

    def reproject(self, dstDomain=None, resamplingAlg=0):
        '''Reproject the object based on the given Domain

        Warp the raw VRT using AutoCreateWarpedVRT() using projection
        from the Domain.
        Modify XML content of the warped vrt using the Domain parameters.
        Generate self.warpedVRT and replace self.vrt to warpedVRT.

        Parameters
        ----------
            dstDomain: domain
                destination Domain where projection and resolution are set

        Modifies
        --------
            self.warpedVRT: VRT dataset
                warped vrt dataset
            self.vrt: VRT dataset
                replaced to warped VRT dataset

        Raises
        ------
            ProjectionError: occurs when the projection of the source data
            is None.
            ProjectionError: occurs when the projection of the target data
            is None.
            OptionError: occures when the option combination is not proper.
            AttributeError: occurs when it is impossible to get warpedVRT.

        See Also
        --------
            http://www.gdal.org/gdalwarp.html

        '''
        # Get source SRS (either Projection or GCPProjection)
        srcWKT = self.rawVRT.GetProjection()
        if srcWKT == '':
            srcWKT = self.rawVRT.GetGCPProjection()

        if srcWKT == '':
            raise ProjectionError("Nansat.reproject(): "
                                  "rawVrt.GetProjection() is None")

        if dstDomain is not None:
            print 'Reprojection with given Domain'
            # warp the Raw VRT onto the coordinate stystem given by
            # input Domain
            rawWarpedVRT = gdal.AutoCreateWarpedVRT(
                                   self.rawVRT, srcWKT,
                                   dstDomain.memDataset.GetProjection(),
                                   resamplingAlg)

        else:
            # Erroneous input options
            raise OptionError("Nansat.reproject(): "
                              "wrong combination of input options")

        # modify extent of the created Warped VRT
        self.warpedVRT = self._modify_warpedVRT(
                                   rawWarpedVRT,
                                   dstDomain.memDataset.RasterXSize,
                                   dstDomain.memDataset.RasterYSize,
                                   dstDomain.memDataset.GetGeoTransform())

        # check created Warped VRT
        if self.warpedVRT is None:
            raise AttributeError("Nansat.reproject():cannot get warpedVRT")

        # set default vrt to be the warped one
        self.vrt = self.warpedVRT

    def reproject_on_gcp(self, gcpImage, resamplingAlg=0):
        ''' Reproject the object onto the input object with gcps
        NB! This is a test function required urgently for the open-wind
        project. It is tesed only on NCEP or GLOBAL DEM and
        RADARSAT2 or MERIS images and should be refined and
        added to Nansat.reproject()

        Parameters
        ----------
            gcpImage: Nansat object of an image with GCPs
            resamplingAlg: integer, option for AutoCreateWarpedVRT

        Modifies
        --------
            self.warpedVRT: VRT dataset
                new warped vrt
            self.vrt: VRT dataset
                replaced warped VRT

        '''
        #name of VRT with 'fake' GCPs
        tmpVRTName = '/vsimem/vsiFileFakeGCP.vrt'

        # prepare pure lat/lon WKT
        proj4string = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs"
        latlongSRS = osr.SpatialReference()
        latlongSRS.ImportFromProj4(proj4string)
        latlongWKT = latlongSRS.ExportToWkt()

        # get source SRS (either Projection or GCPProjection)
        srcWKT = self.vrt.GetProjection()
        if srcWKT == '':
            srcWKT = self.vrt.GetGCPProjection()

        # the transformer converts lat/lon to pixel/line of SRC image
        srcTransformer = gdal.Transformer(
                             self.vrt, None,
                             ['SRC_SRS=' + srcWKT,
                             'DST_SRS=' + latlongWKT])

        # get GCPs from DST image
        gcps = gcpImage.vrt.GetGCPs()

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

        # make copy of the RAW VRT file and replace GCPs
        tmpVRT = self.vrtDriver.CreateCopy(tmpVRTName, self.rawVRT)

        # create 'fake' STEREO projection for 'fake' GCPs of SRC image
        srsString = ("+proj=stere +lon_0=0 +lat_0=0 +k=1 "
                     "+ellps=WGS84 +datum=WGS84 +no_defs ")
        stereoSRS = osr.SpatialReference()
        stereoSRS.ImportFromProj4(srsString)
        stereoSRSWKT = stereoSRS.ExportToWkt()
        tmpVRT.SetGCPs(gcps, stereoSRSWKT)
        tmpVRT.SetProjection('')
        tmpVRT = None

        # remove GeoTransfomr from SRC image
        # read XML content from VSI-file
        vsiFileContent = self._read_write_vsi_file(tmpVRTName)

        # find and remove GeoTransform
        tree = XML(vsiFileContent)
        elemGT = tree.find("GeoTransform")
        tree.remove(elemGT)

        # Write the modified elemements back into VSI-file
        self._read_write_vsi_file(tmpVRTName, tostring(tree))

        # create warped vrt out of tmp vrt
        tmpVRT = gdal.Open(tmpVRTName)
        rawWarpedVRT = gdal.AutoCreateWarpedVRT(tmpVRT, stereoSRSWKT,
                                                stereoSRSWKT,
                                                resamplingAlg)

        # change size and geotransform to fit the DST image
        self.warpedVRT = self._modify_warpedVRT(
                                   rawWarpedVRT,
                                   gcpImage.vrt.RasterXSize,
                                   gcpImage.vrt.RasterYSize,
                                   (0, 1, 0, 0, 0, 1))
        self.vrt = self.warpedVRT

    def write_figure(self, fileName, bands=[1], bandID=None, colormap = cm.jet,
                     cmin=[None], cmax=[None], ratio=[1.0],
                     logarithm=False, gamma=2,
                     numOfTicks=5, fontSize=10,
                     margin=30, pad=60, textwidthMax=600,
                     colorbar_height=10,
                     legend=None,
                     useFullMatrix=False, extension='png'):

        '''Save a raster band to a figure in grapfical format.

        Get numpy array from the band(s) and band information specified
        either by given band number or band id.
        -- If three bands are given, merge them and create PIL image.
        -- If one band is given, adjust the array brightness and contrast
         using the given min/max or histogram ratio.
         If logarithm is True, the array is converted based on tone curve.
         Create PIL images for the image and colorbar.
        Create and output PIL image and paste the image and legend
        (and the colorbar).
        Save the PIL output image in PNG.

        Parameters
        ----------
            fileName: string
                Output file name
            bands : list, optional, default = [1]
                the size of the list has to be 1 or 3.
                if the size is 3, RGB image is created based on the three bands.
                Then the first element is Red, the second is Green,
                and the third is Blue.
            colormap : colormap, optional, default = cm.jet

            cmin, cmax : list, optional
                minimum and maximum pixel values of each band
            ratio : listimageDatatype
                ratio of number of the pixels
                used to compute cmin and cmax.
            logarithm : boolean, optional
                If True, tone curve is used to convert pixel values.
                If False, linear.
            gamma : float ( > 0.0), optional, default = 2.0
                coefficient of the tone curve
            numOfTicks : int, optional, default = 5
                number of ticks
            fontSize : int
                size of font for the legend
            margin: int
                margin of the output file
            pad : int
                height for the legend
            textwidthMax : int
                width for the legend
            colorbar_height : int
                height for the colorbar
            legend : text
            useFullMatrix : boolean
            extension : extension
                extension of the outputfile


            bandName: a list, optional
                (e.g.: bandIdList = {"name":"radiance", "wavelength":"645"})
            thresholdRatio: float (0.0 - 1.0), optional
                e.g. : thresholdRatio = 0.95 means to round off 5%
                        form the both sides (upper and lower sides).
            useFullMatrix: boolean, optional
                if true, the full matrix is used for estimating min/max,
                otherwise only image scaled down to 100x100 (MUCH FASTER)

        Raises
        ------
            OptionError: occurs when number of elements of bands list is not
                          1 or 3
            OptionError: occurs when number of elements of ratio list is not
                          1 or 3
            OptionError: occurs when number of elements of bands list is
                          different from that of ratio list
                          and the element is not 1.0.
            OptionError: occurs when number of elements of cmin list is not
                          1 or 3
            OptionError: occurs when number of elements of cmin list is
                          different from that of cmin list
                          and the element is not None.
            OptionError: occurs when number of elements of cmax list is not
                          1 or 3
            OptionError: occurs when number of elements of cmax list is
                          different from that of cmax list
                          and the element is not None.
            OptionError: occurs when gammma is 0 or negative
            DataError: occurs when the array of the band is empty

        See also
        ------
        http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        '''

        if bandID is not None:
            bands = [self._specify_bandNo(bandID)]

        if not ((len(bands) == 3) or (len(bands) == 1)):
            raise OptionError("number of elements of bands list must be 1 or 3")

        if not ((len(ratio) == 3) or (len(ratio) == 1)):
            raise OptionError("number of elements of ratio list must be 1 or 3")
        elif len(bands) > len(ratio) and ratio[0] == 1.0:
            ratio = [1.0, 1.0, 1.0]
        elif (len(bands) != len(ratio)):
            raise OptionError("number of elements of ratio list must"\
                               " be same with that of bands list")

        if not ((len(cmin) == 3) or (len(cmin) == 1)):
            raise OptionError("number of elements of cmin list must be 1 or 3")
        elif len(bands) > len(cmin) and cmin[0] == None:
            cmin = [None, None, None]
        elif (len(bands) != len(cmin)):
            raise OptionError("number of elements of cmin list must"\
                               " be same with that of bands list")

        if not ((len(cmax) == 3) or (len(cmax) == 1)):
            raise OptionError("number of elements of cmax list must be 1 or 3")
        elif len(bands) > len(cmax) and cmax[0] == None:
            cmax = [None, None, None]
        elif (len(bands) != len(cmax)):
            raise OptionError("number of elements of cmax list must"\
                               " be same with that of bands list")

        if gamma <= 0.0:
            raise OptionError("gamma argument must be a positive float")

        rasterXSize = self.vrt.RasterXSize
        rasterYSize = self.vrt.RasterYSize

        print "Reading figure(s)  size = (", rasterXSize, "x", rasterYSize, ")"
        array = np.array([])

        tic = time.clock()
        for iBand in bands:
            iArray = self[iBand]
            if iArray is None:
                raise DataError("Nansat.write_figure(): "
                                "array of the band is empty")
            array = np.append(array, iArray)
        array = array.reshape(len(bands), rasterYSize, rasterXSize)

        pixvalMinList = []
        pixvalMaxList = []
        pixvalHistBuckets = []
        pixvalHistScale = []
        pixvalHist = np.array([])

        for i, iBand in enumerate(bands):
            bandPixvalueInfo = self.vrt.GetRasterBand(iBand).GetDefaultHistogram()
            pixvalMinList.append(float(bandPixvalueInfo[0]))
            pixvalMaxList.append(float(bandPixvalueInfo[1]))
            pixvalHistBuckets.append(float(bandPixvalueInfo[2]))
            pixvalHistScale.append((pixvalMaxList[i] - pixvalMinList[i]) / pixvalHistBuckets[i])
            pixvalHist = np.append(pixvalHist, (np.array(bandPixvalueInfo[3], dtype=float)))
        pixvalHist = pixvalHist.reshape(len(bands), pixvalHistBuckets[0])

        pixvalMin = pixvalMinList[0]
        pixvalMax = pixvalMaxList[0]

        #------------------------------------------------------------------------------#
        # <OPTION> Set minimum and maximum pixel values
        # -- Get cmin and cmax and clip the pixel values
        numOfPixels = rasterXSize*rasterYSize
        for iBand in range(len(bands)):
            if ratio[iBand] != 1.0:
                accumulateHist = add.accumulate(pixvalHist[iBand, :]) / numOfPixels
                if (accumulateHist[0] > (1-ratio[iBand])) and (accumulateHist[0] == accumulateHist[1]):
                    accumulateHist = accumulateHist[accumulateHist != accumulateHist[0]]
                    cmin[iBand] = pixvalMaxList[iBand] - pixvalHistScale[iBand] * len(accumulateHist)
                    cmax[iBand] = pixvalMaxList[iBand]
                elif accumulateHist[0] == accumulateHist[1]:
                    accumulateHistMin = accumulateHist[accumulateHist != accumulateHist[0]]
                    histMin = accumulateHist[accumulateHist < (accumulateHist[0] + (1 - (ratio[iBand] + accumulateHist[0])) / 2)]
                    histMax = accumulateHist[accumulateHist > (1 - (1 - (ratio[iBand] + accumulateHist[0])) / 2)]
                    cmin[iBand] = pixvalMinList[iBand] + pixvalHistScale[iBand] * len(histMin)
                    cmax[iBand] = pixvalMinList[iBand] + pixvalHistScale[iBand] * (len(accumulateHist) - len(histMax))
                else:
                    histMin = accumulateHist[accumulateHist < ((1 - ratio[iBand]) / 2)]
                    histMax = accumulateHist[accumulateHist > (1 - ((1 - ratio[iBand]) / 2))]
                    cmin[iBand] = pixvalMinList[iBand] + pixvalHistScale[iBand] * len(histMin)
                    cmax[iBand] = pixvalMinList[iBand] + pixvalHistScale[iBand] * (len(accumulateHist) - len(histMax))

            if (cmax[iBand] or cmin[iBand]) is not None:
                if (cmax[iBand] is not None) and (cmax[iBand] < pixvalMaxList[iBand]):
                    pixvalMax = cmax[iBand]
                if (cmin[iBand] is not None) and (cmin[iBand] > pixvalMinList[iBand]):
                    pixvalMin = cmin[iBand]
                array[iBand, :, :] = clip(array[iBand, :, :], pixvalMin, pixvalMax)

# If three bands are given, marge them and crate a PIL image
        if len(bands)==3:
# -- Normalize RGB arrays to the interval [0,255]
            if array[0, :, :].max() == array[0, :, :].min():
                raise OptionError("min. and max. pixel valuses in Red band are same. check ratio, cmin and cmax!!")
            else:
                arrayR = (array[0, :, :] - array[0, :, :].min()) *255/ (array[0, :, :].max()- array[0, :, :].min())

            if array[1, :, :].max() == array[1, :, :].min():
                raise OptionError("min. and max. pixel valuses in Green band are same. check ratio, cmin and cmax!!")
            else:
                arrayG = (array[1, :, :] - array[1, :, :].min()) *255/ (array[1, :, :].max()- array[1, :, :].min())

            if array[2, :, :].max() == array[2, :, :].min():
                raise OptionError("min. and max. pixel valuses in Blue band are same. check ratio, cmin and cmax!!")
            else:
                arrayB = (array[2, :, :] - array[2, :, :].min()) *255/ (array[2, :, :].max()- array[2, :, :].min())

            pilImgR = Image.fromarray(np.uint8(arrayR))
            pilImgG = Image.fromarray(np.uint8(arrayG))
            pilImgB = Image.fromarray(np.uint8(arrayB))
# -- Merge three bands
            pilImgFig = Image.merge("RGB", (pilImgR, pilImgG, pilImgB))
# If one band is given, create its PIL image and color bar PIL image
        else:
# <OPTION> Convart array based on logarithmic tone curve
            if logarithm:
                pixvalPositiveMin = array[0, :, :][array[0, :, :] > 0].min()
                try:
                    pixvalPositiveMin = array[0, :, :][array[0, :, :] > 0].min()
                    if (arrayMax > pixvalPositiveMin) and (pixvalMax > 0):
                        array[0, :, :] = np.power((np.log(array[0, :, :].clip(pixvalPositiveMin, pixvalMax))\
                          - np.log(pixvalPositiveMin)) \
                          / (np.log(pixvalMax) - np.log(pixvalPositiveMin)),
                          (1.0 / gamma))
                except:
                    pass
# Save the array to StringIO object and Open it with PIL
            f1 = cStringIO.StringIO()
            imsave(f1, array[0, :, :], cmap=colormap, format="png")
            f1.reset()
            pilImgFig = Image.open(f1)
            print "(%3.1f sec) " % (time.clock() - tic)

# Create a color bar, save it to StringIO object and Opne it with PIL
            f2 = cStringIO.StringIO()
            bar = outer(np.ones(colorbar_height), np.arange(0, rasterXSize, 1))
            imsave(f2, bar, cmap=colormap, format="png")
            f2.reset()
            pilImgCbar = Image.open(f2)

# Create a new PIL canvas for the color bar
            imgWidth = int(2*margin + rasterXSize)
            imgHeight = colorbar_height*4
            pilImgCbarfig = Image.new('RGB', (imgWidth, imgHeight), "white")
# -- put the colorbar PILimage
            box = (margin, 1)
            pilImgCbarfig.paste(pilImgCbar, box)
            draw = ImageDraw.Draw(pilImgCbarfig)

# -- Add scales into the canvas
# ---- create an array which shows the locations of the values in the colorbar
            scaleTextLocation = np.linspace(pow(10.0, -10), 0.95, numOfTicks)

# ---- Set values of the scale in the colorbar
            if logarithm:
                scaleArray = np.exp(np.power(scaleTextLocation, gamma)\
                             * (np.log(pixvalMax) - np.log(pixvalPositiveMin))\
                             + np.log(pixvalPositiveMin))
            else:
                scaleArray = scaleTextLocation * (pixvalMax - pixvalMin) + pixvalMin

# ---- modify the number of digits of scaleArray for showing
# ------ Get number of digits of the 2nd value in the scaleArray
            scale_temp = scaleArray[1]
            pixvalDecimals = 0
            if abs(scale_temp) > 1:
                while abs(scale_temp) > 1.0:
                    pixvalDecimals += 1
                    scale_temp /= 10
            else:
                while abs(scale_temp) < 1.0:
                    pixvalDecimals -= 1
                    scale_temp *= 10

# ------ Get a number of digits of the scale
            barScale = (pixvalMax - pixvalMin) / (numOfTicks - 1)
            barDecimals = 0
            if barScale > 1:
                while barScale > 1.0:
                    barDecimals += 1
                    barScale /= 10
            else:
                while barScale < 0.1:
                    barDecimals -= 1
                    barScale *= 10

# ------ Round values to proper style
            scaleArray /= pow(10.0, pixvalDecimals)
            if barDecimals < pixvalDecimals:
                barDecimals = (pixvalDecimals - barDecimals) + 2
            else:
                barDecimals = 2
            scaleArray = np.round(scaleArray, barDecimals)

# ------ Convart the scale array to a string list
            scaleText = list(map(str, scaleArray))

# ----  check the scale
            ##print "scaleText : ", scaleText
# ---- draw lines on the color bar and the scales
            for i in range(numOfTicks):
                coordX = int(margin+scaleTextLocation[i]*rasterXSize)
                box = (coordX, 1, coordX, colorbar_height)
                draw.line(box, fill='rgb(0,0,0)')
                box = (coordX-2, int(colorbar_height*1.3))
                draw.text(box, scaleText[i], fill='rgb(0,0,0)')
            if pixvalDecimals != 0:
                    box = (coordX, int(colorbar_height*2.3))
                    draw.text(box, "x10^%d" %pixvalDecimals, fill='rgb(0,0,0)')

# create a new canvas for the output file
        imgWidth = int(2*margin + rasterXSize)
        if len(bands)==3:
            imgHeight = int(2*margin + rasterYSize + pad)
        else:
            imgHeight = int(2*margin + rasterYSize + pad + pilImgCbarfig.size[1])
        pilImg = Image.new('RGB', (imgWidth, imgHeight), "white")

# paste the figure
        box = (margin, margin)
        pilImg.paste(pilImgFig, box)

# paste the color bar
        if len(bands)==1:
            box = (0, margin + rasterYSize + pad)
            pilImg.paste(pilImgCbarfig, box)

# write text
# -- set fonts
        fileName_font = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'fonts/times.ttf')
        font = ImageFont.truetype(fileName_font, fontSize)
        # -- write text each line into pilImgTxt
        if textwidthMax > rasterXSize:
            textwidthMax = rasterXSize
        pilImgTxt = Image.new('L', (textwidthMax, pad-10), "white")
        pilImgDraw = ImageDraw.Draw(pilImgTxt)
        textWidth = 0
        textHeight = 0
        for line in legend.splitlines():
            ext = pilImgDraw.textsize(line, font)
            pilImgDraw.text((0, textHeight), line, font=font, fill= 0)
            textWidth = max(ext[0], textWidth)
            textHeight += ext[1]
        # paste pilImgTxt into pilImg
        pilImg.paste(pilImgTxt, (margin, margin+rasterYSize+10))
        # save the file
        pilImg.save(fileName + "." + extension)

        if len(bands)==1:
            f1.close()
            f2.close()


    """
    def write_figure(self, fileName, bandNo=1,
                     pixelValMin=None, pixelValMax=None,
                     imageDatatype=None, thresholdRatio=1.0,
                     useFullMatrix=False, extension='png',
                     useImsave=False):
        '''Save a raster band to a figure in grapfical format.

        Get numpy array from the band specified either by given band
        number or band id adjust the array brightness and contrast
        using the given min/max or histogram ratio write to file.

        Parameters
        ----------
            fileName: string
                Output file name
            bandNo: int
            bandName: a list, optional
                (e.g.: bandIdList = {"name":"radiance", "wavelength":"645"})
            thresholdRatio: float (0.0 - 1.0), optional
                e.g. : thresholdRatio = 0.95 means to round off 5%
                        form the both sides (upper and lower sides).
            useFullMatrix: boolean, optional
                if true, the full matrix is used for estimating min/max,
                otherwise only image scaled down to 100x100 (MUCH FASTER)

        Raises
        ------
            DataError: occurs when the array of the band is empty

        '''
        # read NumPy array from band
        tic = time.clock()
        print "Writing figure "
        rawArray = self[bandNo]

        if rawArray is None:
            raise DataError("Nansat.write_figure(): "
                            "array of the band is empty")

        #if easy imsave operation is allowed
        if useImsave:
            plt.imsave(fileName + ".png", rawArray)
            print "(%3.1f sec) " % (time.clock() - tic),
            return

        # if value < pixelValMin then replace as value = pixelValMin
        # if value > pixelValMax then replace as value = pixelValMax
        if pixelValMin is None:

            # reduce input matrix to the size 100 x 100 for calculating
            # histogram
            if not useFullMatrix:
                step1 = max(rawArray.shape[0] / 100, 1)
                step2 = max(rawArray.shape[1] / 100, 1)
                histArray = rawArray[::step1, ::step2]
            else:
                histArray = rawArray

            # get minmax from histogram analysis
            pixelValMin, pixelValMax = self._get_pixelValueRange(
                                        histArray, thresholdRatio)
        print "[%f %f]" % (pixelValMin, pixelValMax)
        print "(%3.1f sec) " % (time.clock() - tic)

        # cut away values over limits and save to a PNG
        np.clip(rawArray, pixelValMin, pixelValMax, out=rawArray)
        toimage(rawArray).save(fileName + "." + extension)
        print "(%3.1f sec) " % (time.clock() - tic)
    """
    def get_domain(self):
        ''' Returns: Domain of the Nansat object '''
        return Domain(self.vrt)

    def _get_mapper(self, mapperName, bandList):
        '''Creare VRT file in memory (VSI-file) with variable mapping

        Create a mapperList based on all mappers in the subdir 'mappers'.
        If mapperName is given, it is the first in the mapperList.
        Loop over all availble mappers to get the matching one.
        In the loop:
            If the specific error appears the mapper is not used
            and the next mapper is tested.
            Otherwise the mapper returns VRT.
        If type of the sensor is identified, add mapping variables.
        If all mapper do not fit, simply copy the input DS into a VSI/VRT

        Parameters
        ----------
        mapperName : string, optional
            "ASAR", "hurlam", "merisL1", "merisL2", "ncep", "radarsat2",
            "seawifsL2" are currently available.  (27.01.2012)
        bandList : list, optional
            band numbers to fetch.
            If None is given, all bands in the file are fetched.

        Returns
        -------
            vsiDataset : VRT dataset
                VRT dataset with mapping of variables keeped in memory)

        Raises
        --------
            TypeError: occurs when the given driver type is not registarated
                        in the mappers.

        '''
        # create a mapper list based on the files in the folder "mappers"
        nansatDir = path.dirname(path.realpath(__file__))

        allMapperFiles = listdir(path.join(nansatDir, "winmappers"))
        allMapperFiles = fnmatch.filter(allMapperFiles, 'mapper_*.py')

        # add the given mapper first
        mapperList = ['mapper_' + mapperName]

        # loop through appropriate files and add to the list
        for iFile in allMapperFiles:
            iFile = iFile.replace(".py", "")
            mapperList.append(iFile)

        # try to add path for windows, add for linux otherwise
        try:
            sys.path.append(path.join(unicode(nansatDir, "mbcs"),
                            "mappers"))
        except:
            sys.path.append(path.join(nansatDir, "mappers"))

        # try to import and get VRT datasaet from all mappers. Break on success
        # if none of the mappers worked - None is returned
        vrtDataset = None
        for iMapper in mapperList:
            try:
                mapper_module = __import__(iMapper)
                vrtDataset = mapper_module.Mapper(self.rawVRTFileName,
                                         self.fileName, self.dataset,
                                         self.metadata,
                                         bandList).vsiDataset
                break
            except:
                pass

        # if no mapper fits, make simple copy of the input DS into a VSI/VRT
        if vrtDataset is None:
            print 'No mapper fits!'
            vrtDataset = self.vrtDriver.CreateCopy(self.rawVRTFileName,
                                                   self.dataset)

        return vrtDataset
    """
    def _get_pixelValueRange(self, array, ratio):
        '''Get proper pixel value range for writing a figure in PNG

        Return a proper pixel value range (cmin, cmax)
        to wrige a figure with toimage.
        the argument "ratio" is used to specify the threshold of a pixel value
        that should be counted.

        Parameters
        ----
            array : numpy array
                array of a band
            ratio  : float (0.0 - 1.0)
                1-ratio means round off (1-ratio)x100 %
                form both upper and lower sides

        Returns
        -------
            edge_min : float
                minimum threshold of the pixel value
            edge_max : float
                maximum threshold of the pixel value

        '''
        # exclude zeros from array (wich spoil histo)
        array.flatten()
        array = array[array != 0]

        # try to make histogram
        tic = time.clock()
        try:
            hist, lowerreallimit, binsize,
            extrapoint = cumfreq(array, numbins=15)
        except:
            hist = None

        if hist is None:
            edge_min = 0
            edge_max = 1
        else:
            toc = time.clock()
            hist_eq = hist / max(hist)
            hist_min = hist_eq[hist_eq < (1 - ratio)]
            hist_max = hist_eq[hist_eq > ratio]

            if len(hist_min) == len(hist_eq):
                edge_min = lowerreallimit + (len(hist_eq) - 1.5) * binsize
            elif len(hist_min) == 0:
                edge_min = lowerreallimit + 0.5 * binsize
            else:
                edge_min = lowerreallimit + (len(hist_min) - 0.5) * binsize

            if len(hist_max) == len(hist_eq):
                edge_max = lowerreallimit + (1.0 + 0.5) * binsize
            elif len(hist_eq) == 0:
                edge_max = lowerreallimit + (len(hist_eq) - 0.5) * binsize
            else:
                edge_max = lowerreallimit + \
                         (len(hist_eq) - len(hist_max) + 0.5) * binsize

        return edge_min, edge_max
    """
    def _modify_warpedVRT(self, rawWarpedVRT,
                          rasterXSize, rasterYSize, geoTransform):
        ''' Modify rasterXsize, rasterYsize and geotranforms in the warped VRT

        Parameters
        ----------
            rasterXSize: integer
                desired X size of warped image
            rasterYSize: integer
                desired Y size of warped image
            rasterYSize: tuple of 6 integers
                desired GeoTransform size of the warped image

        Modifies
        --------
            the VRT file which keepes warped vrt is modified

        '''
        # Write the warpedVRT to a VSI-file
        vrtDatasetCopy = self.vrtDriver.CreateCopy(self.warpedVRTFileName,
                                                   rawWarpedVRT)
        # Get XML content from VSI-file
        vsiFileContent = self._read_write_vsi_file(self.warpedVRTFileName)

        # Get element from the XML content and modify some elements
        # using Domain object parameters
        element = XML(vsiFileContent)
        element.set("rasterXSize", str(rasterXSize))
        element.set("rasterYSize", str(rasterYSize))
        tree = ElementTree(element)

        # convert proper string style and set to the GeoTransform element
        geoTransformString = str(geoTransform).\
                        translate(maketrans("", ""), "()")

        # replace GeoTranform
        elem = tree.find("GeoTransform")
        elem.text = geoTransformString
        # replace DstGeoTranform
        elem = tree.find("GDALWarpOptions/Transformer/"
                         "GenImgProjTransformer/DstGeoTransform")
        elem.text = geoTransformString

        elem = tree.find("GDALWarpOptions/Transformer/"
                         "GenImgProjTransformer/DstInvGeoTransform")
        # get inverse geotransform
        invGeotransform = gdal.InvGeoTransform(geoTransform)
        # convert proper string style and set to the DstInvGeoTransform element
        elem.text = str(invGeotransform[1]).\
                        translate(maketrans("", ""), "()")

        # Overwrite element
        element = tree.getroot()

        # Write the modified elemements into VSI-file
        self._read_write_vsi_file(self.warpedVRTFileName, tostring(element))

        newWarpedVRT = gdal.Open(self.warpedVRTFileName)
        return newWarpedVRT

    def _read_write_vsi_file(self, vsiFileName, vsiFileContent=None):
        '''Read or write content of a VSI-file

        If only file name is given then content of the file is read and
        returned, if content is provided then it is written to the file

        Parameters:
            vsiFileName: string
                Name of the VSI file to read/create
            vsiFileContent: string, optional
                Content of the VSI file to write

        Returns:
            vsiFileContent: string
                Content which is read from the VSI file
        '''

        if vsiFileContent is None:
            #read from the vsi-file
            # open
            vsiFile = gdal.VSIFOpenL(vsiFileName, "r")
            # get file size
            gdal.VSIFSeekL(vsiFile, 0, 2)
            vsiFileSize = gdal.VSIFTellL(vsiFile)
             # fseek to start again
            gdal.VSIFSeekL(vsiFile, 0, 0)
            # read
            vsiFileContent = gdal.VSIFReadL(vsiFileSize, 1, vsiFile)
            gdal.VSIFCloseL(vsiFile)
            return vsiFileContent
        else:
            #write to the vsi-file
            vsiFile = gdal.VSIFOpenL(vsiFileName, 'w')
            gdal.VSIFWriteL(vsiFileContent,
                            len(vsiFileContent), 1, vsiFile)
            gdal.VSIFCloseL(vsiFile)
            return 0

    def _specify_bandNo(self, bandID):
        '''Specify a band number based on bandID {'key_name': 'key_value'}

        Compare the key values of the bandID to the values of the
        metadata in each band.
        If matched, append the band number (iBand) into candidate list.
        If not, go to the next band.
        Iterate these steps until all bands are checked.
        If single band is specified at the end, return the band number.
        Otherwise raise OptionError.

        Parameters
        ----------
            bandID: a dictionary
                Parameters to specify single band
                (e.g. {"name":"radiance", "wavelength":"1234"})

        Returns
        -------
            candidate[0]+1 : a band number

        Raises
        ------
            OptionError: occurs when there is no band which satisfies
            the condition (bandID) or when there are several bands chosen
            by the condition.

        '''
        # search for the specific band based on bandID
        candidate = []
        # loop through all bands
        for iBand in range(self.vrt.RasterCount):
            # get metadata from band
            bandMetadata = self.rawVRT.GetRasterBand(iBand + 1).GetMetadata()
            allKeysAreGood = True
            # loop through all input keys
            for bandIDKey in bandID:
                # if band doesn't have key from input or value doesn't match
                if (bandIDKey not in bandMetadata or
                        bandMetadata[bandIDKey] != bandID[bandIDKey]):
                    allKeysAreGood = False

            if allKeysAreGood:
                candidate.append(iBand)

        # if a band is specified, set it to bandNo.
        # if some bands are chosen, give an error message and stop.
        if len(candidate) == 1:
            print "You chose bandNo:", candidate[0] + 1
            return candidate[0] + 1
        elif len(candidate) >= 2:
            raise OptionError("Nansat._specify_bandNo(): "
                              "Cannot specify a single band "
                              "by the given arguments")
        else:
            raise OptionError("Nansat._specify_bandNo(): "
                              "Cannot find any band by the given arguments")
