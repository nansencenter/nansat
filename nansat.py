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

from os import path, listdir

from string import maketrans, ascii_uppercase, digits
import sys
import time
from random import choice

import fnmatch
import Image
import ImageDraw
import ImageFont
import ImageOps
import matplotlib.cm as cm
import numpy as np
from math import log10
from math import floor
from scipy.misc import toimage, pilutil
from scipy.stats import cumfreq
from xml.etree.ElementTree import XML, ElementTree, tostring

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
    '''Error for improper options (arguments) '''
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
        self.mapperList: list of file names
            list of available working mappers

        self.rawVRTFileName : file name
            set '/vsimem/vsiFile.vrt'
        self.warpedVRTFileName : file name
            set '/vsimem/vsi_warped.vrt'
        self.vrtDriver : VRT driver
            set GDAL VRT driver

        self.fileName : file name
            set file name given by the argument
        self.dataset : GDAL dataset
            set GDAL dataset
        self.metadata : metadata
            set metadata of the dataset
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
        # SET SOME CONSTANTS
        # all available mappers
        self.mapperList = [
              'mapper_ASAR.py',
              'mapper_hirlam.py',
              'mapper_merisL1.py',
              'mapper_merisL2.py',
              'mapper_modisL1.py',
              'mapper_ncep.py',
              'mapper_radarsat2.py',
              'mapper_seawifsL2.py',
              #'mapper_MOD44W.py',
              ]

        # names of raw and warped VRT files in memory
        # make random string and append to fileNames
        allChars=ascii_uppercase + digits
        randomChars = ''.join(choice(allChars) for x in range(6))
        self.rawVRTFileName = '/vsimem/rawVRT_%s.vrt' % randomChars
        self.warpedVRTFileName = '/vsimem/warpedVRT_%s.vrt' % randomChars
        # VRT driver
        self.vrtDriver = gdal.GetDriverByName("VRT")

        # set input file name
        self.fileName = fileName

        # set input GDAL dataset
        self.dataset = gdal.Open(self.fileName)
        if (self.dataset is None) or (self.dataset == ""):
            raise GDALError("Nansat._init_(): Cannot get the dataset from "
                            + self.fileName)
        # metadata
        self.metadata = self.dataset.GetMetadata()
        if (self.metadata is None) or (self.metadata == ""):
            raise GDALError("Nansat._init_(): Cannot get the metdadata")

        # get VRT with mapping of variables
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

    def watermask(self, mod44path=None):
        '''Create numpy array with watermask (water=1, land=0)

        250 meters resolution watermask from MODIS 44W Product:
        http://www.glcf.umd.edu/data/watermask/

        Watermask is stored as tiles in TIF(LZW) format and a VRT file
        All files are stored in one directory.
        A tarball with compressed TIF and VRT files should be additionally
        downloaded from the Nansat wiki:
        https://svn.nersc.no/nansat/wiki/Nansat/Data/Watermask

        The method:
            Gets the directory either from input parameter or from environment
            variable MOD44WPATH
            Open Nansat object from the VRT file
            Reprojects the watermask onto the current object using reproject()
            or reproject_on_jcps()
            Returns the reprojected Nansat object

        Parameters:
        -----------
            mod44path : string, optional, default=None
                path with MOD44W Products and a VRT file

        Returns:
        --------
            watermask : Nansat object with water mask in current projection

        See also:
        ---------
            250 meters resolution watermask from MODIS 44W Product:
            http://www.glcf.umd.edu/data/watermask/
        '''
        mod44DataExist = True
        # check if path is given in input param or in environment
        if mod44path is None:
            mod44path = os.getenv('MOD44WPATH')
        if mod44path is None:
            mod44DataExist = False
        # check if VRT file exist
        elif not os.path.exists(mod44path+'/MOD44W.vrt'):
            mod44DataExist = False

        if not mod44DataExist:
            # MOD44W data does not exist generate empty matrix
            watermask = np.zeros(self.vrt.RasterXSize, self.vrt.RasterYSize)
        else:
            # MOD44W data does exist: open the VRT file in Nansat
            watermask = Nansat(mod44path+'/MOD44W.vrt');
            # choose reprojection method
            if self.vrt.GetGCPCount() == 0:
                watermask.reproject(Domain(self.vrt))
            else:
                watermask.reproject_on_gcp(self)

        return watermask

    def write_figure(self, fileName, bands=1, colormapName = "jet",
                     clim = [None, None], ratio=1.0,
                     logarithm=False, gamma=2,
                     numOfTicks=5, fontSize=10,
                     margin=30, pad=60, textWidthMax=600,
                     colorbarHeight=10, barFontSize=8,
                     putLegend = False,
                     titleString = "",
                     legenHeight = 50,
                     extension='png'):

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
            ratio : listimageDatatype, optional
                ratio of number of the pixels
                used to compute cmin and cmax.
            logarithm : boolean, optional
                If True, tone curve is used to convert pixel values.
                If False, linear.
            gamma : float ( > 0.0), optional, default = 2.0
                coefficient of the tone curve
            numOfTicks : int, optional, default = 5
                number of ticks
            fontSize : int, optional, default = 10
                size of font for the legend
            margin: int, optional
                margin of the output file
            pad : int, optional
                height for the legend
            textWidthMax : int, optional
                width for the legend
            colorbarHeight : int, optional
                height for the colorbar
            barFontSize : int, optional, default = 8
                size of font for the colorbar
            titleString : text, optional
            extension : extension, optional
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
        
        # convert band and ratio from a number to list
        if not isinstance(band, list):
            band = [band];
        if not isinstance(ratio, list):
            ratio = [ratio];

        # check if bands is proper (1 or 3)
        if not ((len(bands) == 3) or (len(bands) == 1)):
            raise OptionError("write_figure(): number of elements of bands list must be 1 or 3")

        # check if ratio is proper (1 or 3)
        if not ((len(ratio) == 3) or (len(ratio) == 1)):
            raise OptionError("write_figure(): number of elements of ratio list must be 1 or 3")
        
        #ANTON: TOO COMPLEX IF'S ARE REPLACED BY SIMPLE ONES
        #       AFTER ERROR IS RAISED - NO NEED TO CHECK ELIF
        #       TOO MANY ERRORS ARE RAISED - WE SHOULD GIVE WARNING AND REPLACE
        #       BAD VALUES WITH REASONABLE ONES
        
        # populate list of ratios with the same value if necessary
        if len(ratio) != len(bands):
            ratio = [ratio[0]] * len(bands)
                
        # check if ratios are within possible values
        if any(val > 1.0 for val in ratio) or any(val < 0.0 for val in ratio):
            raise OptionError("write_figure(): elements of ratio list take [0.0 - 1.0]")

        # check if clim consists of two elements
        if len(clim) != 2:
            raise OptionError("write_figure(): clim length should be 2")
        
        # check each element of clim
        # convert [c1, c2] to [[c1], [c2]]
        for i in [0:1]:
            # convert a number to list
            if not isinstance(clim[i], list):
                clim[i] = [clim[i]]
            # compare with length of bands and populate with same value (if not)
            if len(clim[i]) != len(bands):
                clim[i] = clim[i][0] * len(bands)
        
        #ANTON: THE BLOCK BELOW SHOULD BE REFINED
        #       IT IS NOT OBVIOUS WHY cmin.append()... after exept        
        
        # generate vectors with limits of colors
        cmin = []
        cmax = []
        for i in range(len(bands)):
            try:
                val = clim[1][i] - clim[0][i]
                if val>0:
                    cmin.append(clim[0][i])
                    cmax.append(clim[1][i])
                else:
                    raise OptionError("write_figure(): clim is given as [cmin, cmax] or [(cmin1, cmin2, cmin3), (cmax1, cmax2, cmax3)]")
            except:
                cmin.append(clim[0][i])
                cmax.append(clim[1][i])

        if gamma <= 0.0:
            raise OptionError("gamma argument must be a positive float")

        rasterXSize = self.vrt.RasterXSize
        rasterYSize = self.vrt.RasterYSize

        print "Reading figure(s)  size = (", rasterXSize, "x", rasterYSize, ")"
        array = np.array([])

        # read array from the data
        tic = time.clock()
        for iBand in bands:
            iArray = self[iBand]
            if iArray is None:
                raise DataError("Nansat.write_figure(): "
                                "array of the band is empty")
            # if three bands are given,
            # append the second and third arrays to the previous array.
            array = np.append(array, iArray)
        # reshape the array for each band
        array = array.reshape(len(bands), rasterYSize, rasterXSize)
        
        #ANTON: BELOW LINES UNTIL 891 ARE PREFOrmed ONLY IF CMIN OR CMAX ARE NONE
        
        pixvalMinList = []
        pixvalMaxList = []
        histMinList = []
        histMaxList = []
        histBuckets = []
        histScale = []
        histFreq = np.array([])

        # get band infomation
        for i, iBand in enumerate(bands):
            bandPixvalueInfo = self.vrt.GetRasterBand(iBand).\
                               GetDefaultHistogram()
            # minimum pixel value
            histMinList.append(float(bandPixvalueInfo[0]))
            # maximum pixel value
            histMaxList.append(float(bandPixvalueInfo[1]))
            # number of widths of histogram (256?)
            histBuckets.append(float(bandPixvalueInfo[2]))
            # average variation of pixel value for each width
            histScale.append((histMaxList[i] - histMinList[i]) /
                                    histBuckets[i])
            # frequency of each interval (histogram)
            histFreq = np.append(histFreq, (np.array(bandPixvalueInfo[3],
                                                         dtype=float)))
            bandStatistics = self.vrt.GetRasterBand(iBand).GetStatistics(True, True)
            # minimum pixel value
            pixvalMinList.append(float(bandStatistics[0]))
            # maximum pixel value
            pixvalMaxList.append(float(bandStatistics[1]))

        # reshape pixelHist array for each band
        histFreq = histFreq.reshape(len(bands), histBuckets[0])

        # <OPTION> Get cmin and cmax and clip the pixel values
        numOfPixels = rasterXSize*rasterYSize
        for iBand in range(len(bands)):
            pixvalMin = pixvalMinList[iBand]
            pixvalMax = pixvalMaxList[iBand]
            # if ratio is not 1.0, compute cmin and cmax based on the ratio
            if ratio[iBand] != 1.0:
                # percentage of accumulated histogram of iBand
                accumulateHist = np.add.accumulate(histFreq[iBand, :]) / \
                                                   numOfPixels
                # if the frequency (ratio) of the first inverval is
                # larger than the given ratio (1-ratio) and
                # the frequency of the second inverval of the histogram
                # is zero (it means accumulateHist[0] == accumulateHist[1]),
                # cmin is replaced the second minimum pixel value.
                # In short, if there are many black pixels (no value pixels),
                # set the second minimum value to cmin.
                if (accumulateHist[0] > (1-ratio[iBand])) and \
                    (accumulateHist[0] == accumulateHist[1]):
                    accumulateHist = accumulateHist[accumulateHist >
                                                    accumulateHist[0]]
                    cminRatio = histMaxList[iBand] - \
                                  histScale[iBand] * len(accumulateHist)
                    cmaxRatio = pixvalMaxList[iBand]
                # if the frequency (ratio) of the first inverval is
                # less than the given ratio (1-ratio) and
                # the frequency of the second inverval is zero,
                # firstly remove the minimum value and secondly adjust
                # cmin and cmax from both sides (second minimum and
                # maximum values) to satisfy the given ratio
                elif accumulateHist[0] == accumulateHist[1]:
                    histMin = accumulateHist[accumulateHist <
                                             (1 - ratio[iBand] +
                                             accumulateHist[0]) / 2]
                    histMax = accumulateHist[accumulateHist <
                                             (1 + ratio[iBand] +
                                             accumulateHist[0]) / 2]
                    cminRatio = histMinList[iBand] + \
                                  histScale[iBand] * len(histMin)
                    cmaxRatio = histMinList[iBand] + \
                                  histScale[iBand] * len(histMax)
                # if the frequency (ratio) of the last inverval is
                # larger than the given ratio and
                # the frequency of the second to last inverval of the histogram
                # is zero (it means accumulateHist[-2] == accumulateHist[-3]),
                # cmax is replaced the second maximum pixel value.
                # In short, if there are many black pixels (no value pixels),
                # set the second maximum value to cmax.
                elif (accumulateHist[-2] < (ratio[iBand])) and \
                    (accumulateHist[-2] == accumulateHist[-3]):
                    accumulateHist = accumulateHist[accumulateHist <
                                                    accumulateHist[-2]]
                    cminRatio = pixvalMinList[iBand]
                    cmaxRatio = histMinList[iBand] + \
                                  histScale[iBand] * len(accumulateHist)
                # if the frequency (ratio) of the first inverval is
                # less than the given ratio (1-ratio) and
                # the frequency of the second inverval is zero,
                # firstly remove the minimum value and secondly adjust
                # cmin and cmax from both sides (second minimum and
                # maximum values) to satisfy the given ratio
                elif accumulateHist[-2] == accumulateHist[-3]:
                    histMin = accumulateHist[accumulateHist <
                                             (accumulateHist[-2] -
                                             ratio[iBand]) / 2 ]
                    histMax = accumulateHist[accumulateHist <
                                             (accumulateHist[-2] +
                                             ratio[iBand]) / 2 ]
                    cminRatio = histMinList[iBand] + \
                                  histScale[iBand] * len(histMin)
                    cmaxRatio = histMinList[iBand] + \
                                  histScale[iBand] * len(histMax)
                # specify cmin and cmax by removing (1 - ratio) / 2 equally
                # from both sides (minimum and maximum values)
                # to satisfy the given ratio
                else:
                    histMin = accumulateHist[accumulateHist <
                                             ((1 - ratio[iBand]) / 2)]
                    histMax = accumulateHist[accumulateHist >
                                             (1 - ((1 - ratio[iBand]) / 2))]
                    cminRatio = histMinList[iBand] + \
                                  histScale[iBand] * len(histMin)
                    cmaxRatio = histMinList[iBand] + \
                                  histScale[iBand] * \
                                  (len(accumulateHist) - len(histMax))
                if cmin[iBand] is None or cmin[iBand] < cminRatio:
                    cmin[iBand] = cminRatio

                if cmax[iBand] is None or cmax[iBand] > cmaxRatio:
                    cmax[iBand] = cmaxRatio

            # if cmin and cmax are given, clip the array
            if (cmax[iBand] or cmin[iBand]) is not None:
                if (cmax[iBand] is not None) and \
                    (cmax[iBand] < pixvalMax) and (cmax[iBand] > pixvalMin):
                    pixvalMax = cmax[iBand]
                if (cmin[iBand] is not None) and \
                    (cmin[iBand] > pixvalMin) and (cmin[iBand] < pixvalMax):
                    pixvalMin = cmin[iBand]
        
        #ANTON: IT SHOULD NOT BE CLIPPED BUT RATHER CONVERTED TO UINT8
        #       (TO KEEP MEM LOW) AND NOT TO DUPLICATE CODE
        
                array[iBand, :, :] = np.clip(array[iBand, :, :],
                                             pixvalMin, pixvalMax)
        print "903 : ", pixvalMin, " -- ",pixvalMax

        # If three bands are given, marge them and crate a PIL image
        if len(bands)==3:
            # Normalize RGB arrays to the interval [0,255]
            for iBand in range(3):
                if array[iBand, :, :].max() == array[iBand, :, :].min():
                    bandColor = {0: 'red', 1:'green', 2:'blue'}[iBand]

                    raise OptionError("min. and max. pixel valuses in " +
                                      bandColor + " band are same." +
                                      " check ratio, cmin and cmax!!")
                else:
                    array[iBand, :, :] = (array[iBand, :, :] -
                                          array[iBand, :, :].min()) * 255 / \
                                         (array[iBand, :, :].max()-
                                          array[iBand, :, :].min())
        
        #ANTON: NOT A GOOD IDEA
        #           1. RETURN SHOULD BE ONLY ONE
        #           2. WE STILL MIGH WANT TO ADD LEGEND (WITHOUT COLORBAR BUT WITH TITLE)
        
            # convert RGB arrays to PIL image
            # merge three bands and save it to the file
            Image.merge("RGB", (Image.fromarray(np.uint8(array[0, :, :])),
                                Image.fromarray(np.uint8(array[1, :, :])),
                                Image.fromarray(np.uint8(array[2, :, :])))).\
                                save(fileName + "." + extension)
            return

        # If one band is given, create its PIL image and color bar PIL image
        else:
            # <OPTION> Convart array based on logarithmic tone curve
            if logarithm:
                array[0, :, :] = (np.power((array[0, :, :] - pixvalMin) /
                                 (pixvalMax - pixvalMin), (1.0 / gamma)))* \
                                 (pixvalMax - pixvalMin) + pixvalMin

            # create a color palette based on cmapName
            myPalette = self._create_mypalette(colormapName)
        
        #ANTON: 253 SHOULD NOT BE HARDCODED. WHAT HAPPENS WHEN WE WANT TO ADD 
        #       LIGHT GRAY FOR CLOUDS AND DARK GRAY FOR LAND?
        
            # read array with PIL image and set the palette
            array[0, :, :] = (array[0, :, :] -pixvalMin)*253 / \
                             (pixvalMax - pixvalMin)

            if putLegend is False:
                # save the file
                pilImgFig = Image.fromarray(np.uint8(array[0, :, :]))
                pilImgFig.putpalette(myPalette)
                pilImgFig.save(fileName + "." + extension)
                return

            else:
                # create a color bar and set the palette
                array255 = np.ones((legenHeight+colorbarHeight*3,
                                    rasterXSize)) * 255
                pilImgFig = Image.fromarray(np.uint8(np.append(array[0, :, :],
                                            array255, 0)))
                pilImgFig.putpalette(myPalette)

                bar = np.outer(np.ones(colorbarHeight),
                           np.linspace(0, 253, int(rasterXSize*0.8)))
                pilImgCbar = Image.fromarray(np.uint8(bar))
                pilImgCbar.putpalette(myPalette)

                pilImgFig.paste(pilImgCbar,
                                   (int(rasterXSize*0.1),
                                   rasterYSize+legenHeight))

                draw = ImageDraw.Draw(pilImgFig)

                # add scales into the colorbar canvas
                # create an array which indicates the scale values

                # compute the locations of scaleText
                scaleTextLocation = np.linspace(0, 0.95, numOfTicks)
                if logarithm:
                    # create an array for valuse on the colorbar
                    scaleArray = np.array([])
                    for i in range(numOfTicks):
                        if i == 0:
                            scaleArray = np.append(scaleArray, [pixvalMin], 0)
                        else:
                            scaleArray = np.append(scaleArray,
                            [np.power(scaleTextLocation[i], (1.0 / gamma)) *\
                                 (pixvalMax - pixvalMin) + pixvalMin], 0)
                else:
                    scaleArray = scaleTextLocation * (pixvalMax - pixvalMin) \
                                 + pixvalMin

                # specify the description form
                formatList = map(self._create_scaleFormat, scaleArray)
                scaleText = map(self._round_Scale, scaleArray, formatList)
                ##print "1014: ", scaleArrayDigit
                ##print "1015 : ", barScaleDigit, scaleArray, formatList
                ##print "1023 : ", scaleText

                # set fonts for colorbar
                fileName_font = os.path.join(os.path.dirname(
                                             os.path.realpath(__file__)),
                                             'fonts/times.ttf')
                font = ImageFont.truetype(fileName_font, barFontSize)
                # draw lines on the color bar and write the scales
                for i in range(len(scaleArray)):
                    coordX = int(scaleTextLocation[i]*rasterXSize*0.8 +
                                 rasterXSize*0.1)
                    box = (coordX, rasterYSize+legenHeight, coordX,
                           rasterYSize+legenHeight+colorbarHeight-1)
                    draw.line(box, fill= 254)
                    box = (coordX-5, rasterYSize+legenHeight+colorbarHeight)
                    draw.text(box, scaleText[i], font=font, fill= 254)

                # set font size for text
                font = ImageFont.truetype(fileName_font, fontSize)
                longName = self.vrt.GetRasterBand(bands[0]).GetMetadataItem("long_name")
                if longName is None:
                    longName = "no longName"
                units = self.vrt.GetRasterBand(bands[0]).GetMetadataItem("units")
                if units is None:
                    units = "no units"
                caption = longName + " / " + units
                box = (5, rasterYSize+ int(legenHeight*0.7))
                draw.text(box, str(caption), font=font, fill= 254)

                if titleString != "":
                    # write text each line onto pilImgCanvas
                    if textWidthMax > rasterXSize:
                        textWidthMax = rasterXSize
                    textWidth = 0
                    textHeight = rasterYSize + 5
                    for line in titleString.splitlines():
                        text = draw.textsize(line, font)
                        draw.text((0, textHeight), line, font=font, fill= 0)
                        textWidth = max(text[0], textWidth)
                        textHeight += text[1]
                pilImgFig.save(fileName + "." + extension)

    def get_domain(self):
        ''' Returns: Domain of the Nansat object '''
        return Domain(self.vrt)

    def _create_list(self, r, g, b):
        '''Create a list from numbers

        Parameters
        ----------
            r, g, b: float / int

        Returns : list

        '''
        return [r, g, b]


    def _create_scaleFormat(self, value):
        '''Return writing format for scale on the colorbar

        Parameters
        ----------
            values: int / float / exponential

        Returns
        --------
            string

        '''
        if value==0:
            return "%d"
        else:
            
            #ANTON: _get_digit SHOULD NOT BE A METHOD IT IS CALLED ONLY ONCE AND VERY SIMPLE
            
            digit = self._get_digit(value)
            if digit==0:
                return "%4.2f"
            elif digit==1:
                return "%4.1f"
            elif digit==2:
                return "%d"
            elif digit==-1:
                return "%3.1f"
            elif digit==-2:
                return "%4.2f"
            else:
                return "%4.2e"

    def _create_mypalette(self, cmapName, text=False):
        '''Create a palette based on colormap name.

        254 colors are created based on the given colormap name.
        255th color is black and 256th is white.

        Parameters
        ----------
        cmapName : string
            matplotlib colormap name

        Returns
        -------
            palette : a list
                sequence of (256*3) integer

        Raise
        ------
            OptionError : occures when cmapName is not in the matplotlib colormap

        See Also
        --------
        http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        '''
        try:
            colorDic = cm.datad[cmapName]
        except:
            raise OptionError("cmapName is not proper.")

        colorList = [colorDic['red'], colorDic['green'], colorDic['blue']]
        lut = []

        for iColor in range(len(colorList)):
            lut.append([])
            for i in range(len(colorList[iColor])-1):
                spaceNum = int (254 * (colorList[iColor][i+1][0] -
                                       colorList[iColor][i][0]))
                iColorArray = np.array(np.linspace(colorList[iColor][i][2],
                                          colorList[iColor][i+1][1],
                                          num=spaceNum) * 253, dtype = int)
                lut[iColor].extend(iColorArray)

            while len(lut[iColor]) < 254:
                lut[iColor].append(lut[iColor][-1])
            while len(lut[iColor]) > 254:
                del lut[iColor][-1]
        
        #ANTON: _create_list should not be a method. it is used only once and very small
        
        lut = map(self._create_list, lut[0], lut[1], lut[2])
        # 254 is black
        lut.append([0,0,0])
        # 255 is white
        lut.append([255,255,255])
        return self._flatten(lut)

    def _flatten(self, nestList):

        #ANTON: CAN THIS APPROACH REPLACE THIS METHOD:
        #       sum(nestList, [])  ?
        
        '''Create a flat list
        get rid of nests in the list

        Parameters
        ----------
        nestList :  a list

        Returns : abs list
        '''
        if isinstance(nestList, list):
            if nestList == []:
                return []
            else:
                return self._flatten(nestList[0]) + self._flatten(nestList[1:])
        else:
            return [nestList]

    def _get_digit(self, num):
        '''Return a number of digits
        get rid of nests in the list

        Parameters
        ----------
        num :  int / float

        Returns : int
        '''
        if num == 0:
            return 0
        else:
            return floor(log10(abs(num)))

    def _get_mapper(self, mapperName, bandList):
        '''Creare VRT file in memory (VSI-file) with variable mapping

        If mapperName is given, it is added as the first in the self.mapperList.
        Loop over all availble mappers in mapperList to get the matching one.
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

        # add the given mapper first
        self.mapperList = ['mapper_' + mapperName] + self.mapperList

        # try to import and get VRT datasaet from all mappers. Break on success
        # if none of the mappers worked - None is returned
        vrtDataset = None
        for iMapper in self.mapperList:
            try:
                #get rid of .py extension
                iMapper = iMapper.replace('.py', '')
                #import mapper
                mapper_module = __import__(iMapper)
                #create a Mapper object and get VRT dataset from it
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

    def _round_Scale(self, num, formatList):
        ''' return round value

        Parameters
        ----------
            num: int / float / exponential
            formatList : string "%d", "%4.2e" etc...

        Returns
        -------
            int / float / exponential
        '''
        return formatList %num

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
