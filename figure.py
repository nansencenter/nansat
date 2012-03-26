#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      asumak
#
# Created:     13.03.2012
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os
import numpy as np
import Image
import ImageDraw
import ImageFont
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import outer
from math import floor
from math import log10

class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass

class OptionError(Error):
    '''Error for improper options (arguments) '''
    pass

class Figure():
    '''Perform opeartions with graphical files: create, append legend, modify...
    
    Figure instance is created in the Nansat.write_figure method
    The methods below are applied consequently in order to generate a figure
    from one or three bands, estimate min/max, apply logarithmic scaling,
    convert to uint8, append legend, save to a file    
    '''
    
    def __init__(self, array):
        ''' Set attributes

        Parameters
        ----------
        array: numpy array (2D or 3D)
            dataset from Nansat

        Modifies
        --------
        self.sizeX, self.sizeY : int
            X and Y sizes of the image
        self.cmin, self.cmax : float
            minimum and maximum pixel values of the image
        self.gamma : float
            if gamma is None, tone curve is not applied
            if it is replaced, gamma means a coefficient of tone curve
        self.pilImg : PIL image
            figure
        self.pilImgLegend : PIL image
            if pilImgLegend is None, legend is not added to the figure
            if it is replaced, pilImgLegend includes text string, color-bar,
            longName and units.
        self.subsetArraySize : int
            size of the subset array which is used to get the histogram.
        self.numOfColor : int (max 256)
            number of colors for use of the palette.
            default is 250 colors. 254th is black and 255th is white.
        self.cmap : name of Matplotlib colormaps
            see --> http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        '''
        # if 2D array is given, reshape to 3D
        if array.ndim == 2:
            self.array = array.reshape(1, array.shape[0], array.shape[1])
        else:
            self.array = array

        self.sizeX = self.array.shape[2]
        self.sizeY = self.array.shape[1]

        self.cmin = None
        self.cmax = None
        self.gamma = None
        self.pilImg = None
        self.pilImgLegend = None

        self.subsetArraySize = 100000
        self.numOfColor = 250
        self.cmap = "jet"

        self.extensionList = ["png", "PNG", "tif", "TIF", "bmp",
                              "BMP", "jpg", "JPG", "jpeg", "JPEG"]

    def apply_logarithm(self, gamma=2.0):
        '''Apply a tone curve to the array
        After the normalization of the values from 0 to 1, logarithm is applied.
        Then the values are converted to the normal scale.

        Parameters
        ----------
        self.gamma: positive float, default is 2.0, optional
            a coefficient of tone curve. (the coefficient is 1/gamma)
            if 0 < gamma < 1.0, tone curve is convex downward.
            if gamma = 1.0, it is linear.
            if 1.0 < gamma, tone curve is convex upward.

        Modifies
        --------
        self.array : numpy array

        '''
        self.gamma = gamma
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = (np.power((self.array[iBand, :, :] -
                                  self.cmin[iBand]) /
                                  (self.cmax[iBand] - self.cmin[iBand]),
                                  (1.0 / self.gamma))) * ((self.cmax[iBand] -
                                  self.cmin[iBand])) + self.cmin[iBand]

    def clim_from_histogram(self, ratio=1.0):
        '''Get min and max pixel values.

        if ratio=1.0, simply the minimum and maximum values are returned.
        if 0 < ratio < 1.0, get the histogram of the pixel values.
        Then get rid of (1.0-ratio)/2 from the both sides and
        return the minimum and maximum values.

        Parameters
        ----------
        ratio: float, [0.0, 1.0]
            ratio of pixels which are used to write the figure.

        Return
        --------
        clim : numpy array 2D ((3x2) or (1x2))
            minimum and maximum pixel values for each band

        '''
        # create a ratio list for each band
        if isinstance(ratio, float):
            ratioList = np.ones(self.array.shape[0])*ratio
        else:
            ratioList = []
            for iRatio in range(self.array.shape[0]):
                try:
                    ratioList.append(ratio[iRatio])
                except:
                    ratioList.append(ratio[0])

        # create a 2D array and set min and max values
        clim = np.zeros([self.array.shape[0],2])
        for iBand in range(self.array.shape[0]):
            if ratioList[iBand] == 1.0:
                clim[iBand] = (self.array[iBand, :, :].min(),
                               self.array[iBand, :, :].max())
            else:
                hist, bins = self._get_histogram(iBand)
                cumhist = hist.cumsum()
                cumhist /= cumhist[-1]
                if ratioList[iBand] < 0.0 or ratioList[iBand] > 1.0:
                    ratioList[iBand] = 1.0
                clim[iBand] = (bins[len(cumhist[cumhist<
                               (1-ratioList[iBand])/2])],
                               bins[len(cumhist[cumhist<
                               1-((1-ratioList[iBand])/2)])])
        return clim

    def clip(self):
        '''Convert self.array to values between cmin and cmax

        if pixel value < cmin, replaced to cmin.
        if pixel value > cmax, replaced to cmax.

        Modifies
        --------
            self.array : numpy array

        '''
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = np.clip(self.array[iBand, :, :],
                                              self.cmin[iBand],
                                              self.cmax[iBand])

    def convert_palettesize(self, numOfColor=None):
        '''Convert self.array to palette color size in uint8

        Parameters
        ----------
        numOfColor: int, optional
            number of colors which are possible to use for the palette

        Modifies
        --------
            self.array : numpy array (=>uint8)

        '''
        if numOfColor is not None:
            self.numOfColor = numOfColor
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = ((self.array[iBand, :, :] -
                                       self.cmin[iBand]) * (self.numOfColor-1) /
                                       (self.cmax[iBand] - self.cmin[iBand]))
        self.array = self.array.astype(np.uint8)

    def create_legend(self, numOfTicks=5, longName=None, units=None,
                      titleString="", fontSize=10):
        ''' self.legend is replaced from None to PIL image

        PIL image includes colorbar, longname, units and titelString.

        Parameters
        ----------
        numOfTicks : int, optional
            number of ticks for the colorbar
        longName, units : string, optional
            given from nansat WKV
        titleString : string, optional
        fontSize : int
            font size for the longName, units and titleString

        Modifies
        --------
        self.legend : PIL image

        '''
        LEGEND_HEIGHT = 0.1 # percentage of the image height
        CBAR_HEIGHTMIN = 5
        CBAR_HEIGHT = 0.15 # percentage of the image height
        CBAR_WIDTH = 0.8 # percentage of the image width
        CBAR_LOCATION_X = 0.1 # percentage of the legend (=image) width
        CBAR_LOCATION_Y = 0.5 # percentage of the legend height
        CBAR_LOCATION_ADJUST_X = 5 # pixel
        CBAR_LOCATION_ADJUST_Y = 3 # pixel
        TEXT_LOCATION_X = 0.1 # percentage of the legend (=image) width
        TEXT_LOCATION_Y = 0.1 # percentage of the legend height
        NAME_LOCATION_X = 0.1 # percentage of the legend (=image) width
        NAME_LOCATION_Y = 0.3 # percentage of the legend height

        # create a pilImage for the legend
        self.pilImgLegend = Image.new("P", (self.sizeX,
                                      int(self.sizeY * LEGEND_HEIGHT)), 255)
        draw = ImageDraw.Draw(self.pilImgLegend)

        # set fonts for Legend
        fileName_font = os.path.join(os.path.dirname(
                                     os.path.realpath(__file__)),
                                     'fonts/times.ttf')

        # set black color
        if self.array.shape[0] == 1:
            black = 254
        else:
            black = (0,0,0)

        # if 1 band, draw the color bar
        if self.array.shape[0] == 1:
            # make an array for color bar
            bar = np.outer(np.ones(max(int(self.pilImgLegend.size[1]*CBAR_HEIGHT),
                           CBAR_HEIGHTMIN)), np.linspace(0, self.numOfColor,
                           int(self.pilImgLegend.size[0]*CBAR_WIDTH)))
            # create a colorbar pil Image
            pilImgCbar = Image.fromarray(np.uint8(bar))
            # paste the colorbar pilImage on Legend pilImage
            self.pilImgLegend.paste(pilImgCbar,
                               (int(self.pilImgLegend.size[0]*CBAR_LOCATION_X),
                                int(self.pilImgLegend.size[1]*CBAR_LOCATION_Y)))
            # create a scale for the colorbar
            scaleLocation = np.linspace(0, 1, numOfTicks)
            scaleArray = scaleLocation
            if self.gamma is not None:
                scaleArray = (np.power(scaleArray, (1.0/self.gamma)))
            scaleArray = (scaleArray * (self.cmax[0] - self.cmin[0]) + 
                        self.cmin[0])
            scaleArray = map(self._round_number, scaleArray)
            # set fonts size for colorbar
            font = ImageFont.truetype(fileName_font, fontSize)
            # draw scales and lines on the legend pilImage
            for iTick in range(numOfTicks):
                coordX = int(scaleLocation[iTick]*
                             self.pilImgLegend.size[0]*CBAR_WIDTH +
                             self.pilImgLegend.size[0]*((1-CBAR_WIDTH)/2))
                box = (coordX, int(self.pilImgLegend.size[1]*CBAR_LOCATION_Y),
                       coordX, int(self.pilImgLegend.size[1]*(CBAR_LOCATION_Y +
                       CBAR_HEIGHT))-1)
                draw.line(box, fill= black)
                box = (coordX-CBAR_LOCATION_ADJUST_X,
                       int(self.pilImgLegend.size[1]*(CBAR_LOCATION_Y +
                           CBAR_HEIGHT)) + CBAR_LOCATION_ADJUST_Y)
                draw.text(box, scaleArray[iTick], fill= black, font=font)

        # set font size for text
        font = ImageFont.truetype(fileName_font, fontSize)

        # draw longname and units
        if longName is None:
            longName = ""
        if units is None:
            units = ""
        caption = longName + " / " + units
        box = (int(self.pilImgLegend.size[0]*NAME_LOCATION_X),
               int(self.pilImgLegend.size[1]*NAME_LOCATION_Y))
        draw.text(box, str(caption), fill= black, font=font)

        # if titleString is given, draw it
        if titleString != "":
            # write text each line onto pilImgCanvas
            textHeight = int(self.pilImgLegend.size[1]*TEXT_LOCATION_Y)
            for line in titleString.splitlines():
                draw.text((int(self.pilImgLegend.size[0]*TEXT_LOCATION_X),
                           textHeight), line, fill= black, font=font)
                text = draw.textsize(line, font=font)
                textHeight += text[1]

    def create_pilImage(self, cmapName=None):
        ''' self.create_pilImage is replaced from None to PIL image

        If three images are given, create a image with RGB mode.
            if self.pilImgLegend is not None, it is pasted.
        If one image is given, create a image with P(palette) mode.
            if self.pilImgLegend is not None,
            self.array is extended before create the pilImag and
            then paste pilImgLegend onto it.

        Parameters
        ----------
        cmapName : name of Matplotlib colormaps, optional
            see --> http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        Modifies
        --------
        self.pilImg : PIL image
            PIL image with / without the legend
        self.array : replace to None

        '''
        # if three bands are given, create a new PIL image with RGB mode
        if self.array.shape[0] == 3:
            # if self.pilImgLegend is None,
            # create a PILimage whose size is same as image.
            if self.pilImgLegend is None:
                self.pilImg = Image.new("RGB", (self.sizeX, self.sizeY))
            # if self.pilImgLegend is given,
            # create a PILimage whose size is (image + legend).
            else:
                self.pilImg = Image.new("RGB", (self.sizeX,
                                        self.sizeY+self.pilImgLegend.size[1]))
                self.pilImg.paste(self.pilImgLegend, (0, self.sizeY))
            # paste RGB image onto the empty self.pilImgLegend
            self.pilImg.paste(Image.merge("RGB",
                             (Image.fromarray(self.array[0, :, :]),
                              Image.fromarray(self.array[1, :, :]),
                              Image.fromarray(self.array[2, :, :]))), (0, 0))

        else:
            # create a palette
            myPalette = self._create_palette(cmapName)
            if self.pilImgLegend is not None:
                # if self.pilImgLegend is given, self.array is extended
                extentArray = np.append(self.array[0, :, :],
                              np.ones((self.pilImgLegend.size[1],
                              self.pilImgLegend.size[0]), np.uint8), 0)
                self.array = extentArray.reshape(1, extentArray.shape[0],
                                                 extentArray.shape[1])
                extentArray = None
            # create a pilImg from self.array and put color with the palette
            self.pilImg = Image.fromarray(self.array[0, :, :])
            self.pilImg.putpalette(myPalette)
            if self.pilImgLegend is not None:
                # paste the legend on the pilImg
                self.pilImg.paste(self.pilImgLegend, (0, self.sizeY))
                # replace self.sizeX and self.sizeY
                self.sizeX = self.array.shape[2]
                self.sizeY = self.array.shape[1]
        self.array = None

    def save(self, fileName):
        ''' Save self.pilImg to a physical file

        If the given extension is JPG, convert the image mode from Palette to RGB

        Parameters
        ----------
        fileName : string
            name of outputfile

        Modifies
        --------
        self.pilImg : None

        '''
        # set defaults
        DEFAULT_EXTENSION = ".png"
        
        if (fileName.split(".")[-1]=="jpg") or \
           (fileName.split(".")[-1]=="JPG") or \
           (fileName.split(".")[-1]=="jpeg") or \
           (fileName.split(".")[-1]=="JPEG"):
            self.pilImg = self.pilImg.convert("RGB")
            
        if not((fileName.split(".")[-1] in self.extensionList)):
            fileName = fileName + DEFAULT_EXTENSION
        
        self.pilImg.save(fileName)        

    def set_clim(self, clim):
        '''set self.cmin and self.max

        Modifies
        --------
        self.cmin : a list
            minimum values for each image
        self.cmax : a list
            maximum values for each image

        '''
        self.cmin = clim[:, 0]
        self.cmax = clim[:, 1]

    def _create_palette(self, cmapName):
        '''Create a palette based on Matplotlib colormap name

        default number of color palette is 250.
        it means 6 colors are possible to use for other purposes.
        the last palette (255) is white and the second last (254) is black.

        Returns
        -------
        color palette : numpy array (uint8)

        '''
        # create a colorList from matplotlib colormap name
        try:
            colorDic = cm.datad[cmapName]
        except:
            colorDic = cm.datad[self.cmap]

        colorList = [colorDic['red'], colorDic['green'], colorDic['blue']]

        # create a numpyarray for a palette based on the color list
        # default is all values are black (=0)
        lut = np.zeros([3,256])
        # place colors to each number (palette)
        for iColor in range(3):
            iPalette = 0
            for i in range(len(colorList[iColor])-1):
                spaceNum = int(self.numOfColor * (colorList[iColor][i+1][0] -
                               colorList[iColor][i][0]))
                lut[iColor][iPalette:iPalette+spaceNum] = np.array(np.linspace(
                                                  colorList[iColor][i][2],
                                                  colorList[iColor][i+1][1],
                                                  num=spaceNum) *
                                                  (self.numOfColor-1),
                                                  dtype = np.uint8)
                iPalette += (spaceNum)
            # adjust the number of colors on the palette
            while iPalette < self.numOfColor:
                lut[iColor][iPalette] = lut[iColor][iPalette-1]
                iPalette += 1
            while iPalette > self.numOfColor:
                lut[iColor][iPalette] = 0
                iPalette -= 1
        # the last palette color is replaced to white
        for iColor in range(3):
            lut[iColor][-1] = 255
        # return the palette
        return lut.T.flatten().astype(np.uint8)

    def _get_histogram(self, iBand):
        '''Create a subset array and return the histogram.

        Parameters
        ----------
        iBand: int

        Returns
        --------
        hist : numpy array
        bins : numpy array

        '''
        array = self.array[iBand, :, :].flatten()
        array = array[array>array.min()]
        array = array[array<array.max()]
        step = max(int(round(float(len(array)) / float(self.subsetArraySize))), 1.0)
        arraySubset = array[::step]
        hist, bins, patches = plt.hist(arraySubset, bins=100)
        return hist.astype(float), bins

    def _round_number(self, val):
        '''Return writing format for scale on the colorbar

        Parameters
        ----------
        val: int / float / exponential

        Returns
        --------
        string

        '''
        if val==0:
            return str("%d" % (val))
        else:
            digit = floor(log10(abs(val)))
            if digit==0:
                return str("%4.2f" % (val))
            elif digit==1:
                return str("%4.1f" % (val))
            elif digit==2:
                return str("%d" % (val))
            elif digit==-1:
                return str("%3.1f" % (val))
            elif digit==-2:
                return str("%4.2f" % (val))
            else:
                return str("%4.2e" % (val))


