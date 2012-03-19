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
    def __init__(self, array):
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

        self.extensionList = ["png", "PNG", "tif", "TIF", "bmp",
                              "BMP", "jpg", "JPG", "jpeg", "JPEG"]

    def apply_logarithm(self, gamma=2.0):
        '''Modify values with logarithm'''
        self.gamma = gamma
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = (np.power((self.array[iBand, :, :] -
                                  self.cmin[iBand]) /
                                  (self.cmax[iBand] - self.cmin[iBand]),
                                  (1.0 / self.gamma))) * ((self.cmax[iBand] - self.cmin[iBand])) + self.cmin[iBand]

    def clim_from_histogram(self, ratio):
        ratioList = np.ones(3)
        if isinstance(ratio, float) and self.array.shape[0] == 3:
            ratioList *= ratio
        elif isinstance(ratio, list):
            for iRatio in range(min(3, len(ratio))):
                ratioList[iRatio] = ratio[iRatio]
        else:
            ratioList[0] = ratio

        clim = np.zeros([self.array.shape[0],2])
        for iBand in range(self.array.shape[0]):
            if ratioList[iBand] != 1.0:
                hist, bins = self._get_histogram(iBand)
                cumhist = hist.cumsum()
                cumhist /= cumhist[-1]
                if ratioList[iBand] < 0.0 or ratioList[iBand] > 1.0:
                    ratioList[iBand] = 1.0
                clim[iBand] = (bins[len(cumhist[cumhist<(1-ratioList[iBand])/2])], bins[len(cumhist[cumhist<1-((1-ratioList[iBand])/2)])])
            else:
                clim[iBand] = (self.array[iBand, :, :].min(), self.array[iBand, :, :].max())
        return clim

    def clip(self):
        '''Convert array to uint8 with scaling from cmin to cmax'''
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = np.clip(self.array[iBand, :, :], self.cmin[iBand], self.cmax[iBand])

    def convert_palettesize(self, numOfColor=None):
        if numOfColor is not None:
            self.numOfColor = numOfColor
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = ((self.array[iBand, :, :] - self.cmin[iBand]) * self.numOfColor / (self.cmax[iBand] - self.cmin[iBand]))
        self.array = self.array.astype(np.uint8)

    def create_legend(self, numOfTicks=5, barFontSize = 10, longName=None, units=None, titleString="", fontSize=10):
        '''Extend self.array and Create PIL image with legend'''
        # create a pilImage for the legend
        self.pilImgLegend = Image.new("P", (self.sizeX,
                                        int(self.sizeY * 0.1)), 255)
        # make an array for color bar
        bar = np.outer(np.ones(max(int(self.pilImgLegend.size[1]*0.15),
                       5)), np.linspace(0, self.numOfColor,
                       int(self.pilImgLegend.size[0]*0.8)))
        # create a colorbar pil Image
        pilImgCbar = Image.fromarray(np.uint8(bar))
        # paste the colorbar pilImage on Legend pilImage
        self.pilImgLegend.paste(pilImgCbar,
                           (int(self.pilImgLegend.size[0]*0.1),
                            int(self.pilImgLegend.size[1]*0.5)))
        # create a scale for the colorbar
        scaleLocation = np.linspace(0, 1, numOfTicks)
        scaleArray = scaleLocation
        if self.gamma is not None:
            scaleArray = (np.power(scaleArray, (1.0/self.gamma)))
        scaleArray = scaleArray * (self.cmax[0] - self.cmin[0]) + self.cmin[0]
        scaleArray = map(self._round_number, scaleArray)
        # set fonts for colorbar
        fileName_font = os.path.join(os.path.dirname(
                                     os.path.realpath(__file__)),
                                     'fonts/times.ttf')
        font = ImageFont.truetype(fileName_font, barFontSize)
        # draw scales and lines on the legend pilImage
        draw = ImageDraw.Draw(self.pilImgLegend)
        for iTick in range(numOfTicks):
            coordX = int(scaleLocation[iTick]*
                         self.pilImgLegend.size[0]*0.8 +
                         self.pilImgLegend.size[0]*0.1)
            box = (coordX, int(self.pilImgLegend.size[1]*0.5),
                   coordX, int(self.pilImgLegend.size[1]*0.65)-1)
            draw.line(box, fill= 254)
            box = (coordX-5, int(self.pilImgLegend.size[1]*0.65)+3)
            draw.text(box, scaleArray[iTick], fill= 254, font=font)

        # set font size for text
        font = ImageFont.truetype(fileName_font, fontSize)
        if longName is None:
            longName = "no longName"
        if units is None:
            units = "no units"
        caption = longName + " / " + units
        box = (int(self.pilImgLegend.size[0]*0.1),
               int(self.pilImgLegend.size[1]*0.3))
        draw.text(box, str(caption), fill= 254, font=font)

        if titleString != "":
            # write text each line onto pilImgCanvas
            textHeight = int(self.pilImgLegend.size[1]*0.1)
            for line in titleString.splitlines():
                draw.text((int(self.pilImgLegend.size[0]*0.1), textHeight),
                           line, fill= 254, font=font)
                text = draw.textsize(line, font=font)
                textHeight += text[1]

    def create_pilImage(self, cmapName="jet"):
        if self.array.shape[0] == 3:
            self.pilImg = Image.merge("RGB",
                                     (Image.fromarray(self.array[0, :, :]),
                                      Image.fromarray(self.array[1, :, :]),
                                      Image.fromarray(self.array[2, :, :])))
        else:
            myPalette = self._create_palette(cmapName)
            if self.pilImgLegend is not None:
                extentArray = np.append(self.array[0, :, :], np.ones((self.pilImgLegend.size[1], self.pilImgLegend.size[0]), np.uint8), 0)
                self.array = extentArray.reshape(1, extentArray.shape[0], extentArray.shape[1])
                extentArray = None
            self.pilImg = Image.fromarray(self.array[0, :, :])
            self.pilImg.putpalette(myPalette)
            if self.pilImgLegend is not None:
                # paste the legend on the pilImg
                self.pilImg.paste(self.pilImgLegend, (0, self.sizeY))
                self.sizeX = self.array.shape[2]
                self.sizeY = self.array.shape[1]
        self.array = None

    def save(self, fileName):
        if fileName.split(".")[-1] in self.extensionList:
            self.pilImg.save(fileName)
        else:
            self.pilImg.save(fileName + ".png")
        self.pilImg = None

    def set_clim(self, clim):
        '''set self values of clim'''
        self.cmin = clim[:, 0]
        self.cmax = clim[:, 1]

    def _create_palette(self, cmapName):
        try:
            colorDic = cm.datad[cmapName]
        except:
            colorDic = cm.datad["jet"]

        colorList = [colorDic['red'], colorDic['green'], colorDic['blue']]

        lut = np.zeros([3,256])
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
            while iPalette < self.numOfColor:
                lut[iColor][iPalette] = lut[iColor][iPalette-1]
                iPalette += 1
            while iPalette > self.numOfColor:
                lut[iColor][iPalette] = 0
                iPalette -= 1
        # the last palette is white
        for iColor in range(3):
            lut[iColor][-1] = 255
        return lut.T.flatten().astype(np.uint8)

    def _get_histogram(self, iBand):
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
            values: int / float / exponential

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


