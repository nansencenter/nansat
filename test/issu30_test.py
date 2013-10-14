#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      asumak
#
# Created:     14.10.2013
# Copyright:   (c) asumak 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os
from nansat import *
from nansat_tools import *
from mpl_toolkits.basemap import Basemap
import numpy as np
from figure import Figure
from matplotlib.pyplot import imshow, show
from numpy import array
from scipy import rand
from vrt import VRT
#import gdalinfo
import gdal, ogr
import time
import  tarfile
import itertools
import math


def shift_vrt(vrtObj, dataset ,shiftDegree):
    ''' Shift bands and modify geoTransform and return a shifted VRT

    Parameters
    ----------
    xBorderDegree : float
        values of degree to shift.
    xBorderPix : int
        number of pixels to shift.

    Modifies
    -------
    self.vrt
    '''
    # create a new VRT object
    shiftVRT = vrtObj.copy()

    # Copy VRT into self.vrt
    shiftVRT.vrt = vrtObj.copy()

    # if geoTransform is given, modify the 1st value (= top left x )
    if shiftDegree < 0:
        shiftDegree += 360.0

    geoTransform= dataset.GetGeoTransform()

    if geoTransform != '':
        shiftPixel = int(shiftDegree / geoTransform[1])
        geoTransform = list(geoTransform)
        geoTransform[0] = round(geoTransform[0] + shiftDegree, 3)
        if geoTransform[0] + geoTransform[1] * dataset.RasterXSize > 360.0:
            geoTransform[0] -= 360.0
        shiftVRT.dataset.SetGeoTransform(tuple(geoTransform))

    # write dataset content into VRT-file
    dataset.FlushCache()

    # Add bands to self
    for iBand in range(shiftVRT.vrt.dataset.RasterCount):
        src = {'SourceFilename': shiftVRT.vrt.fileName, 'SourceBand': iBand + 1}
        dst = shiftVRT.vrt.dataset.GetRasterBand(iBand+1).GetMetadata()
        shiftVRT._create_band(src, dst)

    # read xml and create the node
    XML = shiftVRT.read_xml()
    node0 = Node.create(XML)

    # divide into two bands and switch the bands
    for i in range(len(node0.nodeList('VRTRasterBand'))):
        # create i-th 'VRTRasterBand' node
        node1 = node0.node('VRTRasterBand', i)
        # modify the 1st band
        node1.node('ComplexSource').node('DstRect').replaceAttribute('xOff', str(shiftPixel))
        # add the 2nd band
        xmlSource = node1.xml()
        dom = xdm.parseString(xmlSource)
        cloneNode = Node.create(dom).node('ComplexSource')
        cloneNode.node('SrcRect').replaceAttribute('xOff', str(shiftVRT.dataset.RasterXSize - shiftPixel))
        cloneNode.node('DstRect').replaceAttribute('xOff', str(0))
        contents = node0.insert(cloneNode.xml(), 'VRTRasterBand', i)
        # overwrite the modified contents and create a new node
        dom = xdm.parseString(contents)
        node0 = Node.create(dom)

    #write XML contents to
    shiftVRT.write_xml(str(node0.rawxml()))

    return shiftVRT


n1 = Nansat('E:/data/NCEP_GRIB/gfs.t06z.master.grbf03')
dataset = n1.vrt.dataset
shiftedVRT = shift_vrt(n1.vrt, dataset, -160.0)
n1.vrt=shiftedVRT
n1.write_figure(True)













