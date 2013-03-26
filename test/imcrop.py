#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:57:55 2012

@author: mag
"""

import matplotlib.pyplot as plt

from numpy import linspace, asarray, sort
from scipy.interpolate import interp2d

def tellme(s):
    print s
    plt.title(s,fontsize=16)
    plt.draw()

def imzoom():
    """
    Crop the image and return the Bbox
    """
    tellme( 'Now do a nested zoom, click to begin' )
    plt.waitforbuttonpress()
    # use thismanager to enable resetting axes when pressing home button
    thismanager = plt.get_current_fig_manager()
    thismanager.toolbar.pan()
    happy = False
    # initiate pt for returning the zoom coords
    pt = []
    # get default axes
    while not happy:
        tellme("Select two corners of zoom, middle mouse button to finish" )
        pts = asarray( plt.ginput(2,timeout=-1, mouse_add=1, \
                        mouse_pop=3, mouse_stop=2) )

        happy = len(pts) < 2
        if happy: break

        pts = sort(pts,axis=0)
        plt.axis( pts.T.ravel() )

        pt = [pt, pts]

    tellme('Done!')
    plt.show()

    return pt[-1] # return only last zoom coords

def imcrop(pts, image):
    x, xx, y, yy = int(pts[0,0]), int(pts[1,0]), int(pts[0,1]), int(pts[1,1])
    imCropped = image[y:yy,x:xx]
    return imCropped
    
#def ptscrop(pts, image, ptsToCrop)
#    """return coords of initial points in a cropped image"""
#    x, xx, y, yy = int(pts[0,0]), int(pts[1,0]), int(pts[0,1]), int(pts[1,1])
#    imCropped = image[y:yy,x:xx]
#    ptsCropped = 

def interpLL(l, pixel, line, gx, gy):
    """
    Interpolating lat/lon to image size for future pcolormesh
    """
    fl = interp2d(pixel[0,:], line[:,0], l, kind='linear')
    l_new = fl(gx, gy)
    return l_new

def intCrpLL(lat, lon, pixel, line, scale=8, ptsCrop=None):

    if ptsCrop is not None:
        # Interpolating lat/lon to image size for future pcolormesh
        RasterXSize = round((-ptsCrop[0,0]+ptsCrop[-1,0])/scale)
        RasterYSize = round((-ptsCrop[0,-1]+ptsCrop[-1,-1])/scale)
        gx = linspace(ptsCrop[0,0], ptsCrop[-1,0], RasterXSize+1)
        gy = linspace(ptsCrop[0,-1], ptsCrop[-1,-1], RasterYSize+1)
    else:
        # Interpolating lat/lon to image size for future pcolormesh
        RasterXSize = round(pixel[0,-1]/scale)
        RasterYSize = round(line[-1,0]/scale)
        gx = linspace(0, pixel[0,-1], RasterXSize+1)
        gy = linspace(0, line[-1,0], RasterYSize+1)

    lat_new = interpLL(lat, pixel, line, gx, gy)
    lon_new = interpLL(lon, pixel, line, gx, gy)

    return lat_new, lon_new



##Another way do to imzoom
#def imzoom():
#    """
#    Zoom the image and return the Bbox
#    """
#    tellme( 'Now do a nested zoom, click to begin' )
#    plt.waitforbuttonpress()
#    plt.draw()
#    # get default axes
#    ax0 = plt.gca()
#    tellme("Use 'p' to Pan/Zoom" )
#    thismanager = plt.get_current_fig_manager()
#    thismanager.toolbar.zoom()
#    plt.draw()
#    ax = plt.gca()
#    plt.show()
#    # for some reason doesn't return pts, only ax object
#    return ax # return zoom coords
#    # that is why should to find pts after the function call
#    pts = plt.gca().viewLim.corners()