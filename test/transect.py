#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 18:06:15 2012

@author: mag
"""

import matplotlib.pyplot as plt
from numpy import mgrid, sqrt, sin, asarray, round, \
                  sort, hypot, linspace, arange, zeros, vstack, float64, flipud
from sympy import Point, Segment, Symbol
from scipy.ndimage import map_coordinates
import gtk

from os import path
from scipy.io import savemat

import imcrop
import distancelib


import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 6, 15)
__modified__ = datetime.datetime(2012, 6, 15)
__version__  = "1.0"
__status__   = "Development"


def subs_point(l, val):
    """Take an arbitrary point and make it a fixed point."""
    t = Symbol('t', real=True)
    ap = l.arbitrary_point()
    #~ return Point(ap[0].subs(t, val), ap[1].subs(t, val))
    return Point(ap.subs(t, val))

def paralline(nl=2, x0=0, y0=100, x1=33, y1=66):
    p1, p2 = Point(x0, y0), Point(x1, y1)
    s1 = Segment(p1, p2)
    # Create a new Line perpendicular to s1 passing through the point p1
    l1 = s1.perpendicular_line(p1)
    l2 = s1.perpendicular_line(p2)
    p1 in l1
    p2 in l2
    s1.is_perpendicular(l2)
    s1.is_perpendicular(l1)
    #p11 = subs_point(l1, s1.length)
    # find coords of parallel nl segments from each side of the transect
    x11, y11 = zeros(2*nl+1), zeros(2*nl+1)
    x22, y22 = zeros(2*nl+1), zeros(2*nl+1)
    j=0
    for i in range(-nl,nl+1):
        p11 = subs_point(l1, 1*i/s1.length) # divide unit segment on its length
        x111, y111 = p11.args
        x11[j], y11[j] = float64(x111), float64(y111)
        p22 = subs_point(l2, 1*i/s1.length) # divide unit segment on its length
        x222, y222 = p22.args
        x22[j], y22[j] = float64(x222), float64(y222)
        j+=1
#    # Checking that segments are parallel and same length
#    s2 = Segment(p11,p22)
#    s2.is_parallel(s1)
#    s1.length - s2.length
#    plt.plot([x0, x1], [y0, y1])
#    plt.plot([x11, x22],[y11, y22], 'k')
    return x11, y11, x22, y22

def tellme(s):
    print s
    plt.title(s, fontsize=16)
    plt.draw()

def imagePts():
    """Function to extract transect from a figure"""
    # get the screen resolution
    window = gtk.Window()
    screen = window.get_screen()
    width = screen.get_width()
    height = screen.get_height()
    # resize the figure
    mng = plt.get_current_fig_manager()
    mng.resize(width, height)
    tellme( 'Now do a transect, click to begin' )
    plt.waitforbuttonpress()
    happy = False
    # get initial axes
    ax = plt.gca().get_xbound()
    ay = plt.gca().get_ybound()
    # initiate pt for returning the zoom coords
    pt = []
    # get default axes
    while not happy:
        tellme("Select two point of transect, Enter or middle mouse button to finish" )
        pts = asarray( plt.ginput(2, timeout=15, mouse_add=1, \
                        mouse_pop=3, mouse_stop=2) )
        happy = len(pts) < 2
        if happy: break
        try:
            pl0.pop(0).remove() # remove previous transect line and text
            pl1.remove()
            pl2.remove()
        except NameError:
            pl0, pl1, pl2 = [], [], []
        plt.draw()
        # new points
        x = pts[:, 0]
        y = pts[:, 1]
        # plot transect line
        plt.hold(True)
        pl0 = plt.plot(x, y, '-o')
        pl1 = plt.text(x[0],y[0],"A",bbox=dict(facecolor='c', alpha=0.3),stretch='expanded',fontsize=27)
        pl2 = plt.text(x[1],y[1],"B",bbox=dict(facecolor='c', alpha=0.3),stretch='expanded',fontsize=27)
        # make the original axes current again
        plt.gca().set_xbound(ax)
        plt.gca().set_ybound(ay)
        plt.draw()
        pts = sort(pts, axis=0)
        pt = [pt, pts]
    tellme('') # Finish
    plt.show()
    pts = pt[-1]  # return only last zoom coords
    pts[:,1] = flipud(pts[:,1]) # flip Y-coords after sort
    return pts

def dataPts(data=None):
    """Function to extract transect from inputted data"""
    if data is None:
        #-- Generate some data...
        x, y = mgrid[-5:5:0.01, -5:5:0.01]
        data = sqrt(x**2 + y**2) + sin(x**2 + y**2)
    # get the screen resolution
    window = gtk.Window()
    screen = window.get_screen()
    width = screen.get_width()
    height = screen.get_height()
    # plot the figure
    plt.figure(333)
    mng = plt.get_current_fig_manager()
    mng.resize(width, height)
    ims = plt.imshow(data)
    ims.axis = 'tight'
    plt.gray()
    plt.colorbar()
    # set the points for transect on image
    pts = imagePts()
    return pts

def transect(data=None, pts=None, nl=0, method=None):

    if data is None:
        #-- Generate some data...
        x, y = mgrid[-5:5:0.01, -5:5:0.01]
        data = sqrt(x**2 + y**2) + sin(x**2 + y**2)

    if pts is None:
        pts = dataPts(data)

    # Make a line with "num" points...
    x0, y0 = pts[0,:] # These are in _pixel_ coordinates!!!
    x1, y1 = pts[1,:]

    # Return the Euclidean norm, sqrt(x*x + y*y). This is the length
    # of the vector from the origin to point (x, y).
    length = int(hypot(x1-x0, y1-y0))
    x, y = linspace(x0, x1, length), linspace(y0, y1, length)

    # Parallel transect lines
    if nl is not 0:
        x11, y11, x22, y22 = paralline(nl, x0, y0, x1, y1)
        x = zeros((len(x11),length))
        y = zeros((x.shape))
        for i in range(0,len(x11)):
            x[i,:], y[i,:] = linspace(x11[i], x22[i], length), \
                             linspace(y11[i], y22[i], length)
        # x and y coords must be bounded by data size, rounding to 0
        x[x<0] = 0
        y[y<0] = 0
        x[x>=data.shape[1]] = data.shape[1]-1
        y[y>=data.shape[0]] = data.shape[0]-1
        # Extract the values along the line
        x = x.astype(int)
        y = y.astype(int)
#        x_ = diff(x, axis=0)
        trn = data[y, x]
        trnMean = trn.mean(axis=0)
        return trn, trnMean
    else:
        return data[y.astype(int), x.astype(int)]

    # use scipy.ndimage.map_coordinates for interpolation of "order="
    if method is not None:
        trn = zeros((x.shape))
        if method=='bilinear':
            for i in range(0,len(x11)-1):
                trn[i,:] = map_coordinates(data, vstack((y[i,:],x[i,:])), order=2)
        elif method=='cubic':
            for i in range(0,len(x11)-1):
                trn[i,:] = map_coordinates(data, vstack((y[i,:],x[i,:])), order=3)
        return trn, trn.mean(axis=0)



def plotTrnsImage(trn, pts, lat, lon, data, pixel, line,  pn='/home/mag/', \
             label='Wb contribution [linear units]', clm=(0,0.03), \
             ptsCrop=None, m=None, scale=1, fign=None):

    # Scaling and cropping the image and lat/lon
    lat_new, lon_new = imcrop.intCrpLL(lat, lon, pixel, line, scale, ptsCrop)
    if ptsCrop is not None:
        data = imcrop.imcrop(ptsCrop, data)
    
    data = data[::scale,::scale]
    pts = pts/scale
    pts = pts.astype(int)

    # 1. Plot the data image with transect
    
    if m is None:
        plt.figure()
        plt.imshow(data)
        plt.gray()
        plt.clim(clm)
        plt.colorbar()
        ax = plt.gca().get_xbound()
        ay = plt.gca().get_ybound()
        x = pts[:, 0]
        y = pts[:, 1]
        plt.plot(x, y, '-o')
        plt.text(x[0],y[0],"A",bbox=dict(facecolor='w', alpha=0.7),stretch='expanded',fontsize=27)
        plt.text(x[1],y[1],"B",bbox=dict(facecolor='w', alpha=0.7),stretch='expanded',fontsize=27)
        plt.gca().set_xbound(ax)
        plt.gca().set_ybound(ay)
        plt.draw()
        # resize the figure window
        mng = plt.get_current_fig_manager()
        mng.resize(1920,1080)
        plt.draw()
        a = label.split(' ')[0] # split the name if it has spaces
        plt.savefig(pn+a+'-with-trns.tiff', facecolor='w', edgecolor='w', \
                     dpi=300, bbox_inches="tight", pad_inches=1.75)
    elif type(m).__name__ is 'Basemap':
        print ( "Using existing Basemap \n" )
#    else:
#        print ( "Using default NSPER Basemap \n")
#        # Lat/Lon coords of image corners
#        ll_lat = lat_new.min()
#        ur_lat = lat_new.max()
#        ll_lon = lon_new.min()
#        ur_lon = lon_new.max()
#        cent_lat = lat_new.mean()
#        cent_lon = lon_new.mean()
#        m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
#                    llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
#                    resolution='i', projection='nsper', \
#                    satellite_height=798000, \
#                    lat_0=cent_lat,lon_0=cent_lon)

        color = 'black'
        
        # compute native map projection coordinates of lat/lon grid.
        x, y = m(lon_new,lat_new)
        
        mng = plt.get_current_fig_manager()
        mng.resize(1920,1080)
    
        cs = m.pcolormesh(x,y,data)
        cs.axis='tight'
        plt.gray()
        plt.clim(clm)
        cb = plt.colorbar(cs, shrink=0.8, extend='both', \
                          orientation='horizontal', pad=0.1, aspect=33)
        # A working example (for any value range) with 5 ticks along the bar is:
        m0=(clm[0])            # colorbar min value
        m5=(clm[1])             # colorbar max value
        m1=round((1*(m5-m0)/5.0 + m0),2)               # colorbar mid value 1
        m2=round((2*(m5-m0)/5.0 + m0),2)               # colorbar mid value 2
        m3=round((3*(m5-m0)/5.0 + m0),2)               # colorbar mid value 3
        m4=round((4*(m5-m0)/5.0 + m0),2)
        cb.set_ticks([m0,m1,m2,m3,m4,m5])
        cb.set_ticklabels([m0,m1,m2,m3,m4,m5])
        cb.update_ticks()
        cb.set_label(label)
        plt.draw()
    
        # set the step of Lat/Lon to plot
        stepLon = round((lon_new.max()-lon_new.min())/3, 1)
        stepLat = round((lat_new.max()-lat_new.min())/1.5, 1)
        # fool proofing, so that ronded value is not 0, but 0.05 in case of small area
        if stepLat == 0: stepLat=0.05
        if stepLon == 0: stepLon=0.05
        m.drawmeridians(arange(round(lon_new.min(),1),round(lon_new.max(),1), stepLon), \
                        labels=[0,0,0,1], color=color, dashes=[5,5], linewidth=0)
        m.drawparallels(arange(round(lat_new.min(),1),round(lat_new.max(),1), stepLat), \
                        labels=[1,0,0,0], color=color, dashes=[5,5], linewidth=0, rotation=90)
    
        # plot the transect
        lattrns = lat_new[pts[:,1], pts[:,0]]
        lontrns = lon_new[pts[:,1], pts[:,0]]
        # compute native map projection coordinates of lat/lon trns
        xtrns, ytrns = m(lontrns,lattrns)
        delta = xtrns[0]/30 # to plot text shifted several points from the marker
        # plot A-B label of the transect
        m.plot(xtrns, ytrns, 'w-o', linewidth=2)
        plt.text(xtrns[0]-delta,ytrns[0]+delta,"A",bbox=dict(facecolor='w', alpha=1),stretch='expanded',fontsize=27)
        plt.text(xtrns[1]+delta,ytrns[1]+delta,"B",bbox=dict(facecolor='w', alpha=1),stretch='expanded',fontsize=27)
        if fign is not None:
            xfign, yfign = m(lon_new.max(),lat_new.max())
            plt.text(xfign*0.1, yfign*0.9,fign,bbox=dict(facecolor='w', alpha=1),stretch='expanded',fontsize=27)
        plt.draw()
        a = label.split(' ')[0] # split the name if it has spaces
        plt.savefig(pn+a+'-with-trns.tiff', facecolor='w', edgecolor='w', \
                     dpi=300, bbox_inches="tight", pad_inches=0.1)


def plotTrns(trn, pts, lat, lon, data, pn='/home/mag/', \
             label='Wb contribution [linear units]', fign=None):
    """
    Plotting the transect
    """

    # Get the length of the transect. Input coords in (row,col)
    trnsLength = distancelib.getDistancePx([pts[0,1], pts[0,0]], [pts[1,1], pts[1,0]], lat, lon, data)
    dist = linspace(0, trnsLength, len(trn))

    # 2. Plot the transect
    plt.figure(figsize=(7, 5))
    plt.plot(dist, trn)
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.title(label)
    plt.axis('tight')
    plt.text(dist[15],(trn.max()+trn.min())/2,"A",bbox=dict(facecolor='w', alpha=0.7),stretch='expanded',fontsize=27)
    plt.text(dist[-55],(trn.max()+trn.min())/2,"B",bbox=dict(facecolor='w', alpha=0.7),stretch='expanded',fontsize=27)
    if fign is not None:
        plt.text(dist.max()*0.9,trn.max()-(trn.max()-trn.min())*0.1,fign,bbox=dict(facecolor='w', alpha=1),stretch='expanded',fontsize=27)
#    mng = plt.get_current_fig_manager()
#    mng.resize(1920,1080)
    plt.draw()
    a = label.split(' ')[0] # split the name if it has spaces
    plt.savefig(pn+a+'-trns.tiff', facecolor='w', edgecolor='w', \
                 dpi=300, bbox_inches="tight", pad_inches=0.1)

#    plt.close('all')
