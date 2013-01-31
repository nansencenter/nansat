#-------------------------------------------------------------------------------
# Name:        sample_nansatBasemap
# Purpose:
#
# Author:      asumak
#
# Created:     28.01.2013
# Copyright:   (c) asumak 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import scipy.ndimage as ndimage
from nansat import *
import math

def Map(lon, lat, dpi=300, projection='cyl', resolution='l',
        continent=True, geoCoordinates=True, mapboundary = False):

    m = Basemap(projection=projection, llcrnrlon=lon.min(),
                llcrnrlat=lat.min(), urcrnrlon=lon.max(),
                urcrnrlat=lat.max(), resolution=resolution)
    m.lat = lat
    m.lon = lon
    m.x, m.y = m(lon, lat)
    m.dpi = dpi
    m.cs = None
    m.continent = continent
    m.geoCoordinates = geoCoordinates
    m.colorbar = False
    m.mapboundary = mapboundary

    # create figure
    plt.close()
    figureSize = 3,3
    m.fig = plt.figure(num=1, figsize=figureSize, dpi=dpi)
    return m

def contour_plots(m, data, style="line", level=8, linewidths=0.2, cmap=plt.cm.jet,
                  clabel=False, inline=1, fontsize=3, smoothing=False, mode="gaussian"):
    if smoothing:
        data = image_process(data, mode)
    if style == "fill":
        m.cs = m.contourf(m.x, m.y, data, level, cmap=cmap)
        m.colorbar=True
    else:
        # draw contour lines
        m.cs = m.contour(m.x, m.y, data, level,
                       linewidths=linewidths, cmap=cmap)
        # add values for the contour lines
        if clabel:
            plt.clabel(m.cs, inline=inline, fontsize=fontsize)
    if m.continent:
            draw_continent(m)
    if m.geoCoordinates:
            draw_geoCoordinates(m)

def quiver_plots(m, dataX, dataY, quivectors=30):
    # subsample for quiver plot
    step0, step1 = dataX.shape[0]/quivectors, dataX.shape[1]/quivectors
    dataX2 = dataX[::step0, ::step1]
    dataY2 = dataY[::step0, ::step1]
    lon2 = m.lon[::step0, ::step1]
    lat2 = m.lat[::step0, ::step1]

    x, y = m(lon2, lat2)
    im2 = m.quiver(x, y, dataX2, dataY2)
    if m.continent:
            draw_continent(m)
    if m.geoCoordinates:
            draw_geoCoordinates(m)

def put_color(m, data, shading='flat',cmap=plt.cm.jet):
    m.pcolormesh(m.x, m.y, data, shading=shading, cmap=cmap)

def image_process(data, mode="gaussian", sigma=2.5, order=0, weight=None,
                  weightMtxSize=7, convMode="constant", cval=0.0,
                  splineOrder=1.0):
    if mode=="convolve":
        # if weight is None, create a weight matrix
        if weight is None:
            weight = np.ones((weightMtxSize, weightMtxSize))
            center = (weightMtxSize - 1) / 2
            for i in range(-(center), center+1, 1):
                for j in range(-(center), center+1, 1):
                    weight[i][j] /= math.pow(2.0, max(abs(i),abs(j)))
        return ndimage.convolve(data, weight, mode=convMode, cval=cval)
    elif mode=="fourier_gaussian":
        return ndimage.fourier_gaussian(data, sigma=sigma)
    elif mode=="spline":
        return ndimage.spline_filter1d(data, order=splineOrder)
    else:
        if mode!="gaussian":
            print "apply Gaussian filter in image_process()"
        return ndimage.gaussian_filter(data, sigma=sigma, order=order)

def add_legend(m, orientation='horisontal', pad=0.01,
               tickFontSize=4,
               title="", titleFontSize=5):
    # add colorbar and reduce font size
    if m.colorbar:
        cbar = m.fig.colorbar(m.cs, orientation=orientation, pad=pad)
        imaxes = plt.gca()
        plt.axes(cbar.ax)
        plt.xticks(fontsize=tickFontSize)
        plt.axes(imaxes)
    # add title
    if title != "":
        plt.title(title, fontsize=titleFontSize)

def draw_continent(m, continentColor='#cc9966', lakeColor='#99ffff'):
    m.fillcontinents(color=continentColor, lake_color=lakeColor)

def draw_geoCoordinates(m,
                        latNum=5, latFontsize=4,
                        latLabels=[True, False, False, False],
                        lonNum=5, lonFontsize=4,
                        lonLabels=[False, False, True, False]):
    # draw lat and lon
    m.drawparallels(np.arange(m.lat.min(), m.lat.max(),
                              (m.lat.max()-m.lat.min())/latNum),
                    labels=latLabels, fontsize=latFontsize)
    m.drawmeridians(np.arange(m.lon.min(), m.lon.max(),
                              (m.lon.max()-m.lon.min())/lonNum),
                    labels=lonLabels, fontsize=lonFontsize)

def draw_mapboudary(m, lineWidth=1, color="k", fillColor='0.3'):
    m.drawmapboundary(linewidth=lineWidth, color=color, fill_color=fillColor)

def save_map(m, fileName):
    m.fig.savefig(fileName, dpi=m.dpi)

#------------------------------------------------------------------------------#

# file with wind data
iFileName = 'c:/Users/asumak/Data/input/NCEP_GRIB/gfs.t06z.master.grbf03'
n = Nansat(iFileName)
##print n

# Norwegain and Barents Seas
d = Domain('+proj=longlat', '-te 0 60 30 80 -ts 300 300')

n.reproject(d)
lon,lat = n.get_geolocation_grids()

u = n[1]
v = n[2]
w = n[3]

# Create map
nMap = Map(lon, lat)

# add image
put_color(nMap, w)

# add contour1 (line)
contour_plots(nMap, v, level=8, clabel=True, smoothing=True)

# add contour2 (fill)
##contour_plots(nMap, w, style="fill", cmap=plt.cm.winter, smoothing=True, mode ='convolve')

# add quiver
##quiver_plots(nMap, u, v)

# add colorbar and title
add_legend(nMap, title='NCEP wind speed and direction')

# save to file
save_map(nMap, 'c:/Users/asumak/Data/output/basemap02.png')

