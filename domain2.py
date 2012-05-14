# Name:     domain2.py
# Purpose:  Optional methods for Domain. Depend on Basemap, Polygon, Matplotlib
#
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad
#
# Created:     01.05.2012
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

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np


def write_map(self, outputFileName, lonBorder=10., latBorder=10.,
                                    figureSize=(6, 6), dpi=50,
                                    projection='cyl', resolution='c',
                                    continetsColor='coral',
                                    meridians=10, parallels=10,
                                    pColor='r', pLine='k', pAlpha=0.5,
                                    padding=0.):
    ''' Create an image with a map of the domain

    Uses Basemap to create a World Map
    Adds a semitransparent patch with outline of the Domain
    Writes to an image file

    Parameters
    ----------
        outputFileName : string
            name of the output file name
        lonBorder : float
            10, horisontal border around patch (degrees of longitude)
        latBorder : float
            10, vertical border around patch (degrees of latitude)
        figureSize : tuple of two integers
            (6, 6), size of the generated figure in inches
        dpi : integer
            50, resolution of the output figure (size 6,6 and dpi 50
            produces 300 x 300 figure)
        projection : string, one of Basemap projections
            'cyl', projection of the map
        resolution : string, 'c', 'h', ...
            'c', crude, resolution of the map.
        continetsColor : string or any matplotlib color representation
            'coral', color of continets
        meridians : int
            10, number of meridians to draw
        parallels : int
            10, number of parallels to draw
        pColor : string or any matplotlib color representation
            'r', color of the Domain patch
        pLine :  string or any matplotlib color representation
            'k', color of the Domain outline
        pAlpha : float 0 - 1
            0.5, transparency of Domain patch
        padding : float
            0., width of white padding around the map
    '''
    # get vectors with lat/lon values
    lonVec, latVec = self.get_border()
    lonVec = np.array(lonVec)
    latVec = np.array(latVec)

    # estimate mean/min/max values of lat/lon of the shown area
    # (real lat min max +/- latBorder) and (real lon min max +/- lonBorder)
    minLon = max(-180, lonVec.min() - lonBorder)
    maxLon = min(180, lonVec.max() + lonBorder)
    minLat = max(-90, latVec.min() - latBorder)
    maxLat = min(90, latVec.max() + latBorder)
    meanLon = lonVec.mean()
    meanLat = latVec.mean()

    # generate template map (can be also tmerc)
    f = plt.figure(num=1, figsize=figureSize, dpi=dpi)
    bmap = Basemap(projection=projection,
                   lat_0=meanLat, lon_0=meanLon,
                   llcrnrlon=minLon, llcrnrlat=minLat,
                   urcrnrlon=maxLon, urcrnrlat=maxLat,
                   resolution=resolution)

    # add content: coastline, continents, meridians, parallels
    bmap.drawcoastlines()
    bmap.fillcontinents(color=continetsColor)
    bmap.drawmeridians(np.linspace(minLon, maxLon, meridians))
    bmap.drawparallels(np.linspace(minLat, maxLat, parallels))

    # convert lat/lons to map units
    mapX, mapY = bmap(list(lonVec.flat), list(latVec.flat))

    # from x/y vectors create a Patch to be added to map
    boundary = Polygon(zip(mapX, mapY), alpha=pAlpha, ec=pLine, fc=pColor)

    # add patch to the map
    plt.gca().add_patch(boundary)
    plt.gca().set_aspect('auto')

    # save figure and close
    plt.savefig(outputFileName, bbox_inches='tight',
                                dpi=dpi,
                                pad_inches=padding)
    plt.close('all')
